#include "ElectronIonCollisionManager.hpp"

using namespace std;
using namespace chrono;

ElectronIonCollisionManager::ElectronIonCollisionManager(shared_ptr<Loader> ldr,
                                               shared_ptr<GridManager> gridMnr,
                                               shared_ptr<Pusher> pshr,
                                               shared_ptr<ClosureManager> cm):
                                               loader(move(ldr)), 
                                               gridMgr(move(gridMnr)),
                                               pusher(move(pshr)),
                                               closureMng(move(cm))
{
    logger.reset(new Logger());
    initialize();
    logger->writeMsg("[ElectronIonCollisionManager] create...OK", DEBUG);
}

void ElectronIonCollisionManager::initialize(){
    int xSize = loader->resolution[0];
    int ySize = loader->resolution[1];
    int zSize = loader->resolution[2];
    
    int xSizeG2 = xSize+2;
    int ySizeG2 = ySize+2;
    int zSizeG2 = zSize+2;
    
    int G2nodesNumber = xSizeG2*ySizeG2*zSizeG2;
    int numOfSpecies = loader->getNumberOfSpecies();

    ionTrPbeforeCollision = new double[G2nodesNumber*numOfSpecies];
}

ElectronIonCollisionManager::~ElectronIonCollisionManager(){
    delete [] ionTrPbeforeCollision;
}

void ElectronIonCollisionManager::collideElectronWithIons(int phase){
    auto start_time = high_resolution_clock::now();
    string msg0 ="[ElectronIonCollisionManager] start to collide electrons with ions ";
    logger->writeMsg(msg0.c_str(), DEBUG);

    int velShift = 0;
    if( phase == CORRECTOR ){
        velShift = 3;
    }

    int numOfSpecies = loader->getNumberOfSpecies();
    Particle** particles = pusher->getParticles();
    int totalPrtclNumber = pusher->getTotalParticleNumber();

    VectorVar** epres = gridMgr->getVectorVariableOnG2(PRESSURE);
    VectorVar** edens = gridMgr->getVectorVariableOnG2(DENSELEC);

    int type, type2, idx, idxG2, i, j, k ;
   
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];

    double G2shift = 0.5;// in pixels

    int xSize = loader->resolution[0];
    int ySize = loader->resolution[1];
    int zSize = loader->resolution[2];
    
    int xSizeG2 = xSize+2;
    int ySizeG2 = ySize+2;
    int zSizeG2 = zSize+2;
    
    int G2nodesNumber = xSizeG2*ySizeG2*zSizeG2;
    int* particlesNumber = new int[G2nodesNumber*numOfSpecies];
    map<int, map<int, vector<int>>> particlesInEachCell;

    for( idx = 0; idx < G2nodesNumber; idx++ ){
        map<int, vector<int>> particlesOfTheGivenType;
        vector<int> particleIndecies;
        for( type = 0; type < numOfSpecies; type++ ){
            particlesNumber[numOfSpecies*idx+type] = 0.0;
            particlesOfTheGivenType[type] = particleIndecies;
            ionTrPbeforeCollision[numOfSpecies*idx+type] = 0.0;
        }
        particlesInEachCell[idx] = particlesOfTheGivenType;
    }
  
    double neighbourhood[8][3] = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0},
                                  {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};

    double x, y, z;
    double alpha0, betta0, gamma0;
    double alpha, betta, gamma;
    double* pos;
    double* vel;
    double pw, mass;
    double vx, vy, vz, weight;
    double trPe, pxx, pyy, pzz;
    int idx_x, idx_y, idx_z;
    map<int, VectorVar**> dens_vel;
    for( type = 0; type < numOfSpecies; type++ ){
        dens_vel[type] = gridMgr->getVectorVariableOnG2(gridMgr->DENS_VEL(type));
    }


    for( idx = 0; idx < totalPrtclNumber; idx++ ){        
        pos  = particles[idx]->getPosition();
        vel  = particles[idx]->getVelocity();
        type = particles[idx]->getType();
        pw   = pusher->getParticleWeight4Type(type);
        mass = pusher->getParticleMass4Type(type);
        if( pusher->getIfParticleTypeIsFrozen(type) == 1 ) continue;
        
        x = pos[0+velShift];
        y = pos[1+velShift];
        z = pos[2+velShift];

        x = (x - domainShiftX)/dx+G2shift;
        y = (y - domainShiftY)/dy+G2shift;
        z = (z - domainShiftZ)/dz+G2shift;
        
        alpha0 = x-int(x); betta0 = y-int(y); gamma0 = z-int(z);
        
        double alphas[8] = {1.0-alpha0, alpha0,     1.0-alpha0, alpha0,
                            1.0-alpha0, alpha0,     1.0-alpha0, alpha0};
        double bettas[8] = {1.0-betta0, 1.0-betta0, betta0,     betta0,
                            1.0-betta0, 1.0-betta0, betta0,     betta0};
        double gammas[8] = {1.0-gamma0, 1.0-gamma0, 1.0-gamma0, 1.0-gamma0,
                                gamma0,     gamma0,     gamma0,     gamma0};
        
        x = pos[0+velShift];
        y = pos[1+velShift];
        z = pos[2+velShift];

        i = int((x - domainShiftX)/dx+G2shift);// G2 index
        j = int((y - domainShiftY)/dy+G2shift);
        k = int((z - domainShiftZ)/dz+G2shift);
        
        idxG2 = IDX(i, j, k, xSizeG2, ySizeG2, zSizeG2);
        particlesNumber[numOfSpecies*idxG2+type] += 1;
        particlesInEachCell[idxG2][type].push_back(idx);

        for( int neigh_num = 0; neigh_num < 8; neigh_num++ ){            
            idx_x = i + neighbourhood[neigh_num][0];
            idx_y = j + neighbourhood[neigh_num][1];
            idx_z = k + neighbourhood[neigh_num][2];
            idxG2 = IDX(idx_x ,idx_y ,idx_z, xSizeG2, ySizeG2, zSizeG2);
            
            vx = dens_vel[type][idxG2]->getValue()[1];
            vy = dens_vel[type][idxG2]->getValue()[2];
            vz = dens_vel[type][idxG2]->getValue()[3];
            
            alpha  = alphas[neigh_num];
            betta  = bettas[neigh_num];
            gamma  = gammas[neigh_num];
            weight = alpha*betta*gamma;

            pxx = pw*mass*weight*pow(vel[0+velShift] - vx, 2);
            pyy = pw*mass*weight*pow(vel[1+velShift] - vy, 2);
            pzz = pw*mass*weight*pow(vel[2+velShift] - vz, 2);
            ionTrPbeforeCollision[numOfSpecies*idxG2+type] += (pxx+pyy+pzz);
       	}
    }

    int ptclIdx, ionidx;
    double dens, etemp, itemp;
    
    double pw1, pw2;
    double ts = loader->getTimeStep();

    for( int i = 1; i < xSize+1; i++ ){
        for( int j = 1; j < ySize+1; j++ ){
            for( int k = 1; k < zSize+1; k++ ){
                int idxG2 = IDX(i, j, k, xSizeG2, ySizeG2, zSizeG2);
                double ne = edens[idxG2]->getValue()[0];
                double trPe  = epres[idxG2]->getValue()[0] +
                               epres[idxG2]->getValue()[3] +
                               epres[idxG2]->getValue()[5];
                double etemp = trPe / 3.0;                    
                etemp *= ne < EPS8 ? 0.0 : edgeProfile(ne)/ne;

                for( int type = 0; type < numOfSpecies; type++ ){
                    if( pusher->getIfParticleTypeIsFrozen(type) ) continue;
                    double massi = pusher->getParticleMass4Type(type);
                    int numOfPartclsOfGvnType = particlesNumber[numOfSpecies*idxG2 + type];
                    if( numOfPartclsOfGvnType == 0 ) continue;

                    double dens  = dens_vel[type][idxG2]->getValue()[0];
                    double itemp = ionTrPbeforeCollision[numOfSpecies*idxG2+type] / (3.0 * dens);

                    double deltaT = etemp - itemp;
                    if( fabs(deltaT) < EPS8 || trPe < 1e-3 ) continue;

                    double nuei = loader->eiCollFreq;
                    double mi = edens[idxG2]->getValue()[3];//cell ion mass
                    double me = mi*1e-1;
                    if( nuei == 0.0 ){
                        double z  = edens[idxG2]->getValue()[2];//cell charge
                        double logC = 10.0;
                        double denom = pow(mi*etemp+me*itemp, 1.5);
                        nuei =  sqrt(me*mi) * pow(z,2) * dens * logC / denom;
                    }
                    
                    double tau_ie = 1.0 / (nuei + EPS8);

                    double deltaE_max = me / mi * itemp;
                    double itemp_target = 0.0;
                    if( etemp < 1e-2 ){  //electrons are too cold to absorb full energy
                        itemp_target = max(itemp - deltaE_max, itemp * 0.9);
                    }else{
                        //relaxation equation: dTi/dt = (Te - Ti)/tau_ie
                        itemp_target = etemp - deltaT * exp(-ts / tau_ie);
                    }

                    double vth = sqrt(itemp_target / massi);
                    double itemp_current = itemp;
                    if( itemp_current == 0.0 ){
                        double r1, r2;
                        double vpb[3];
                        
                        for( int ptclIdx = 0; ptclIdx < numOfPartclsOfGvnType; ptclIdx++ ) {
                            int ionidx = particlesInEachCell[idxG2][type][ptclIdx];
                            if( particles[ionidx]->getCharge() < 1.0 ) continue;
                            double vx = particles[ionidx]->getVelocity()[0 + velShift];
                            double vy = particles[ionidx]->getVelocity()[1 + velShift];
                            double vz = particles[ionidx]->getVelocity()[2 + velShift];
                            double modv = sqrt(vx*vx+vy*vy+vz*vz);
                            if( modv == 0.0 ){
                                r1 = RNM;
                                r2 = RNM;
                                r1 = fabs(r1 - 1.0) < EPS8 ? r1 - EPS8 : r1;
                                r1 = r1 > EPS8 ? r1 : r1 + EPS8;
                                vpb[0] = sqrt(-2*log(r1))*vth * cos(2*PI*r2);
                                vpb[1] = sqrt(-2*log(r1))*vth * sin(2*PI*r2);

                                r1 = RNM;
                                r2 = RNM;
                                r1 = fabs(r1 - 1.0) < EPS8 ? r1 - EPS8 : r1;
                                r1 = r1 > EPS8 ? r1 : r1 + EPS8;
                                vpb[2] = sqrt(-2*log(r1))*vth * cos(2*PI*r2);
    
                                for( int c = 0; c < 3; c++ ){
                                    if( abs(vpb[c]) < 3.0*vth ){
                                        pusher->setParticleVelocity(ionidx, c + velShift, vpb[c]);
                                    }else{
                                        pusher->setParticleVelocity(ionidx, c + velShift, 3.0*vth*copysign(1.0, vpb[c]));
                                    }
                                }
                            }
                        }
                        continue;
                    }

                    double lambda = sqrt(itemp_target / itemp_current);
                    if( lambda > 3.0 ) lambda = 3.0;
                    if( lambda < 0.3 ) lambda = 0.3;

                    double v_mean[3] = {0.0, 0.0, 0.0};
                    
                    for( int ptclIdx = 0; ptclIdx < numOfPartclsOfGvnType; ptclIdx++ ){
                        int ionidx = particlesInEachCell[idxG2][type][ptclIdx];
                        for( int c = 0; c < 3; c++ )
                            v_mean[c] += particles[ionidx]->getVelocity()[c + velShift];
                    }
                    
                    for(int c = 0; c < 3; c++)
                        v_mean[c] /= numOfPartclsOfGvnType;

                    for( int ptclIdx = 0; ptclIdx < numOfPartclsOfGvnType; ptclIdx++ ){
                        int ionidx = particlesInEachCell[idxG2][type][ptclIdx];
                        for( int c = 0; c < 3; c++ ){
                            double v_old = particles[ionidx]->getVelocity()[c + velShift];
                            double v_th  = v_old - v_mean[c];
                            double v_new = v_mean[c] + lambda * v_th;  
                            pusher->setParticleVelocity(ionidx, c + velShift, v_new);
                        }
                    }

                    if( numOfPartclsOfGvnType > 1 ){                    
                        random_shuffle(particlesInEachCell[idxG2][type].begin(), 
                                       particlesInEachCell[idxG2][type].end());                    
                        double collisionFrequency = 0.1;
                        for( int ptclIdx = 0; ptclIdx < numOfPartclsOfGvnType - 1; ptclIdx += 2 ){                    
                            int ion1idx = particlesInEachCell[idxG2][type][ptclIdx];
                            int ion2idx = particlesInEachCell[idxG2][type][ptclIdx + 1];                            
                            scatterVelocities(velShift, ion1idx, ion2idx, type, type,
                                              particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(), 
                                              collisionFrequency);
                        }

                        if( numOfPartclsOfGvnType % 2 != 0 ){
                            int ion1idx = particlesInEachCell[idxG2][type][numOfPartclsOfGvnType - 1];
                            int ion2idx = particlesInEachCell[idxG2][type][RNM * (numOfPartclsOfGvnType - 1)];                    
                            scatterVelocities(velShift, ion1idx, ion2idx, type, type,
                                              particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                              collisionFrequency);
                        }
                    }
                }
            }
        }
    }

    delete [] particlesNumber;
    auto end_time = high_resolution_clock::now();
    auto msg ="[ElectronIonCollisionManager] collideElectronsWithIons()... duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}

void ElectronIonCollisionManager::applyElectronPressure(int phase){
    auto start_time = high_resolution_clock::now();
    string msg0 ="[ElectronIonCollisionManager] applyElectronPressure... ";
    logger->writeMsg(msg0.c_str(), DEBUG);
    VectorVar** edens = gridMgr->getVectorVariableOnG2(DENSELEC);
    int velShift = 0;
    if( phase == CORRECTOR ){
        velShift = 3;
    }

    int numOfSpecies = loader->getNumberOfSpecies();
    Particle** particles = pusher->getParticles();
    int totalPrtclNumber = pusher->getTotalParticleNumber();
    
    int type, type2, idx, idxG2, i, j, k ;
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];

    double G2shift = 0.5;// in pixels

    int xSize = loader->resolution[0];
    int ySize = loader->resolution[1];
    int zSize = loader->resolution[2];
    
    int xSizeG2 = xSize+2;
    int ySizeG2 = ySize+2;
    int zSizeG2 = zSize+2;
    
    int G2nodesNumber = xSizeG2*ySizeG2*zSizeG2;
    double* ionTrPafterCollision = new double[G2nodesNumber*numOfSpecies];

    for( idx = 0; idx < G2nodesNumber; idx++ ){
        closureMng->setAlphas(idx, 0, 0.0);
        closureMng->setAlphas(idx, 3, 0.0);
        closureMng->setAlphas(idx, 5, 0.0);

        for( type = 0; type < numOfSpecies; type++ ){
            ionTrPafterCollision[numOfSpecies*idx+type] = 0.0;
        }

    }
  
    double neighbourhood[8][3] = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0},
                                  {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};

    double x, y, z;
    double alpha0, betta0, gamma0;
    double alpha, betta, gamma;
    double* pos;
    double* vel;
    double pw, mass;
    double vx, vy, vz, weight;
    double pxx, pyy, pzz;
    int idx_x, idx_y, idx_z;
    map<int, VectorVar**> dens_vel;
    for( type = 0; type < numOfSpecies; type++ ){
        dens_vel[type] = gridMgr->getVectorVariableOnG2(gridMgr->DENS_VEL(type));
    } 

    for( idx = 0; idx < totalPrtclNumber; idx++ ){
        
        pos  = particles[idx]->getPosition();
        vel  = particles[idx]->getVelocity();
        type = particles[idx]->getType();
        pw   = pusher->getParticleWeight4Type(type);
        mass = pusher->getParticleMass4Type(type);
        
        x = pos[0+velShift];
        y = pos[1+velShift];
        z = pos[2+velShift];

        x = (x - domainShiftX)/dx+G2shift;
        y = (y - domainShiftY)/dy+G2shift;
        z = (z - domainShiftZ)/dz+G2shift;
        
        alpha0 = x-int(x); betta0 = y-int(y); gamma0 = z-int(z);
        
        double alphas[8] = {1.0-alpha0, alpha0,     1.0-alpha0, alpha0,
                            1.0-alpha0, alpha0,     1.0-alpha0, alpha0};
        double bettas[8] = {1.0-betta0, 1.0-betta0, betta0,     betta0,
                            1.0-betta0, 1.0-betta0, betta0,     betta0};
        double gammas[8] = {1.0-gamma0, 1.0-gamma0, 1.0-gamma0, 1.0-gamma0,
                                gamma0,     gamma0,     gamma0,     gamma0};
        
        x = pos[0+velShift];
        y = pos[1+velShift];
        z = pos[2+velShift];

        i = int((x - domainShiftX)/dx+G2shift);// G2 index
        j = int((y - domainShiftY)/dy+G2shift);
        k = int((z - domainShiftZ)/dz+G2shift);
        
        for( int neigh_num = 0; neigh_num < 8; neigh_num++ ){
            idx_x = i + neighbourhood[neigh_num][0];
            idx_y = j + neighbourhood[neigh_num][1];
            idx_z = k + neighbourhood[neigh_num][2];
            idxG2 = IDX(idx_x ,idx_y ,idx_z, xSizeG2, ySizeG2, zSizeG2);
            
            vx = dens_vel[type][idxG2]->getValue()[1];
            vy = dens_vel[type][idxG2]->getValue()[2];
            vz = dens_vel[type][idxG2]->getValue()[3];
            
            alpha = alphas[neigh_num];
            betta = bettas[neigh_num];
            gamma = gammas[neigh_num];
            weight = alpha*betta*gamma;

            pxx = pw*mass*weight*pow(vel[0+velShift] - vx, 2);
            pyy = pw*mass*weight*pow(vel[1+velShift] - vy, 2);
            pzz = pw*mass*weight*pow(vel[2+velShift] - vz, 2);
            ionTrPafterCollision[numOfSpecies*idxG2+type] += (pxx+pyy+pzz);
        }
    }

    double deltaP;
    VectorVar** epres = gridMgr->getVectorVariableOnG2(PRESSURE);
    double trPe  = epres[idxG2]->getValue()[0] +
                                   epres[idxG2]->getValue()[3] +
                                   epres[idxG2]->getValue()[5];
    
    for( int type = 0; type < numOfSpecies; type++ ){
        for( idx = 0; idx < G2nodesNumber; idx++ ){
            double pion_before = ionTrPbeforeCollision[numOfSpecies*idx+type];
            double pion_after  = ionTrPafterCollision[ numOfSpecies*idx+type];
            deltaP = pion_before - pion_after;
            if (deltaP > 0.0 || (trPe + deltaP) > EPS8) {
                double pressure2Pump = deltaP / 3.0;
                closureMng->addAlphas(idx, 0, pressure2Pump);
                closureMng->addAlphas(idx, 3, pressure2Pump);
                closureMng->addAlphas(idx, 5, pressure2Pump);
            }
        }
    }

    delete [] ionTrPafterCollision;
    auto end_time = high_resolution_clock::now();
    auto msg ="[ElectronIonCollisionManager] applyElectronPressure()... duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}

void ElectronIonCollisionManager::scatterVelocities(
    int velShift,
    int ion1idx,
    int ion2idx,
    int type1,
    int type2,
    double* ion1Vel,
    double* ion2Vel,
    double collisionFrequency){

    const double m1 = pusher->getParticleMass4Type(type1);
    const double m2 = pusher->getParticleMass4Type(type2);
    const double invM = 1.0 / (m1 + m2);
    double v1[3], v2[3];
    for (int i = 0; i < 3; ++i) {
        v1[i] = ion1Vel[i + velShift];
        v2[i] = ion2Vel[i + velShift];
    }
    double Vcm[3];
    for( int i = 0; i < 3; ++i )
        Vcm[i] = (m1 * v1[i] + m2 * v2[i]) * invM;
    double g[3];
    for( int i = 0; i < 3; ++i )
        g[i] = v1[i] - v2[i];

    double gmag = sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
    if( gmag < EPSILON ) return;

    double e1[3] = { g[0]/gmag, g[1]/gmag, g[2]/gmag };

    double e2[3];
    if( fabs(e1[0]) < 0.9 ){
        e2[0] = 0.0;
        e2[1] = -e1[2];
        e2[2] =  e1[1];
    } else {
        e2[0] = -e1[1];
        e2[1] =  e1[0];
        e2[2] = 0.0;
    }

    double norm = sqrt(e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]);
    e2[0] /= norm; e2[1] /= norm; e2[2] /= norm;

    double e3[3] = {
        e1[1]*e2[2] - e1[2]*e2[1],
        e1[2]*e2[0] - e1[0]*e2[2],
        e1[0]*e2[1] - e1[1]*e2[0]
    };

    double r1 = RNM;
    r1 = -2.0 * log((r1 > EPSILON) ? r1 : EPSILON);
    double r2 = 2.0 * PI * RNM;

    double s = sqrt(collisionFrequency * r1) * cos(r2);

    double cosTheta = (1.0 - s*s) / (1.0 + s*s);
    double sinTheta = 2.0 * s / (1.0 + s*s);
    double phi = 2.0 * PI * RNM;

    double gnew[3];
    for( int i = 0; i < 3; ++i ){
        gnew[i] = gmag * (
            cosTheta * e1[i]
          + sinTheta * (cos(phi) * e2[i] + sin(phi) * e3[i])
        );
    }
    for( int i = 0; i < 3; ++i ){
        double v1new = Vcm[i] + (m2 * invM) * gnew[i];
        double v2new = Vcm[i] - (m1 * invM) * gnew[i];

        pusher->setParticleVelocity(ion1idx, i + velShift, v1new);
        pusher->setParticleVelocity(ion2idx, i + velShift, v2new);
    }
}
