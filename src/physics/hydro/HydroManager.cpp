#include "HydroManager.hpp"

using namespace std;
using namespace chrono;


HydroManager::HydroManager(shared_ptr<Loader> ldr,
                           shared_ptr<GridManager> gridMnr,
                           shared_ptr<Pusher> pshr):loader(move(ldr)),
                           gridMgr(move(gridMnr)), pusher(move(pshr)){
                                
    logger.reset(new Logger());
    initialize();
    logger->writeMsg("[HydroManager] create...OK", DEBUG);
}

void HydroManager::initialize(){
    gatherMoments(PREDICTOR);
    gatherMoments(CORRECTOR);
}


void HydroManager::calculateAvgFluidVelocity4AllSpiecies(int phase){
   
    int numOfSpecies = loader->getNumberOfSpecies();
    
    Particle** particles = pusher->getParticles();
    
    int totalPrtclNumber = pusher->getTotalParticleNumber();
    
    string msg1 ="[HydroManager] totalPrtclNumber = "+to_string(totalPrtclNumber);
    logger->writeMsg(msg1.c_str(), DEBUG);
    
    
    int posShift = 0, velShift = 0;
    if( phase == CORRECTOR ){
        posShift = 3;
        velShift = 3;
    }
    
    double x, y, z;
    double lx, ly, lz;
    int idx, idxG2, idx_x, idx_y, idx_z, i, j, k, coord, spn;
    
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    int xSize = loader->resolution[0];
    int ySize = loader->resolution[1];
    int zSize = loader->resolution[2];
    
    int xSizeG2 = xSize+2;
    int ySizeG2 = ySize+2;
    int zSizeG2 = zSize+2;
    
    int G2nodesNumber = xSizeG2*ySizeG2*zSizeG2;
    double alpha, betta, gamma, alpha0, betta0, gamma0, weight;
    
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    
    double* weightsCharge = new double[G2nodesNumber*numOfSpecies];
    double* weights       = new double[G2nodesNumber*numOfSpecies];
    double* velocityWeighted = new double[G2nodesNumber*numOfSpecies*3];
    
    for( idx = 0; idx < G2nodesNumber; idx++ ){
        for( spn = 0; spn < numOfSpecies; spn++ ){
            weightsCharge[numOfSpecies*idx+spn] = 0.0;
            weights[numOfSpecies*idx+spn] = 0.0;
            for( coord = 0; coord < 3; coord++){
                velocityWeighted[(numOfSpecies*idx+spn)*3+coord] = 0.0;
            }
        }
    }


    double G2shift = 0.5;// in pixels
    
    
    //
    //
    //             FIRST ORDER SHAPE FUNCTION
    //
    //         |----x-------------0---------x-----| G2
    //              i <---------> xp <----> i+1
    //                     D           1-D
    //  i = int(xp)
    //  D = xp - i
    //  for i   weight = 1 - D
    //  for i+1 weight =     D
    
    
    double neighbourhood[8][3] = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0},
                                  {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};
    
    double* pos;
    double* vel;
    int type;
    double charge;
    for( idx = 0; idx < totalPrtclNumber; idx++ ){        
        pos  = particles[idx]->getPosition();
        vel  = particles[idx]->getVelocity();
        type = particles[idx]->getType();
        charge = particles[idx]->getCharge();
        
        if ( charge < MIN_CHARGE ){
            charge = MIN_CHARGE;
        }
        x = 0.5*(pos[0]+pos[3]);
        y = 0.5*(pos[1]+pos[4]);
        z = 0.5*(pos[2]+pos[5]);
        
        x = (x - domainShiftX)/dx+G2shift;
        y = (y - domainShiftY)/dy+G2shift;
        z = (z - domainShiftZ)/dz+G2shift;
        
        alpha0 = x-int(x);
        betta0 = y-int(y);
        gamma0 = z-int(z);
        
        double alphas[8] = {1.0-alpha0, alpha0    , 1.0-alpha0, alpha0,
                            1.0-alpha0, alpha0    , 1.0-alpha0, alpha0};
        double bettas[8] = {1.0-betta0, 1.0-betta0, betta0    , betta0,
                            1.0-betta0, 1.0-betta0, betta0    , betta0};
        double gammas[8] = {1.0-gamma0, 1.0-gamma0, 1.0-gamma0, 1.0-gamma0,
                            gamma0    , gamma0    , gamma0    , gamma0};
        
        x = pos[0+posShift];
        y = pos[1+posShift];
        z = pos[2+posShift];
        
        i = int((x - domainShiftX)/dx+G2shift);// G2 index
        j = int((y - domainShiftY)/dy+G2shift);
        k = int((z - domainShiftZ)/dz+G2shift);
        
        idxG2 = IDX(i ,j ,k, xSizeG2, ySizeG2, zSizeG2);
        
        for( int neigh_num = 0; neigh_num < 8; neigh_num++ ){
             
             idx_x = i + neighbourhood[neigh_num][0];
             idx_y = j + neighbourhood[neigh_num][1];
             idx_z = k + neighbourhood[neigh_num][2];
             
             idxG2 = IDX(idx_x, idx_y, idx_z, xSizeG2, ySizeG2, zSizeG2);
             
             alpha = alphas[neigh_num];
             betta = bettas[neigh_num];
             gamma = gammas[neigh_num];
             
             weight = alpha*betta*gamma;
             
             weightsCharge[numOfSpecies*idxG2+type] += weight*charge;
             weights[numOfSpecies*idxG2+type]       += weight;
             for( coord = 0; coord < 3; coord++ ){
                 velocityWeighted[(numOfSpecies*idxG2+type)*3+coord]
                 += weight*vel[coord+velShift];
            }
        }
    }
    
    for( spn = 0; spn < numOfSpecies; spn++ ){
        for( idx = 0; idx < G2nodesNumber; idx++ ){            
            gridMgr->setVectorVariableForNodeG2(idx, gridMgr->DENS_AUX(spn), WEIGHTS,
                                                weights[numOfSpecies*idx+spn]);
            gridMgr->setVectorVariableForNodeG2(idx, gridMgr->DENS_AUX(spn), WEIGHTS_Q,
                                                weightsCharge[numOfSpecies*idx+spn]);
            for( coord = 0; coord < 3; coord++ ){
                gridMgr->setVectorVariableForNodeG2(idx, gridMgr->DENS_VEL(spn), 1+coord,
                                                    velocityWeighted[(numOfSpecies*idx+spn)*3+coord]);
            }
        }
        
        gridMgr->gatherBoundaryUsingNeighbor(gridMgr->DENS_VEL(spn));
        gridMgr->applyBC(gridMgr->DENS_VEL(spn));
        gridMgr->gatherBoundaryUsingNeighbor(gridMgr->DENS_AUX(spn));
        gridMgr->applyBC(gridMgr->DENS_AUX(spn));
    }
   
    delete [] weightsCharge;
    delete [] weights;
    delete [] velocityWeighted;
}


void HydroManager::gatherMoments(int phase){
    
    auto start_time = high_resolution_clock::now();
    
    string msg0 ="[HydroManager] start to gather moments ";
    logger->writeMsg(msg0.c_str(), DEBUG);
    
    calculateAvgFluidVelocity4AllSpiecies(phase);
    
    auto end_time1 = high_resolution_clock::now();
    string msg1 ="[HydroManager] ------ calculateAvgFluidVelocity4AllSpiecies() ------ duration = "
                    +to_string(duration_cast<milliseconds>(end_time1 - start_time).count())+" ms";
    logger->writeMsg(msg1.c_str(), DEBUG);
    
    int xSize = loader->resolution[0];
    int ySize = loader->resolution[1];
    int zSize = loader->resolution[2];
    
    int G2nodesNumber = (xSize+2)*(ySize+2)*(zSize+2);
    int idx, coord, spn;
    
    double particleSpatialWeigth, particleFullWeigthQ, pw;
    
    int numOfSpecies = loader->getNumberOfSpecies();
    
    map<int, VectorVar**> dens_vel, dens_aux;
    
    for( spn = 0; spn < numOfSpecies; spn++ ){
            dens_vel[spn] = gridMgr->getVectorVariableOnG2(gridMgr->DENS_VEL(spn));
            dens_aux[spn] = gridMgr->getVectorVariableOnG2(gridMgr->DENS_AUX(spn));
    }
    
    double* densOfIons = new double[numOfSpecies*G2nodesNumber];
    double* densOfIonsNoCharge = new double[numOfSpecies*G2nodesNumber];
    double* massOfIonsNoCharge = new double[numOfSpecies*G2nodesNumber];
    
    double vel;
    
    for( idx = 0; idx < G2nodesNumber; idx++ ){
        for( spn = 0; spn < numOfSpecies; spn++ ){
            
            pw = pusher->getParticleWeight4Type(spn);
    
            particleSpatialWeigth = dens_aux[spn][idx]->getValue()[WEIGHTS];
            
            if( particleSpatialWeigth < EPS8 ){
                for( coord = 0; coord < 3; coord++ ){
                    vel = 0.0;
                    gridMgr->setVectorVariableForNodeG2(idx, gridMgr->DENS_VEL(spn), 1+coord, vel);
                }
            }else{
                for( coord = 0; coord < 3; coord++ ){
                    vel = dens_vel[spn][idx]->getValue()[1+coord]/particleSpatialWeigth;
                    gridMgr->setVectorVariableForNodeG2(idx, gridMgr->DENS_VEL(spn), 1+coord, vel);
                }
            }
            
            // multiply spatial weight by elementary particle weight for a certain sort of spicies
            particleFullWeigthQ = dens_aux[spn][idx]->getValue()[WEIGHTS_Q]*pw;
            
            gridMgr->setVectorVariableForNodeG2(idx, gridMgr->DENS_VEL(spn), 0, particleFullWeigthQ);

            gridMgr->setVectorVariableForNodeG2(idx, gridMgr->DENS_AUX(spn), WEIGHTS, dens_aux[spn][idx]->getValue()[WEIGHTS]*pw);

            densOfIons[numOfSpecies*idx+spn] = 0.0;
            densOfIonsNoCharge[numOfSpecies*idx+spn] = 0.0;
            massOfIonsNoCharge[numOfSpecies*idx+spn] = 0.0;
        }
    }
    
    double mass = 0.0;
    for( spn = 0; spn < numOfSpecies; spn++ ){
        mass = pusher->getParticleMass4Type(spn);
        for( idx = 0; idx < G2nodesNumber; idx++ ){
            densOfIons[numOfSpecies*idx+spn]         = dens_vel[spn][idx]->getValue()[0];
            densOfIonsNoCharge[numOfSpecies*idx+spn] = dens_aux[spn][idx]->getValue()[WEIGHTS];
            massOfIonsNoCharge[numOfSpecies*idx+spn] = dens_aux[spn][idx]->getValue()[WEIGHTS]*mass;
        }
    }

    double* totIonDensNoCharge = new double[G2nodesNumber];
    double* totIonMassNoCharge = new double[G2nodesNumber];
    double* densEle  = new double[G2nodesNumber];
    double* fluidVel = new double[3*G2nodesNumber];
    
    for( idx = 0; idx < G2nodesNumber; idx++ ){
        densEle[idx] = 0.0;
        totIonDensNoCharge[idx] = 0.0;
        totIonMassNoCharge[idx] = 0.0;
        for( coord = 0; coord < 3; coord++ ){
            fluidVel[3*idx+coord] = 0.0;
        }
    }

    for( idx = 0; idx < G2nodesNumber; idx++ ){
        for( spn = 0; spn < numOfSpecies; spn++ ){
            densEle[idx]            += densOfIons[numOfSpecies*idx+spn];
            totIonDensNoCharge[idx] += densOfIonsNoCharge[numOfSpecies*idx+spn];
            totIonMassNoCharge[idx] += massOfIonsNoCharge[numOfSpecies*idx+spn];
        }
        gridMgr->setVectorVariableForNodeG2(idx, DENSELEC, 0, densEle[idx]);
        gridMgr->setVectorVariableForNodeG2(idx, DENSELEC, 1, totIonDensNoCharge[idx]);//cell dens no charge
        if ( totIonDensNoCharge[idx] > EPS8 ){
            gridMgr->setVectorVariableForNodeG2(idx, DENSELEC, 2, densEle[idx]/totIonDensNoCharge[idx]);//cell charge
        }else{
            gridMgr->setVectorVariableForNodeG2(idx, DENSELEC, 2, 0.0);
        }
        if ( densEle[idx] > EPS8 ){
            gridMgr->setVectorVariableForNodeG2(idx, DENSELEC, 3, totIonMassNoCharge[idx]/densEle[idx]);//cell mass
        }else{
            gridMgr->setVectorVariableForNodeG2(idx, DENSELEC, 3, 0.0);
        }
    }
    
    for( idx = 0; idx < G2nodesNumber; idx++ ){
        for( spn = 0; spn < numOfSpecies; spn++ ){
            double densEleRevert = densEle[idx] < EPS8 ? 0.0 : edgeProfile(densEle[idx])/densEle[idx];
            double ionDens = densOfIons[numOfSpecies*idx+spn];
            for( coord = 0; coord < 3; coord++ ){
                vel = dens_vel[spn][idx]->getValue()[1+coord];
                fluidVel[3*idx+coord] += ionDens*vel*densEleRevert;
            }
        }
        for( coord = 0; coord < 3; coord++ ){
            gridMgr->setVectorVariableForNodeG2(idx, VELOCION, coord, fluidVel[3*idx+coord]);
        }
    }

    gridMgr->smoothDensAndIonVel();
    gridMgr->applyBC(DENSELEC);
    gridMgr->applyBC(VELOCION);
    
    delete [] densOfIons;
    delete [] densOfIonsNoCharge;
    delete [] massOfIonsNoCharge;
    delete [] densEle;
    delete [] totIonDensNoCharge;
    delete [] totIonMassNoCharge;
    delete [] fluidVel;
    
    auto end_time = high_resolution_clock::now();
    string msg ="[HydroManager] gatherMoments()... duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
    
}



void HydroManager::setIonPressureTensor(){

    auto start_time = high_resolution_clock::now();
    
    string msg0 ="[HydroManager] start to calculate ion pressure tensor ";
    logger->writeMsg(msg0.c_str(), DEBUG);
    
    int numOfSpecies = loader->getNumberOfSpecies();        
    double pxx, pxy, pxz, pyy, pyz, pzz;
    double vx, vy, vz;
    
    Particle** particles = pusher->getParticles();
    int totalPrtclNumber = pusher->getTotalParticleNumber();
    
    int posShift = 0, velShift = 0;
    
    double x, y, z;
    double lx, ly, lz;
    int idx, idxG2, idx_x, idx_y, idx_z, i, j, k, coord;
    
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    int xSize = loader->resolution[0];
    int ySize = loader->resolution[1];
    int zSize = loader->resolution[2];
    
    int xSizeG2 = xSize+2;
    int ySizeG2 = ySize+2;
    int zSizeG2 = zSize+2;
    
    int G2nodesNumber = xSizeG2*ySizeG2*zSizeG2;
    double alpha, betta, gamma, alpha0, betta0, gamma0, weight;
    
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    
    double G2shift = 0.5;// in pixels
    
    double neighbourhood[8][3] = {{0,0,0}, {1,0,0}, {0,1,0}, {1,1,0},
        {0,0,1}, {1,0,1}, {0,1,1}, {1,1,1}};
    
    double* pos;
    double* vel;
    int type;
    double pw, mass;
    for( idx = 0; idx < G2nodesNumber; idx++ ){
        for( type = 0; type < numOfSpecies; type++ ){
            for( int comp = 0; comp < 6; comp++ ){
                gridMgr->setVectorVariableForNodeG2(idx, gridMgr->ION_PRESSURE(type), comp, 0.0);
            }
        }
    }
    map<int, VectorVar**> dens_vel;
    for( type = 0; type < numOfSpecies; type++ ){
            dens_vel[type] = gridMgr->getVectorVariableOnG2(gridMgr->DENS_VEL(type));
    }

    for( idx=0; idx < totalPrtclNumber; idx++ ){
        
        pos  = particles[idx]->getPosition();
        vel  = particles[idx]->getVelocity();
        type = particles[idx]->getType();
        pw   = pusher->getParticleWeight4Type(type);
        mass = pusher->getParticleMass4Type(type);
        
        x = 0.5*(pos[0]+pos[3]);
        y = 0.5*(pos[1]+pos[4]);
        z = 0.5*(pos[2]+pos[5]);
        
        x = (x - domainShiftX)/dx+G2shift;
        y = (y - domainShiftY)/dy+G2shift;
        z = (z - domainShiftZ)/dz+G2shift;
        
        alpha0 = x-int(x);
        betta0 = y-int(y);
        gamma0 = z-int(z);
        
        double alphas[8] = {1.0-alpha0, alpha0    , 1.0-alpha0, alpha0,
            1.0-alpha0, alpha0    , 1.0-alpha0, alpha0};
        double bettas[8] = {1.0-betta0, 1.0-betta0, betta0    , betta0,
            1.0-betta0, 1.0-betta0, betta0    , betta0};
        double gammas[8] = {1.0-gamma0, 1.0-gamma0, 1.0-gamma0, 1.0-gamma0,
            gamma0    , gamma0    , gamma0    , gamma0};
        
        x = pos[0+posShift];
        y = pos[1+posShift];
        z = pos[2+posShift];
        
        i = int((x - domainShiftX)/dx+G2shift);// G2 index
        j = int((y - domainShiftY)/dy+G2shift);
        k = int((z - domainShiftZ)/dz+G2shift);
        
        idxG2 = IDX(i ,j ,k, xSizeG2, ySizeG2, zSizeG2);
       

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

            pxx = pw*mass*weight*(vel[0] - vx)*(vel[0] - vx);
            pxy = pw*mass*weight*(vel[0] - vx)*(vel[1] - vy);
            pxz = pw*mass*weight*(vel[0] - vx)*(vel[2] - vz);
            pyy = pw*mass*weight*(vel[1] - vy)*(vel[1] - vy);
            pyz = pw*mass*weight*(vel[1] - vy)*(vel[2] - vz);
            pzz = pw*mass*weight*(vel[2] - vz)*(vel[2] - vz);


            gridMgr->addVectorVariableForNodeG2(idxG2, gridMgr->ION_PRESSURE(type), 0, pxx);
            gridMgr->addVectorVariableForNodeG2(idxG2, gridMgr->ION_PRESSURE(type), 1, pxy);
            gridMgr->addVectorVariableForNodeG2(idxG2, gridMgr->ION_PRESSURE(type), 2, pxz);
            gridMgr->addVectorVariableForNodeG2(idxG2, gridMgr->ION_PRESSURE(type), 3, pyy);
            gridMgr->addVectorVariableForNodeG2(idxG2, gridMgr->ION_PRESSURE(type), 4, pyz);
            gridMgr->addVectorVariableForNodeG2(idxG2, gridMgr->ION_PRESSURE(type), 5, pzz);
        }
    }
    
    for( type = 0; type < numOfSpecies; type++ ){
        gridMgr->gatherBoundaryUsingNeighbor(gridMgr->ION_PRESSURE(type));
        gridMgr->applyBC(gridMgr->ION_PRESSURE(type));
    }

    auto end_time = high_resolution_clock::now();
    string msg ="[HydroManager] setIonPressureTensor()... duration = "
    +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}
