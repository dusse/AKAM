#include "IonIonCollisionManager.hpp"

using namespace std;
using namespace chrono;

IonIonCollisionManager::IonIonCollisionManager(shared_ptr<Loader> ldr,
                                               shared_ptr<GridManager> gridMnr,
                                               shared_ptr<Pusher> pshr):
                                               loader(move(ldr)), 
                                               gridMgr(move(gridMnr)),
                                               pusher(move(pshr))
{


    logger.reset(new Logger());
    initialize();
    logger->writeMsg("[IonIonCollisionManager] create...OK", DEBUG);
}

void IonIonCollisionManager::initialize(){
}

IonIonCollisionManager::~IonIonCollisionManager(){
}

void IonIonCollisionManager::collideIons(int phase){
    auto start_time = high_resolution_clock::now();
    string msg0 ="[IonIonCollisionManager] start to collide ions ";
    logger->writeMsg(msg0.c_str(), DEBUG);

    int velShift = 0;
    if( phase == CORRECTOR ){
        velShift = 3;
    }

    int numOfSpecies = loader->getNumberOfSpecies();
    Particle** particles = pusher->getParticles();
    int totalPrtclNumber = pusher->getTotalParticleNumber();

    int type, type2, idx, idxG2;
    
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
            particlesNumber[numOfSpecies*idx+type] = 0;
            particlesOfTheGivenType[type] = particleIndecies;
        }
        particlesInEachCell[idx] = particlesOfTheGivenType;
    }
    
    double* pos;
    double* vel;
    double pw, mass;
    map<int, VectorVar**> dens_vel;
    for( type = 0; type < numOfSpecies; type++ ){
            dens_vel[type] = gridMgr->getVectorVariableOnG2(gridMgr->DENS_VEL(type));
    }

    int ptclIdx, ion1idx, ion2idx;
    double dens1, dens2;
    double pw1, pw2;
    int group1Idx, group2Idx, restIdx;
    for( idxG2 = 0; idxG2 < G2nodesNumber; idxG2++ ){

        /*** shuffle particles ***/
        for( type = 0; type < numOfSpecies; type++ ){
            random_shuffle(particlesInEachCell[idxG2][type].begin(), particlesInEachCell[idxG2][type].end());
        }
        /** based on work of Nicolas Loic (2017)
         * Effects of collisions on the magnetic streaming instability 
         * (Doctoral dissertation, Paris 6).**/        
        for( type = 0; type < numOfSpecies; type++ ){
            if( pusher->getIfParticleTypeIsFrozen(type) == 1 ) continue;

            dens1 = dens_vel[type][idxG2]->getValue()[0];

            int numOfPartclsOfGvnType = particlesNumber[numOfSpecies*idxG2+type];
            if( numOfPartclsOfGvnType == 0 ) continue;
            /** intra-species collisions
                Takizuka, T. and H. Abe (1977). “A binary collision model for plasma simulation with a
                particle code” Journal of Computational Physics 25 **/
            if( numOfPartclsOfGvnType%2 == 0 ){
                for( ptclIdx = 0; ptclIdx < numOfPartclsOfGvnType/2; ptclIdx++ ){
                    ion1idx = particlesInEachCell[idxG2][type][2*ptclIdx];  
                    ion2idx = particlesInEachCell[idxG2][type][2*ptclIdx+1];
                    scatterVelocities(velShift, ion1idx, ion2idx, type, type, dens1, dens1, dens1,
                                      particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                      1);
                }
            }else{
                for( ptclIdx = 0; ptclIdx < (numOfPartclsOfGvnType/2)-1; ptclIdx++ ){
                    ion1idx = particlesInEachCell[idxG2][type][2*ptclIdx];
                    ion2idx = particlesInEachCell[idxG2][type][2*ptclIdx+1];
                    scatterVelocities(velShift, ion1idx, ion2idx, type, type, dens1, dens1, dens1,
                                      particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                      1);
                }
                if( numOfPartclsOfGvnType >= 3 ){
                    ion1idx = particlesInEachCell[idxG2][type][numOfPartclsOfGvnType-2];
                    ion2idx = particlesInEachCell[idxG2][type][numOfPartclsOfGvnType-1];
                    scatterVelocities(velShift, ion1idx, ion2idx, type, type, dens1, dens1, dens1,
                                      particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                      0.5);

                    ion1idx = particlesInEachCell[idxG2][type][numOfPartclsOfGvnType-3];
                    ion2idx = particlesInEachCell[idxG2][type][numOfPartclsOfGvnType-1];
                    scatterVelocities(velShift, ion1idx, ion2idx, type, type, dens1, dens1, dens1,
                                      particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                      0.5);

                    ion1idx = particlesInEachCell[idxG2][type][numOfPartclsOfGvnType-3];
                    ion2idx = particlesInEachCell[idxG2][type][numOfPartclsOfGvnType-2];
                    scatterVelocities(velShift, ion1idx, ion2idx, type, type, dens1, dens1, dens1,
                                      particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                      0.5);
                }
            }
            
            pw1 = pusher->getParticleWeight4Type(type);
            /** inter-species collisions
                Miller, R. H. and M. R. Combi (1994). “A Coulomb collision algorithm for weighted particle
                simulations” Geophysical Research Letters 21 **/
            for( type2 = type+1; type2 < numOfSpecies; type2++ ){
                if( pusher->getIfParticleTypeIsFrozen(type2) == 1 ) continue;

                int numOfPartclsOfGvnType2 = particlesNumber[numOfSpecies*idxG2+type2];
                if( numOfPartclsOfGvnType2 == 0 ) continue;

                pw2 = pusher->getParticleWeight4Type(type2);
                dens2 = dens_vel[type2][idxG2]->getValue()[0];
                
                if( numOfPartclsOfGvnType >= numOfPartclsOfGvnType2 ){

                    double weightFactor = (pw1 >= pw2) ? pw1/pw2 : 1;

                    int    quotient  = floor((double)numOfPartclsOfGvnType/(double)numOfPartclsOfGvnType2);
                    double remainder = ((double)numOfPartclsOfGvnType/(double)numOfPartclsOfGvnType2)-quotient;

                    int firstGroupSpecie2 = (int)round(remainder*double(numOfPartclsOfGvnType2));

                    for( group1Idx = 0; group1Idx < firstGroupSpecie2; group1Idx++ ){
                        ion1idx = particlesInEachCell[idxG2][type2][group1Idx];
                        for (restIdx = 0; restIdx < quotient+1; restIdx++) {
                            ptclIdx = group1Idx*(quotient+1)+restIdx;
                            ion2idx = particlesInEachCell[idxG2][type][ptclIdx];
                            scatterVelocities(velShift, ion1idx, ion2idx, type2, type, dens2, dens1,
                                         min(dens1,dens2), particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                         weightFactor);
                        }
                    }

                    for( group2Idx = firstGroupSpecie2; group2Idx < numOfPartclsOfGvnType2; group2Idx++ ){
                        ion1idx = particlesInEachCell[idxG2][type2][group2Idx];
                        for( restIdx = 0; restIdx < quotient; restIdx++ ){
                            ptclIdx = firstGroupSpecie2*(quotient+1)+(group2Idx-firstGroupSpecie2)*quotient+restIdx;
                            ion2idx = particlesInEachCell[idxG2][type][ptclIdx];
                            scatterVelocities(velShift, ion1idx, ion2idx, type2, type, dens2, dens1,
                                         min(dens1,dens2), particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                         weightFactor);
                        }
                    }
                }else{
                    double weightFactor = (pw1 >= pw2) ? 1 : pw2/pw1;

                    int    quotient  = floor((double)numOfPartclsOfGvnType2/(double)numOfPartclsOfGvnType);
                    double remainder = ((double)numOfPartclsOfGvnType2/(double)numOfPartclsOfGvnType)-quotient;
                    int firstGroupSpecie1 = (int)round(remainder*numOfPartclsOfGvnType);

                    for( group1Idx = 0; group1Idx < firstGroupSpecie1; group1Idx++ ){
                        ion1idx = particlesInEachCell[idxG2][type][group1Idx];
                        for( restIdx = 0; restIdx < quotient+1; restIdx++ ){
                            ptclIdx = group1Idx*(quotient+1)+restIdx;
                            ion2idx = particlesInEachCell[idxG2][type2][ptclIdx];
                            scatterVelocities(velShift, ion1idx, ion2idx, type, type2, dens1, dens2,
                                         min(dens1,dens2), particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                         weightFactor);
                        }
                    }

                    for( group2Idx = firstGroupSpecie1; group2Idx < numOfPartclsOfGvnType; group2Idx++ ){
                        ion1idx = particlesInEachCell[idxG2][type][group2Idx];
                        for( restIdx = 0; restIdx < quotient; restIdx++ ){
                            ptclIdx = firstGroupSpecie1*(quotient+1)+(group2Idx-firstGroupSpecie1)*quotient+restIdx;
                            ion2idx = particlesInEachCell[idxG2][type2][ptclIdx];
                            scatterVelocities(velShift, ion1idx, ion2idx, type, type2, dens1, dens2,
                                         min(dens1,dens2), particles[ion1idx]->getVelocity(), particles[ion2idx]->getVelocity(),
                                         weightFactor);
                        }
                    }
                }
            }
        }
    }


    delete [] particlesNumber;

    auto end_time = high_resolution_clock::now();
    auto msg ="[IonIonCollisionManager] collideIons()... duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}


void IonIonCollisionManager::scatterVelocities(
    int velShift,
    int ion1idx,
    int ion2idx,
    int type1,
    int type2,
    double dens1,
    double dens2,
    double lowestDensity,
    double* ion1Vel,
    double* ion2Vel,
    double factor){

    const double q1 = pusher->getParticleCharge4Type(type1);
    const double q2 = pusher->getParticleCharge4Type(type2);
    const double m1 = pusher->getParticleMass4Type(type1);
    const double m2 = pusher->getParticleMass4Type(type2);
    const double mu = m1 * m2 / (m1 + m2);
    const double w1 = pusher->getParticleWeight4Type(type1);
    const double w2 = pusher->getParticleWeight4Type(type2);

    double v1[3], v2[3];
    for (int i = 0; i < 3; ++i) {
        v1[i] = ion1Vel[i + velShift];
        v2[i] = ion2Vel[i + velShift];
    }

    double g[3] = {
        v1[0] - v2[0],
        v1[1] - v2[1],
        v1[2] - v2[2]
    };

    double gmag = sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
    if (gmag < EPSILON) return;

    const double defaultFactor =
        factor * loader->getCollisionFrequencyFactor() * loader->getTimeStep();

    const double coulombLog = loader->getDefaultCoulombLogarithm();

    double nu =
        (defaultFactor * q1*q1 * q2*q2 * lowestDensity * coulombLog) /
        (mu*mu * gmag*gmag*gmag);

    double r1 = RNM;
    r1 = -2.0 * log((r1 > EPSILON) ? r1 : EPSILON);
    double r2 = 2.0 * PI * RNM;

    double sinTheta, cosTheta;

    if (nu < 1.0) {
        double s = sqrt(nu * r1) * cos(r2);
        cosTheta = (1.0 - s*s) / (1.0 + s*s);
        sinTheta = 2.0 * s / (1.0 + s*s);
    } else {
        double theta = acos(1.0 - 2.0 * RNM);
        cosTheta = cos(theta);
        sinTheta = sin(theta);
    }

    double phi = 2.0 * PI * RNM;

    double e1[3] = { g[0]/gmag, g[1]/gmag, g[2]/gmag };

    double e2[3];
    if (fabs(e1[0]) < 0.9) {
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

    double gnew[3];
    for (int i = 0; i < 3; ++i) {
        gnew[i] = gmag * (
            cosTheta * e1[i]
          + sinTheta * (cos(phi) * e2[i] + sin(phi) * e3[i])
        );
    }

    double dg[3] = {
        gnew[0] - g[0],
        gnew[1] - g[1],
        gnew[2] - g[2]
    };

    double rnd = RNM;
    double alpha = (rnd < w2 / w1) ? 1.0 : 0.0;
    double betta = (rnd < w1 / w2) ? 1.0 : 0.0;

    for (int i = 0; i < 3; ++i) {
        pusher->setParticleVelocity(
            ion1idx,
            i + velShift,
            v1[i] + alpha * (mu / m1) * dg[i]
        );

        pusher->setParticleVelocity(
            ion2idx,
            i + velShift,
            v2[i] - betta * (mu / m2) * dg[i]
        );
    }
}
