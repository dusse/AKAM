#include "EnvelopeEquationSolver.hpp"

using namespace std;
using namespace chrono;


EnvelopeEquationSolver::EnvelopeEquationSolver(shared_ptr<Loader> load,
                           shared_ptr<GridManager> grid,
                           shared_ptr<Pusher> pshr,
                           shared_ptr<ClosureManager> cm):loader(move(load)),
                            gridMgr(move(grid)), pusher(move(pshr)), closureMng(move(cm)){

    logger.reset(new Logger());    
    initialize();
    logger->writeMsg("[EnvelopeEquationSolver] create...OK", DEBUG);
}

void EnvelopeEquationSolver::initialize(){

    int xRes = loader->resolution[0],
        yRes = loader->resolution[1],
        zRes = loader->resolution[2];
    
    int xResG2 = xRes+2, yResG2 = yRes+2, zResG2 = zRes+2;
    int xResG4 = xRes+4, yResG4 = yRes+4, zResG4 = zRes+4;
    int totG2 = xResG2*yResG2*zResG2;
    int totG4 = xResG4*yResG4*zResG4;

    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];

    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];

    potentialN        = new complex<double>[totG2];
    potentialNxtended = new complex<double>[totG4];
    potentialNminus1  = new complex<double>[totG2];

    double x,y,z;
    int idx, idxOnG2, i,j,k,n;

    for( idx = 0; idx < totG2; idx++ ){
        potentialN[idx]     = 0.0;
        potentialNminus1[idx] = 0.0;
    }

    complex<double> imaj(0.,1.);
    double G2shift = 0.0;
    
    for( i = 1; i < xResG2; i++ ){
        for( j = 1; j < yResG2; j++ ){
            for( k = 1; k < zResG2; k++ ){
                
                idxOnG2 = IDX(i, j, k, xResG2, yResG2, zResG2);
                
                x = (i + G2shift)*dx + domainShiftX;
                y = (j + G2shift)*dy + domainShiftY;
                z = (k + G2shift)*dz + domainShiftZ;
                
                complex<double> pot = loader->getLaserBeamEnvelopeProfile(x,y,z,x);// x-ct, t=0
                if( abs(pot) > 1e-2 ){
                    gridMgr->setVectorVariableForNodeG2(idxOnG2, MPOTENTIAL, 0,  real(pot));
                    gridMgr->setVectorVariableForNodeG2(idxOnG2, MPOTENTIAL, 1,  real(pot));
                    gridMgr->setVectorVariableForNodeG2(idxOnG2, MPOTENTIAL, 3,  imag(pot));                
                    gridMgr->setVectorVariableForNodeG2(idxOnG2, MPOTENTIAL, 4,  imag(pot));
                }
                complex<double> pot_shifted = loader->getLaserBeamEnvelopeProfile(x,y,z,x+loader->getTimeStep());// x-c(t-dt), t=0
                if( abs(pot_shifted) > 1e-2 ){
                    gridMgr->setVectorVariableForNodeG2(idxOnG2, MPOTENTIAL, 2,  real(pot_shifted));
                    gridMgr->setVectorVariableForNodeG2(idxOnG2, MPOTENTIAL, 5,  imag(pot_shifted));
                }
            }
        }
    }

    gridMgr->sendBoundary2Neighbor(MPOTENTIAL);
    gridMgr->applyBC(MPOTENTIAL);

    VectorVar** mpotential = gridMgr->getVectorVariableOnG2(MPOTENTIAL);
   
    for( idx = 0; idx < totG2; idx++ ){
        potentialN[idx]       = mpotential[idx]->getValue()[0]+imaj*mpotential[idx]->getValue()[3];
        potentialNminus1[idx] = mpotential[idx]->getValue()[2]+imaj*mpotential[idx]->getValue()[5];
    }


    double* potG2 = new double[2*totG2];
    double* potG4 = new double[2*totG4];

    for( i = 0; i < xRes+2; i++ ){
        for( j = 0; j < yRes+2; j++ ) {
            for( k = 0; k < zRes+2; k++ ){
                idxOnG2 = IDX(i, j, k, xResG2, yResG2, zResG2);
                potG2[2*idxOnG2]   = mpotential[idxOnG2]->getValue()[0];
                potG2[2*idxOnG2+1] = mpotential[idxOnG2]->getValue()[3];
            }
        }
    }
    gridMgr->fillG4variable(potG2, potG4, 2);

    for( idx = 0; idx < totG4; idx++ ){
        potentialNxtended[idx] = potG4[2*idx]+imaj*potG4[2*idx+1];
    }

    delete[] potG2;
    delete[] potG4;


    logger->writeMsg("[EnvelopeEquationSolver] initialize...OK", DEBUG);
}


int EnvelopeEquationSolver::solve(int phase, int i_time){
    auto start_time = high_resolution_clock::now();
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    int xSize = loader->resolution[0];
    int ySize = loader->resolution[1];
    int zSize = loader->resolution[2];

    int totG2 = (xSize+2)*(ySize+2)*(zSize+2);

    int xResG2 = xSize+2, yResG2 = ySize+2, zResG2 = zSize+2;
    int xResG4 = xSize+4, yResG4 = ySize+4, zResG4 = zSize+4;
    int totG4 = xResG4*yResG4*zResG4;

    double ts = loader->getTimeStep();
    complex<double> tsimaj(ts,0.);
    double ts2 = pow(ts,2); 
    complex<double> imaj(0.,1.);
    
    complex<double> potentialNplus1, laplcianOfPotential, gradOfPotential;
    complex<double> chi (0.0,0);
    complex<double> omega(10.0,0);
    complex<double> k0 = omega;

    double delta = ( 1.-pow(ts/dx,2) ) / 3.;
    complex<double> ratio = 1.0;

    int idx, left, rigt, left2, rigt2, up, down, back, front, idxG2pluse, idxG2minus;
    int idxG2,idxG4,i,j,k,neighbour;

    VectorVar** mpotential = gridMgr->getVectorVariableOnG2(MPOTENTIAL);
    VectorVar** density    = gridMgr->getVectorVariableOnG2(DENSELEC);

    if( i_time > loader->laserPulseDuration_tsnum ){
        for( idx = 0; idx < totG2; idx++ ){
            for( int h = 0; h < 6; h++ ){
                gridMgr->setVectorVariableForNodeG2(idx, MPOTENTIAL, h, 0.0);
            }
        }
        return 0;
    }

    double ne;
    double ncritical = loader->criticalDensity;
    for( idx = 0; idx < totG2; idx++ ){
        potentialN[idx] = mpotential[idx]->getValue()[0]+imaj*mpotential[idx]->getValue()[3];
    }

    for( i = 1; i < xSize+1; i++ ){
        for( j = 1; j < ySize+1; j++ ){
            for( k = 1; k < zSize+1; k++ ){                
                idxG2 = IDX(i,  j,  k,  xSize+2,ySize+2,zSize+2);
                idxG4 = IDX(i+1,j+1,k+1,xResG4,yResG4,zResG4);
                left  = IDX(i  ,j+1,k+1,xResG4,yResG4,zResG4);
                rigt  = IDX(i+2,j+1,k+1,xResG4,yResG4,zResG4);
                left2 = IDX(i-1,j+1,k+1,xResG4,yResG4,zResG4);
                rigt2 = IDX(i+3,j+1,k+1,xResG4,yResG4,zResG4);

                back  = IDX(i+1,j  ,k+1,xResG4,yResG4,zResG4);
                front = IDX(i+1,j+2,k+1,xResG4,yResG4,zResG4);
                up    = IDX(i+1,j+1,k+2,xResG4,yResG4,zResG4);
                down  = IDX(i+1,j+1,k  ,xResG4,yResG4,zResG4);
 
                laplcianOfPotential = 0.0;
                laplcianOfPotential += (1.+delta)*(potentialNxtended[left]
                                              -2.0*potentialNxtended[idxG4]
                                                  +potentialNxtended[rigt])/pow(dx,2);
                laplcianOfPotential -= 0.25*delta*(potentialNxtended[left2]
                                              -2.0*potentialNxtended[idxG4]
                                                  +potentialNxtended[rigt2])/pow(dx,2);
                laplcianOfPotential += (potentialNxtended[back]
                                   -2.0*potentialNxtended[idxG4]
                                       +potentialNxtended[front])/pow(dy,2);
                laplcianOfPotential += (potentialNxtended[down]
                                   -2.0*potentialNxtended[idxG4]
                                       +potentialNxtended[up])/pow(dz,2);

                gradOfPotential  = (1.0+delta)*(potentialNxtended[rigt]-potentialNxtended[left ])/2.0/dx;
                gradOfPotential -=  0.5*delta*(potentialNxtended[rigt2]-potentialNxtended[left2])/2.0/dx;

                ne = density[idxG2]->getValue()[0];
                chi = ne / ncritical;
                ratio = (1.0+imaj*k0*tsimaj)/(1.0+pow(k0,2)*ts2);
                potentialNplus1 = ((-chi*potentialN[idxG2]+laplcianOfPotential+2.0*imaj*k0*gradOfPotential)*ts2
                               +2.0*potentialN[idxG2]-(1.0+(imaj*k0)*tsimaj)*potentialNminus1[idxG2])*ratio;

                pumpElectronPressureViaInverseBremsstrahlung( phase, abs(omega), idxG2 );

                gridMgr->setVectorVariableForNodeG2(idxG2, MPOTENTIAL, 0, real(potentialNplus1));
                gridMgr->setVectorVariableForNodeG2(idxG2, MPOTENTIAL, 3, imag(potentialNplus1));
            }
        }
    }

    gridMgr->sendBoundary2Neighbor(MPOTENTIAL);
    gridMgr->applyBC(MPOTENTIAL);


    switch (phase) {
        case PREDICTOR:
            for( idx = 0; idx < totG2; idx++ ){
                double extrapolationReal = -mpotential[idx]->getValue()[1]+2.0*mpotential[idx]->getValue()[0];
                double extrapolationImaj = -mpotential[idx]->getValue()[4]+2.0*mpotential[idx]->getValue()[3];
                gridMgr->setVectorVariableForNodeG2(idx, MPOTENTIAL, 2, mpotential[idx]->getValue()[0]);
                gridMgr->setVectorVariableForNodeG2(idx, MPOTENTIAL, 0, extrapolationReal);
                gridMgr->setVectorVariableForNodeG2(idx, MPOTENTIAL, 5, mpotential[idx]->getValue()[3]);
                gridMgr->setVectorVariableForNodeG2(idx, MPOTENTIAL, 3, extrapolationImaj);
            }
            break;
        case CORRECTOR:
            for( idx = 0; idx < totG2; idx++ ){
                double correctionReal = 0.5*(mpotential[idx]->getValue()[2]+mpotential[idx]->getValue()[0]);
                double correctionImaj = 0.5*(mpotential[idx]->getValue()[5]+mpotential[idx]->getValue()[3]);
                gridMgr->setVectorVariableForNodeG2(idx, MPOTENTIAL, 0, correctionReal);
                gridMgr->setVectorVariableForNodeG2(idx, MPOTENTIAL, 1, correctionReal);
                gridMgr->setVectorVariableForNodeG2(idx, MPOTENTIAL, 3, correctionImaj);
                gridMgr->setVectorVariableForNodeG2(idx, MPOTENTIAL, 4, correctionImaj);
            }
            break;
        default:
            throw runtime_error("no phase");
    }

    for( idx = 0; idx < totG2; idx++ ){
        potentialNminus1[idx] = potentialN[idx];
        potentialN[idx] = mpotential[idx]->getValue()[0]+imaj*mpotential[idx]->getValue()[3];
        double efield = real( (potentialN[idx]-potentialNminus1[idx])/ts - imaj*omega*potentialN[idx] );
        gridMgr->setVectorVariableForNodeG2(idx, EFIELDLASER, 0, efield);
    }

    double* potG2 = new double[2*totG2];
    double* potG4 = new double[2*totG4];

    for( idx = 0; idx < totG2; idx++ ){
        potG2[2*idx]   = mpotential[idx]->getValue()[0];
        potG2[2*idx+1] = mpotential[idx]->getValue()[3];
    }
    gridMgr->fillG4variable(potG2, potG4, 2);

    for( idx = 0; idx < totG4; idx++ ){
        potentialNxtended[idx] = potG4[2*idx]+imaj*potG4[2*idx+1];
    }

    delete[] potG2;
    delete[] potG4;

    auto end_time = high_resolution_clock::now();
    string msg ="[EnvelopeEquationSolver] solve() duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
    
    return 0;
}

void EnvelopeEquationSolver::pumpElectronPressureViaInverseBremsstrahlung(int phase, double omega0, int idxG2){
    VectorVar** density = gridMgr->getVectorVariableOnG2(DENSELEC);
    VectorVar** mpotential = gridMgr->getVectorVariableOnG2(MPOTENTIAL);
    if( abs(mpotential[idxG2]->getValue()[0]) < loader->potentialThreshold4Interaction ){
        return;
    }

    double ne  = density[idxG2]->getValue()[0];
    double ni = ne;
    VectorVar** epres = gridMgr->getVectorVariableOnG2(PRESSURE_SMO);
    double trPe = epres[idxG2]->getValue()[0] + epres[idxG2]->getValue()[3] + epres[idxG2]->getValue()[5];
    double ne_loc = max(ne, EPS8);
    double Te = trPe / (3.0 * ne_loc);
    double z  = max(density[idxG2]->getValue()[2], MIN_CHARGE);//cell charge
    double mi = density[idxG2]->getValue()[3];//cell ion mass
    double mElectron = 0.1;
    double me = mi*mElectron;//cell electron mass
    double logC = 10.0;
    double denom = pow(mi*Te, 1.5);
    denom = max(denom, EPS8);
    double nuei = sqrt(me*mi) * pow(z,2) * ni * logC / denom;
    double ts = loader->getTimeStep();
    double e2 = 0.5*pow(omega0*mpotential[idxG2]->getValue()[0], 2);
    double denomIB = (omega0*omega0 + nuei*nuei);
    double Q_ib = 1.0/mElectron * ne * nuei / denomIB * e2;
    double dE = Q_ib;
    double dPtrace = (2.0 / 3.0) * dE;
    double pressure2Pump = dPtrace / 3.0;

    closureMng->addAlphas(idxG2, 0, pressure2Pump);
    closureMng->addAlphas(idxG2, 3, pressure2Pump);
    closureMng->addAlphas(idxG2, 5, pressure2Pump);
}

EnvelopeEquationSolver::~EnvelopeEquationSolver(){
    finilize();
}

void EnvelopeEquationSolver::finilize(){
    delete[] potentialN;
    delete[] potentialNminus1;
    delete[] potentialNxtended;
}
