#include "IonizationManager.hpp"

using namespace std;
using namespace chrono;


IonizationManager::IonizationManager(shared_ptr<Loader> ldr,
                           shared_ptr<GridManager> gridMnr,
                           shared_ptr<Pusher> pshr):loader(move(ldr)),
                           gridMgr(move(gridMnr)), pusher(move(pshr)){
                                
    logger.reset(new Logger());
    initialize();
    logger->writeMsg("[IonizationManager] create...OK", DEBUG);
}

void IonizationManager::initialize(){
    calculateIonizationLevel(PREDICTOR);
    calculateIonizationLevel(CORRECTOR);
}

void IonizationManager::calculateIonizationLevel(int phase) {
    
    auto start_time = high_resolution_clock::now();    
    string msg0 ="[IonizationManager] start to calculate ionization level ";
    logger->writeMsg(msg0.c_str(), DEBUG);

    Particle** particles = pusher->getParticles();
    int totalPrtclNumber = pusher->getTotalParticleNumber();
    
    int numOfSpecies = loader->getNumberOfSpecies();
    
    int posShift = 0;
    if( phase == CORRECTOR ){
        posShift = 3;
    }

    VectorVar** epres = gridMgr->getVectorVariableOnG2(PRESSURE);
    VectorVar** edens = gridMgr->getVectorVariableOnG2(DENSELEC);
    
    double dx = loader->spatialSteps[0];
    double dy = loader->spatialSteps[1];
    double dz = loader->spatialSteps[2];
    
    double domainShiftX = loader->boxCoordinates[0][0];
    double domainShiftY = loader->boxCoordinates[1][0];
    double domainShiftZ = loader->boxCoordinates[2][0];
    
    int xSize = loader->resolution[0];
    int ySize = loader->resolution[1];
    int zSize = loader->resolution[2];
    
    int xSizeG2 = xSize + 2;
    int ySizeG2 = ySize + 2;
    int zSizeG2 = zSize + 2;
    
    double ts = loader->getTimeStep();

    int i, j, k, idxG2;
    double x, y, z;
    double G2shift = 0.5;
    int type;
    double maxCharge4type = 0.0;
    for( int idx = 0; idx < totalPrtclNumber; idx++ ) {
        
        Particle* p = particles[idx];
        type = p->getType();
        maxCharge4type = pusher->getParticleCharge4Type(type);
        if( pusher->getIfParticleTypeIsFrozen(type) == 1 ) continue;
        double* pos = p->getPosition();
        x = pos[0+posShift];
        y = pos[1+posShift];
        z = pos[2+posShift];    
        i = int((x - domainShiftX)/dx + G2shift);
        j = int((y - domainShiftY)/dy + G2shift);
        k = int((z - domainShiftZ)/dz + G2shift);
        
        if( i<0 || j<0 || k<0 || i>=xSizeG2 || j>=ySizeG2 || k>=zSizeG2 ) continue;
        
        idxG2 = IDX(i, j, k, xSizeG2, ySizeG2, zSizeG2);
        
        double ne = edens[idxG2]->getValue()[0];
        if( ne < EPS8 ) continue;
        VectorVar* pres = epres[idxG2];
        double trPe = pres->getValue()[0] + pres->getValue()[3] + pres->getValue()[5];
        double pe = trPe / 3.0;
        double Te = ne < EPS8 ? 0.0 : pe*edgeProfile(ne)/ne;
        
        double ionRate    = computeIonizationRate(pe);
        double recombRate = computeRecombinationRate(p, Te, ne);
        
        double Pion = 1.0 - exp(-ionRate * ts);
        double Prec = 1.0 - exp(-recombRate * ts);
        
        int deltaZ_step = 0;
        
        double r = RNM;
        
        if( p->getCharge() < maxCharge4type && r < Pion ){
            deltaZ_step = +1;
        } else if( p->getCharge() > MIN_CHARGE && r < Prec ){
            deltaZ_step = -1;
        }
        
        double newZ = p->getCharge() + double(deltaZ_step);
        newZ = round(newZ);
        if( newZ < MIN_CHARGE ) {
            newZ = MIN_CHARGE;
        }
        if( newZ > maxCharge4type ){
            newZ = maxCharge4type;
        }
        if( pe > EPS8 ){//do not ionize vacuum
            p->setCharge(newZ);
        }
    }

    auto end_time = high_resolution_clock::now();
    string msg ="[IonizationManager] calculateIonizationLevel()... duration = "
                +to_string(duration_cast<milliseconds>(end_time - start_time).count())+" ms";
    logger->writeMsg(msg.c_str(), DEBUG);
}


double IonizationManager::computeIonizationRate(double pe) {
    const double nu0 = 1e3;
    double rate = pe*nu0;
    return rate;
}


double IonizationManager::computeRecombinationRate(Particle* p, double Te, double ne) {
    double Z = p->getCharge();
    if( Z <= MIN_CHARGE ) return 0.0;
    double Trecomb =  0.1 * (Z + 1.0);
    double rate = ne * Z * exp(-Te / Trecomb);
    const double nu0 = 0.05;
    rate *= nu0;
    return rate;
}