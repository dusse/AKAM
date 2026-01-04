#ifndef LaserMockManager_hpp
#define LaserMockManager_hpp

#include <stdio.h>

#include <chrono>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <string>
#include <memory>
#include<algorithm>

#include "../../grid/GridManager.hpp"
#include "../../particles/Particle.hpp"
#include "../pusher/Pusher.hpp"
#include "../../input/Loader.hpp"
#include "../../misc/Misc.hpp"
#include "../../common/variables/VectorVar.hpp"
#include "../pressure-closure/ClosureManager.hpp"
#include "beam/EnvelopeEquationSolver.hpp"

class LaserMockManager{
    
private:
    std::unique_ptr<Logger> logger;
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMgr;
    std::shared_ptr<Pusher> pusher;
    std::shared_ptr<EnvelopeEquationSolver> envelopeEqSolver;
    
    int PARTICLE_TYPE2LOAD = 1;
    
    double* electronPressureProfile;
    double* targetIonDensityProfile;
    double* ionThermalVelocityProfile;
    double* ionFluidVelocityProfile;
    
    
    void initialize();

    
public:
    LaserMockManager(std::shared_ptr<Loader>,
                     std::shared_ptr<GridManager>,
                     std::shared_ptr<Pusher>,
                     std::shared_ptr<ClosureManager>);

    ~LaserMockManager();
    
    void addIons(int);
    void accelerate(int);
    void moveBeamEnvelope(int, int);
    

};


#endif /* LaserMockManager_hpp */
