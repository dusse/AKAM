#ifndef ElectronIonCollisionManager_hpp
#define ElectronIonCollisionManager_hpp

#include <stdio.h>

#include <chrono>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <string>
#include <memory>
#include <random>

#include "../../grid/GridManager.hpp"
#include "../pusher/Pusher.hpp"
#include "../../input/Loader.hpp"
#include "../../misc/Misc.hpp"
#include "../pressure-closure/ClosureManager.hpp"


class ElectronIonCollisionManager{
    
private:
    std::unique_ptr<Logger> logger;
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMgr;
    std::shared_ptr<Pusher> pusher;
    std::shared_ptr<ClosureManager> closureMng;
    
    double* ionTrPbeforeCollision;
    void initialize();

    void scatterVelocities(int, int, int, int, int,
                           double*, double* , double);
    
public:
    ElectronIonCollisionManager(std::shared_ptr<Loader>,
                           std::shared_ptr<GridManager>,
                           std::shared_ptr<Pusher>,
                           std::shared_ptr<ClosureManager>);


    ~ElectronIonCollisionManager();
    void collideElectronWithIons(int);
    void applyElectronPressure(int);
};


#endif /* ElectronIonCollisionManager_hpp */
