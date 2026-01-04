#ifndef IonizationManager_hpp
#define IonizationManager_hpp

#include <stdio.h>

#include <chrono>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <string>
#include <memory>
#include <algorithm>

#include "../../grid/GridManager.hpp"
#include "../../particles/Particle.hpp"
#include "../pusher/Pusher.hpp"
#include "../../input/Loader.hpp"
#include "../../misc/Misc.hpp"
#include "../../common/variables/VectorVar.hpp"


class IonizationManager{
    
private:
    std::unique_ptr<Logger> logger;
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMgr;
    std::shared_ptr<Pusher> pusher;
    
    void initialize();
    double computeIonizationRate(double);
    double computeRecombinationRate(Particle*, double, double);
    
public:
    IonizationManager(std::shared_ptr<Loader>, std::shared_ptr<GridManager>, std::shared_ptr<Pusher>);
    void calculateIonizationLevel(int);
    

};
#endif /* IonizationManager_hpp */
