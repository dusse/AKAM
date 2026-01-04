//
//  EnvelopeEquationSolver.hpp

#ifndef EnvelopeEquationSolver_hpp
#define EnvelopeEquationSolver_hpp

#include <chrono>
#include <stdio.h>
#include <iostream>
#include <string>
#include <memory>
#include <complex>

#include "../../../grid/GridManager.hpp"
#include "../../../input/Loader.hpp"
#include "../../../misc/Misc.hpp"
#include "../../../physics/pusher/Pusher.hpp"
#include "../../../physics/pressure-closure/ClosureManager.hpp"

class EnvelopeEquationSolver
{
    
private:
    
    std::unique_ptr<Logger> logger;
    std::shared_ptr<Loader> loader;
    std::shared_ptr<GridManager> gridMgr;
    std::shared_ptr<Pusher> pusher;
    std::shared_ptr<ClosureManager> closureMng;

    std::complex<double>*  potentialN;
    std::complex<double>*  potentialNxtended;
    std::complex<double>*  potentialNminus1; 
    
    void pumpElectronPressureViaInverseBremsstrahlung(int, double, int);

    
public:
    EnvelopeEquationSolver(std::shared_ptr<Loader>, 
                           std::shared_ptr<GridManager>, 
                           std::shared_ptr<Pusher>,
                           std::shared_ptr<ClosureManager>);
        
    void initialize();

    int solve(int, int);
    void finilize();
    ~EnvelopeEquationSolver();
};
#endif
