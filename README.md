
________________________________________________
[0] [1] ****                                                   ****.
[0] [1] ****       O        0  A       A         M          M  ****.
[0] [1] ****      O O       O  A      A A        M M      M M  ****.
[0] [1] ****     O   O      O A      A   A       M M      M M  ****.
[0] [1] ****    OOOOOOO     0A      AAAAAAA      M   M   M  M  ****.
[0] [1] ****   O       0    O A    A       A     M     M    M  ****.
[0] [1] ****  O         O   O  A  A         A    M          M  ****.
[0] [1] ****   Arbitrary  Kinetic  Algorithm           M       ****.
[0] [1] ***********************************************************.
________________________________________________

# AKAM
AKAM is a three-dimensional hybrid PIC–fluid code for multiphysics plasma modeling. Ions are treated kinetically, while electrons are described by a ten-moment fluid model, enabling efficient simulation of ion-scale dynamics, pressure anisotropy, and non-Maxwellian effects. Laser–plasma interaction is modeled using a laser envelope formulation that captures energy deposition, electrons ponderomotive heating, and ions collisional absorption. The framework is designed for scalable, end-to-end simulations of laser-driven and high-energy-density plasmas.

 stack: C++11, MPI, HDF5, python2-3

**Paper reference:** [Andrey Sladkov, *Hybrid PIC-fluid model for numerical simulation of laser-plasma interaction*, 2026]


_______________________

#       HOWTO
_______________________
1. before 'make' need to set in makefile  
    HDF5_PATH= path to hdf5 lib (last well used 1.10.5)  
    MPI_PATH= path to mpi lib (last well used openmpi 9.0.0)  
    PYTHON_INC= path to python include  
    PYTHON_LIB= path to python lib  

2. for running default example from src/input/Initializer.py  
    mpirun -n 2 akam.exe

3. normally need to set input file path containing Initializer.py  
    mpirun -n 2 akam.exe PATH/TO/PYTHON/INPUT/FILE

4. before running need to create output folder and set in Initializer.py  

5. for visualization use python notebook in folder NOTEBOOK/  

_______________________

# physics included:
_______________________
* space and time quantities are normalized on:  
  - ion inertia length d0 = c/omega_pi  
  - inverse ion gyrofrequency 1/Omega_ci  

* electro-magnetic fields are calculated on two staggered grid (G1 and G2)  
  using predictor-corrector scheme (see EleMagManager.cpp)  

  E = -[VixB] + [JxB]/n - divPe/n - eta ΔJ   

  B' = -rotE  

  J = rotB  

* ions are described in a kinetic way,   
  dynamics is solved using first order interpolation of em fields  

* electrons are described in a fluid way by:   
  - density (equals to ion density ne=ni=n)   
  - bulk velocity Ve = Vi - J/n  
  - six-component pressure tensor Pij  

* six-component pressure tensor P is integrated in time   
  using subcycling explicit scheme  

  P' = - Ve.∇P - P∇.Ve - P.∇Ve - (P.∇Ve)^T - q/m [ PxB + (PxB)^T ]   

  where Ve - electron flow velocity  
        q - electron charge  
        m - electron mass  
        B - magnetic field  
 

* laser–plasma interaction is done via envelope model with ponderomotive heating and inverse Bremsstrahlung absorption

* electron–ion collisions are included through a Landau–Spitzer relaxation operator, enabling temperature equilibration and pressure feedback between ions and the electron fluid

* ionisation processes are included via rate-based evolution of particle charge states driven by local electron density and temperature.



_______________________
#     FEATURES
_______________________
1. AKAM is 3D, parallel (MPI), multispecies code 

2. BC type for em fields and hydro quantities: 1 - periodic, 0 - ideal (see GridManager.cpp)

3. BC type for particle properties: 2 - reflect 1 - periodic 0 - outflow (see BoundaryManager.cpp)
   *outflow BC - reaching border particle leaves domain forever

4. for small scale dissipation use resistivity (eta) parameter

5. for pressure tensor integration need to set:
   * electron mass
   * relaxation factor for izotropization operator (see ClosureManager.cpp)
   * smooth stride for pressure tensor smoothing

6. use 'make FLAGS=-DLOG' to set debug log level in Logger.hpp
