

HDF5_PATH=./libs/hdf5/
MPI_PATH=./libs/openmpi/

PYTHON_INC=./libs/python3.11ssl/include/python3.11/
PYTHON_LIB=./libs/python3.11ssl/lib/


LIBS=-lpython3.11 -lhdf5
INCLUDES=-I$(PYTHON_INC) -I$(HDF5_PATH)/include/ -I$(MPI_PATH)/include/

DSRC = ./src
DEXE = ./

LD_LIBRARY_PATH=$(HDF5_PATH)/lib/:$(PYTHON_LIB)


export LIBRARY_PATH=$LIBRARY_PATH:$(LD_LIBRARY_PATH)

CXX = $(MPI_PATH)/bin/mpicxx

CXXFLAGS  =   -DWRITE_PARTICLES -DGET_ION_PRESSURE -DUSE_EDGE_FACTOR -DLOG -O3 -Wall -c -std=c++11 -Wno-sign-compare -Wno-unused-variable

_SRCS =  $(DSRC)/core/SimulationManager.cpp \
               $(DSRC)/grid/GridManager.cpp \
               $(DSRC)/grid/boundary/BoundaryManager.cpp \
               $(DSRC)/input/Loader.cpp \
               $(DSRC)/particles/Particle.cpp \
               $(DSRC)/misc/Logger.cpp \
               $(DSRC)/misc/Misc.cpp \
               $(DSRC)/output/Writer.cpp \
               $(DSRC)/physics/pusher/Pusher.cpp \
               $(DSRC)/physics/hydro/HydroManager.cpp \
               $(DSRC)/physics/electro-magnetic/EleMagManager.cpp \
               $(DSRC)/physics/pressure-closure/ClosureManager.cpp \
               $(DSRC)/physics/laser/LaserMockManager.cpp \
               $(DSRC)/physics/laser/beam/EnvelopeEquationSolver.cpp \
               $(DSRC)/physics/collisions/IonIonCollisionManager.cpp \
               $(DSRC)/physics/collisions/ElectronIonCollisionManager.cpp \
               $(DSRC)/physics/ionization/IonizationManager.cpp \
               $(DSRC)/common/variables/VectorVar.cpp \
               $(DSRC)/solvers/Solver.cpp \
               $(DSRC)/solvers/ModelInitializer.cpp \
               $(DSRC)/AKA.cpp \

_OBJS            = $(_SRCS:.cpp=.o)

_EXEN            = $(DEXE)/akam.exe

all : $(_EXEN)

$(_EXEN) : $(_OBJS)
	@echo 'Building target: $@'
	$(CXX) -o $@ $^  $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

%.o : %.cpp
	$(CXX) $(INCLUDES) -o $@ $< $(CXXFLAGS) $(FLAGS)

clean :
	rm -f $(_OBJS)
