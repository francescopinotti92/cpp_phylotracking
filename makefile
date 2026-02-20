# change compiler according to OS
CXX:=clang++
# this is valid on macOS. Remove -undefined dynamic_lookup on Ubuntu
CXXFLAGS:=-O3 -Wall -shared -std=c++11 -undefined dynamic_lookup -fPIC
# type python3 -m pybind11 --includes in the terminal and paste its output here:
INC:=-I/Users/francesco_pinotti/.pyenv/versions/3.10.17/include/python3.10 -I/Users/francesco_pinotti/.pyenv/versions/cmdstanpy310_env/lib/python3.10/site-packages/pybind11/include
# type python3-config --extension-suffix in the terminal and paste its output here:
EXT:=.cpython-310-darwin.so

DEPS:=src/random.cpp src/simulator.cpp src/pysimBD.cpp src/bind_pysimBD.cpp

lib/pysimBD$(EXT): $(DEPS)
	$(CXX) $(CXXFLAGS) $(INC) -I/src $(DEPS) -o lib/pysimBD$(EXT)

clean:
	@rm lib/*$(EXT)
    
.PHONY=clean

