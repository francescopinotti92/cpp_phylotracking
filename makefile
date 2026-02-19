#==Local

CXX:=clang++
CXXFLAGS:=-O3 -Wall -shared -std=c++11 -undefined dynamic_lookup -fPIC
# use output of python3 -m pybind11 --includes
INC:=-I/Users/francesco_pinotti/.pyenv/versions/3.10.17/include/python3.10 -I/Users/francesco_pinotti/.pyenv/versions/cmdstanpy310_env/lib/python3.10/site-packages/pybind11/include -I/src
#INC:=-I/Users/user/miniconda3/envs/elfi_env/include/python3.12 -I/Users/user/miniconda3/envs/elfi_env/lib/python3.12/site-packages/pybind11/include -I/src
# use output of python3-config --extension-suffix
EXT:=.cpython-310-darwin.so

DEPS:=src/random.cpp src/simulator.cpp src/pysimBD.cpp src/bind_pysimBD.cpp

lib/pysimBD$(EXT): $(DEPS)
	$(CXX) $(CXXFLAGS) $(INC) $(DEPS) -o lib/pysimBD$(EXT)

clean:
	@rm lib/*$(EXT)
    
.PHONY=clean

