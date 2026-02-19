//
//  pysimBD.hpp
//  BDmodel
//
//  Created by Francesco PINOTTI on 25/09/2025.
//

#ifndef pysimBD_hpp
#define pysimBD_hpp

#include "simulator.hpp"
#include <stdio.h>
#include <string>

// simulated a birth-death model with basic reproduction number R0, duration of infection dI and sampling probability rho
// max_samples is the desired number of lineages
// max_cases sets a further stopping condition depending on the total number of cases: just set it to a very large number
std::string simulate_BD( int seed, int max_cases, int max_samples, double R0, double dI, double rho ) ;


#endif /* pysimBD_hpp */
