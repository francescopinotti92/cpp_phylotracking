//
//  pysimBD.cpp
//  BDmodel
//
//  Created by Francesco PINOTTI on 25/09/2025.
//

#include "pysimBD.hpp"


std::string simulate_BD( int seed, int max_cases, int max_samples, double R0, double dI, double rho ) {
    
    m_mt.seed( seed ) ;
    
    Simulator simulator = Simulator( R0, dI, rho ) ;
    
    simulator.set_max_cases( max_cases ) ;
    simulator.set_max_samples( max_samples ) ;
    
    simulator.initialise_single_infection() ;
    
    bool success = simulator.simulate() ;
    if ( success ) {
        
        // The following lines explain how to extract a phylogenetic tree
        LineageTree<int,int>* tree_mngr = simulator.get_tree() ;
        LineageTreeNode<int,int>* rtree = tree_mngr->subSampleTree()[0] ; // obtained the reduced transmission tree (removes all nodes that are not necessary to construct a phylogenetic tree given sampled nodes.
        PhyloNode<int,int>* atree = getAncestralTree( rtree ) ; // extracts the phylogenetic tree
        std::string nwk = getSimpleNewick( atree ) ; // newick string representation of phylogenetic tree
        return nwk ;
        
    }
    else // a simulation may fail due to early extinction
        return "" ;
    
}
