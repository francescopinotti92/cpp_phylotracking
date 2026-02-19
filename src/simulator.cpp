//
//  simulator.cpp
//  BDmodel
//
//  Created by Francesco PINOTTI on 23/09/2025.
//

#include "simulator.hpp"

void rmv_element( std::vector<int>& v, int ix ) {
    
    std::swap( v[ix], v.back() ) ;
    v.pop_back() ;
    
}

Simulator::Simulator(  double R0_, double dI_, double rho_ ): R0( R0_ ), dI( dI_ ), rho( rho_ ){
    
    I = 0 ;
    t = 0. ;
    
    mu = 1 / dI ;
    beta    = R0 * mu ;
    
    next_lng = 1 ;
    I_lngs = {} ;
    I_lngs.reserve( 10000 ) ;
    n_sampled = 0 ;
    max_cases = 100000000 ;
    max_samples = 10 ;
    
    tree_mngr = new LineageTree<int,int>() ;

}

void Simulator::initialise_single_infection() {
    
    // must call addExtantLineageExternal whenever an introduction event occurs. 'next_lng' is the infected lineage and 't' is infection time. The third entry is just optional metadata: I simply set it to 0 because I am not interested in metadata.
    tree_mngr->addExtantLineageExternal( t, next_lng, 0 ) ;
    I_lngs.push_back( next_lng ) ;
    next_lng++ ;
    I++ ;
    
}

void Simulator::set_max_cases( int max_cases_ ) {
    
    max_cases = max_cases_ ;
}

void Simulator::set_max_samples( int max_samples_ ) {
    
    max_samples = max_samples_ ;
}


bool Simulator::simulate() {
        
    while ( true ) {
        
        double tot_rate = ( beta + mu ) * I ;
        
        if ( tot_rate == 0 ) {
            return false ;
        }
        
        double dt = getExpo( tot_rate ) ;
        
        t += dt ;
        
        double u = getUni() * tot_rate ;
        
        if ( u <= beta * I ) {
            // transmission event
            apply_infection() ;
            
        }
        else {
            // removal event
            apply_removal( rho ) ;
        }
        
        // stopping conditions
        if ( ( next_lng > max_cases ) ) {
            return false ;
        }
        
        if ( n_sampled >= max_samples ) {
            
            return true ;
            
        }
        
    }
        
    return true ;
    
}

void Simulator::apply_infection() {
    
    // update tree by selecting infector from I_lngs
    int ix_infector = getUniInt( I - 1 ) ;
    int lng_infector = I_lngs[ix_infector] ;
    tree_mngr->addExtantLineage( t, next_lng, 0, lng_infector ) ; // must call addExtantLineage whenever a transmission event occurs. 'next_lng' is the name of the lineage created during the transmission event, 'lng_infector' is the parent lineage, 't' is the time of infection. The third entry is just optional metadata: I simply set it to 0 because I am not interested in metadata.
        
    I_lngs.push_back( next_lng ) ;
    next_lng++ ;
    I++ ;
    
}

void Simulator::apply_removal( double prob_sampling ) {
    
    int ix = getUniInt( I - 1 ) ;
    int lng = I_lngs[ix] ;
    
    if ( getBool( prob_sampling ) ) {

        tree_mngr->sampleExtantLineage( lng, t ) ; // marks lineage 'lng' as sampled at time 't'
        n_sampled++ ;

    }
    
    tree_mngr->removeExtantLineage( lng ) ; // must call removeExtantLineage whenever a lineage (lng) is removed from the simulation
    rmv_element( I_lngs, ix ) ;
    I-- ;
    
}
