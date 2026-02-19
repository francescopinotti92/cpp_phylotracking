//
//  simulator.hpp
//  BDmodel
//
//  Created by Francesco PINOTTI on 23/09/2025.
//

#ifndef simulator_hpp
#define simulator_hpp

#include "tree.hpp"
#include "random.hpp"
#include <stdio.h>
#include <vector>

void rmv_element( std::vector<int>& v, int ix ) ;

class Simulator {
public:
    Simulator( double R0, double dI, double rho ) ;
    void initialise_single_infection() ;
    void set_max_cases( int max_cases ) ;
    void set_max_samples( int max_samples ) ;
    
    bool simulate() ;
    void apply_infection() ;
    void apply_removal( double prob_sampling ) ;
    
   
    LineageTree<int,int>* get_tree() { return tree_mngr ; }

private:
   
    int I ;
    double t ;
    double R0 ;
    double dI ;
    double rho ;
    double mu ;
    double beta ;
    
    int next_lng ;
    std::vector<int> I_lngs ;
    int n_sampled ;
    int max_cases ;
    int max_samples ;

    /*
     LineageTree<T,U> manages the transmission tree
     T is the type associated with lineage identifiers (int here)
     U is the type associated with lineage metadata (int here, but not used)
     
     If T is a simple type like int, then no need to do anything else.
     If T is a more complicated type, e.g. a custom class, you must write some more code (see Example below)
     */
    LineageTree<int,int>* tree_mngr ;


} ;

/*

 EXAMPLE: CUSTOM STRUCTURE FOR LINEAGE IDENTITY
 
 struct LineageInfo uses 2 integers to identify a lineage (here it was strain ID and chicken ID).
 Note that I overloaded the == operator so that c++ knows how to tell two lineages apart.
 
 hash<LineageInfo> creates a hash function so that c++ knows how to put LineageInfo objects in maps.
 
 Must overload the operator << as well so that c++ knows how to generate a string from LineageInfo.
 
 
 
struct LineageInfo { // goes into tree tracking

    LineageInfo( const int& host_id_ = -1, const int& strain_id_ = -1 ) : host_id( host_id_ ), strain_id( strain_id_ ) {} ;
    int host_id ;
    int strain_id ;
    
    bool operator==(const LineageInfo& other) const { return host_id == other.host_id && strain_id == other.strain_id ; }
    
} ;


// Define hash specialization inside std namespace
namespace std {
    template<>
    struct hash<LineageInfo> {
        std::size_t operator()(const LineageInfo& li) const noexcept {
            std::size_t h1 = std::hash<int>()(li.host_id);
            std::size_t h2 = std::hash<int>()(li.strain_id);
            return h1 ^ (h2 << 1);
        }
    };
}

 std::ostream& operator<<( std::ostream& os, LineageInfo& li ) {
    os << li.host_id << "-" << li.strain_id ;
    return os ;
 }
 
 
 */
#endif /* simulator_hpp */
