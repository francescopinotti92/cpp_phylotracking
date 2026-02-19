//
//  tree.hpp
//  BDmodel
//
//  Created by Francesco PINOTTI on 23/09/2025.
//

#ifndef tree_h
#define tree_h

#include <cassert>
#include <sstream>
#include <ostream>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>


//====== LineageTreeNode ======//

template <typename T, typename U, class Hash = std::hash<T>>
struct LineageTreeNode {
    
    LineageTreeNode( const T& lng, const U& data, const double& t, const bool& extant, LineageTreeNode<T,U>* parent = nullptr  ): t( t ), tSample( 0. ), tBranchParent( t ), locSample("NA"), lng( lng ), data( data ), extant( extant ), parent( parent ), needed( false ), sampled( false ) {
        
        children = {} ;
        //children_branching_times = {} ;
        
    }
    
    void eraseChild( LineageTreeNode<T,U>* child ) { // removes child from children (does not perform further updates though)
        
        auto it = children.begin() ;
        while( it != children.end() ) {
            
            if ( (*it)->lng == child->lng ) {
                std::swap( *it, children.back() ) ;
                children.pop_back() ;
                break ;
            }
            ++it ;
            
        }
        
    } ;
    
    uint getSizeChildren() {
        
        return static_cast<uint>( children.size() ) ;
        
    }
    double t ; // birth time
    double tSample ; // sampling time
    double tBranchParent ; // time at which lineage branched from parent node (not necessarily the true parent)
    std::string locSample ; // sampling location
    T lng ;
    U data ;
    LineageTreeNode<T,U>* parent ;
    std::vector<LineageTreeNode<T,U>*> children ;
    bool extant ; // true if still around in simulation
    bool needed ; // true if required in reduced transmission tree
    bool sampled ; // true if sampled
    //std::unordered_map<T, double, Hash > children_branching_times ;
    
    
} ;

/*
 Frees memory allocated to a LineageTreeNode<T,U> tree.
 */

template<typename T, typename U>
void deleteLineageTreeNodeTree( LineageTreeNode<T,U>* root ) {
    
    while( not root->children.empty()  ) {
        
        deleteLineageTreeNodeTree( root->children.back() ) ;
        root->children.pop_back() ;
        
    }
    
    delete root ;
    
}


//====== LineageTree ======//

template <typename T, typename U, class Hash = std::hash<T>>
class LineageTree {
public:
    /*
     
     Constructor (creates an empty tree).
     
     */
    LineageTree(): nnodes( 0 ) {
        
        extantLngs = {} ;
        roots.clear() ;
        sampled_lineages.clear() ;
        //parent_info = {} ;
        
    } ;
    
    
    /*
     Reset function.
     
     SHOULD BE USED WHEN THE SAME LINEAGE TREE INSTANCE
     IS RE-USED ACROSS INDEPENDENT SIMULATIONS LAUNCHED
     SEQUENTIALLY.
     
     */
    void reset() {
               
        /*

        auto extantLngsCopy = extantLngs ; // Use copy to avoid errors during lineage removal
        
        for ( auto leaves : extantLngsCopy ) {
            
            T lng = leaves.first ;
            removeExtantLineage( lng, true ) ;
            
        }

        sampled_lineages.clear() ; */
        
        for ( auto root : this->roots )
            deleteLineageTreeNodeTree( root ) ;
        
        roots.clear() ;
        extantLngs.clear() ;
        sampled_lineages.clear() ;
        //parent_info.clear() ;
        
        nnodes = 0 ;
      
        // ??? persistent reminder to avoid memory leaks
    } ;
    
    /*
     
     Adds a lineage 'lng' born at time 't' with parent 'lngParent' with
     metadata 'data'.
     
     THIS FUNCTION SHOULD BE CALLED AFTER A TRANSMISSION EVENT.
     
     */
    void addExtantLineage( const double& t, const T& lng, const U& data, const T& lngParent )  {
        
        LineageTreeNode<T,U>* lngParentNode = extantLngs[lngParent] ;
        LineageTreeNode<T,U>* lngNode = new LineageTreeNode<T,U>( lng, data, t, true, lngParentNode ) ;
        lngParentNode->children.push_back( lngNode ) ;
        //lngParentNode->children_branching_times[lng] = t ;
        extantLngs[lng] = lngNode ;
        ++nnodes ;
        
    }
    
    /*
     
     Adds a lineage 'lng' born at time 't' with metadata 'data' and no parent.
     
     THIS FUNCTION SHOULD BE USED WHEN AN EXTERNAL INTRODUCTION OCCURS.
     
     */
    void addExtantLineageExternal( const double& t, const T& lng, const U& data )  {
                
        LineageTreeNode<T,U>* lngNode = new LineageTreeNode<T,U>( lng, data, t, true, nullptr ) ;
        extantLngs[lng] = lngNode ;
        roots.insert( lngNode ) ;
        ++nnodes ;
        
    } ;
    
    /*
     
     Removes lineage 'lng' that became extinct.
     
     'lng' could be EXTREMAL (i.e. with no children) or INTERNAL within the tree.
     
     'lng' could be SAMPLED or UNSAMPLED.
     
     If 'lng' is NOT SAMPLED, do the following:
     
        Check if 'lng' is EXTREMAL or INTERNAL:
     
            a. if EXTREMAL, check if lng is ROOT
                i.. if YES, remove 'lng' from the root list
                ii. if NOT, notify its parent about 'lng' removal (call 'notifyParent').
            b. if INTERNAL with exactly one child, prune 'lng'
               and attach its only child to its parent (call 'mergeParent').
     
    */
    void removeExtantLineage( const T& lng, bool ignore_sampled = false )  {
        
        LineageTreeNode<T,U>* lngNode = extantLngs[lng] ;
        lngNode->extant = false ;
        
        bool proceed = true ;
        if ( lngNode->sampled )
            proceed = false ;
        
        if ( proceed ) { // remove node only if not sampled
            
            if ( ( lngNode->children ).size() == 0 ) { // has no extant children
                
                if ( lngNode->parent != nullptr ) // if parent is not ROOT, broadcast removal upstream
                    notifyParent( lngNode->parent, lngNode, ignore_sampled ) ;
                else
                    roots.erase( lngNode ) ; // remove from root
                
                delete lngNode ;
                --nnodes;
                
            }
            
            // check if merge is possible
            else if ( ( lngNode->children ).size() == 1 )
                mergeParentChild( lngNode ) ;
            
            // else, must keep node
            
        }
        
        extantLngs.erase( lng ); // lng is not extant anymore, hence update extant list
        
    } ;
    
    
    /*
     Marks lineage 'lng' as SAMPLED.
     
     Also adds additional info about the time of sampling 't' and the location
     of sampling 'locSample' (optional).
     
     Returns 'true' if the lineage is sampled successfully. Returns 'false' if not,
     i.e. if it had been sampled already. This prevents a lineage from being sampled twice
     and is relevant in models where sampled lineages are not removed right after sampling.
     
     */
    
    bool sampleExtantLineage( const T& lng, const double& t, const std::string& locSample = "@" ) {
        
        if ( extantLngs[lng]->sampled ) // lng has already been sampled
            return false ;
        else { // lng has not been sampled already
            
            extantLngs[lng]->sampled = true ;
            extantLngs[lng]->tSample = t ;
            extantLngs[lng]->locSample = locSample ;
            
            sampled_lineages.insert( lng ) ;
            return true ;
            
        }

    } ;
    
    
    /*
     
     Retrieves all SAMPLED lineages descending from 'rootNode', i.e. all lineages
     that were marked as SAMPLED via 'sampleExtantLineage'.
     
     N.B.The result may include extinct sampled lineages too.
     
     */

    std::vector<T> getSampledLineages( LineageTreeNode<T,U>* rootNode )  {
        
        // return empty vector if not root
        if ( rootNode->parent != nullptr )
            return {} ;
        
        // find all sampled lineages recursively
        std::vector<T> sampledLngs = {} ;
        getSampledLineagesRecursive( rootNode, sampledLngs ) ;
        
        return sampledLngs ;
        
    } ;
    
    /*
     
     Returns 'true' if 'lng' has already been sampled.
     It is useful for external code in order to avoid sampling the same lineage twice
     
     */
    bool is_lineage_sampled( const T& lng ) {
        
        if ( sampled_lineages.find( lng ) != sampled_lineages.end() )
            return true ;
        else
            return false ;
        
    }


    //std::vector<LineageTreeNode<T,U>*> subSampleTree( const std::unordered_map<T,DataLineageSampling,Hash>& sampledLngsInfo ) ;
    
    
    /*
     
     Yields the reduced transmission tree from the full transmission tree.
     
     Keeps only sampled EXTREMAL (and INTERNAL) nodes that were marked as
     SAMPLED or that are instrumental to draw the ancestry of SAMPLED
     lineages.
     
     Returns a vector with multiple trees corresponding to disjoint trees
     (happens when sampled lineages descend from distinct introductions)
     
     */
    
    std::vector<LineageTreeNode<T,U>*> subSampleTree()  {
        
        // loop over root nodes
        std::vector<LineageTreeNode<T,U>*> res = {} ;
        for ( LineageTreeNode<T,U>* rootNode : roots ) {
            
            // find sampled lngs descending from root
            std::vector<T> selectedLngs = getSampledLineages( rootNode );
            
            uint nSampledLngs = static_cast<uint>( selectedLngs.size() ) ;
            
            if ( nSampledLngs > 0 ) {
                
                // extract subtree..
                markNodeNeeded( rootNode, selectedLngs ) ;
                LineageTreeNode<T, U>* subTreeRoot = extractSubTree( rootNode, nullptr ) ;
                subTreeRoot = eliminateRedundantNodes( subTreeRoot, selectedLngs ) ;
                //printEdgesFromNode( subTreeRoot ) ;

                res.push_back( subTreeRoot ) ;
                
            }

        }
        
        return res ;
        
    } ;

    /*
     
     Returns the ROOT node of the tree containing 'lngNode'.
     
     N.B. if 'lngNode' is ROOT, it returns itself.
     
     */
    
    LineageTreeNode<T,U>* getRootNode( LineageTreeNode<T,U>* lngNode )  {
        
        if ( lngNode->parent == nullptr )
            return lngNode ;
        else
            return getRootNode( lngNode->parent ) ;
        
    } ;

    /*
     
     Returns the number of extant lineages.
     
     */
    
    uint getSizeExtantLineages() {
        
        return static_cast<uint>( extantLngs.size() ) ;
        
    }
    
    /*
     
     Returns the number of nodes.
     
     */
    
    uint getSizeNodes() { return nnodes ; }
    
private:
    uint nnodes ;
    std::unordered_map<T, LineageTreeNode<T,U>*, Hash > extantLngs ; // list of extant lineages
    std::unordered_set<LineageTreeNode<T,U>*> roots ; // list of roots, i.e. trees
    std::unordered_set<T> sampled_lineages ;
    //std::unordered_map<LineageTreeNode<T,U>*, std::pair<LineageTreeNode<T,U>*, double>> parent_info ;
    
    /*
          
     Notifies 'parent' lineage that 'child' lineage went extinct.
     Uses a recursive strategy to broadcast the signal.
     Here 'parent' is the current node.
     
     Decides whether 'parent' is made redundant after removing 'child'.
     The outcome depends on whether 'child' was SAMPLED and on whether
     'parent' is EXTANT and on the number of its children.
     
     - Remove 'child' if it is NOT SAMPLED.
     - If 'parent' is EXTANT, stop (it cannot be removed).
     - If 'parent' is EXTINCT but SAMPLED, stop (it cannot be removed).
     - If 'parent' is EXTINCT but NOT SAMPLED, check how many children M.
        + If M = 0 (no children), check if grandparent node is ROOT:
            > if NOT ROOT, notify grandparent.
            > if ROOT, remove from root list.
            > In either case, free memory & update 'nnodes'.
        + if M = 1 (exactly one child), make merge move ('mergeParentChild').
        + if M > 1 do nothing.
     
        1. If YES, do nothing
        2. If NOT, check number of children:
            a. If childless, check if root:
                i.  If NOT, notify grandparent node
                ii. If YES, remove from root list
                iii. In both cases, free memory and update nnodes
            b. If only one child, do merge move
            c. Else, do nothing
     
     */
    
    void notifyParent( LineageTreeNode<T,U>* parent, LineageTreeNode<T,U>* child, bool ignore_sampled = false ) {
        
        bool parentExtinct = !parent->extant   ;
        bool parentSampled = parent->sampled   ;
        bool childSampled  = child->sampled    ;
        
        if ( !childSampled ) { // erase child only if unsampled
            // this should be a full removal
            //parent->children_branching_times.erase( child->lng ) ;
            parent->eraseChild( child ) ;
        }
        
        if ( parentExtinct ) { // parent is extinct, check if it was also sampled
            
            /*
            bool proceed = true ;
            if ( parentSampled )
                proceed = false ;*/
            
            if ( !parentSampled ) {
                
                // parent was not sampled
                
                uint nChildren = parent->getSizeChildren() ;
                
                if ( nChildren == 0 ) {
                    
                    // after removing child, parent is a redundant extinct leaf: remove it
                    bool parentRoot = parent->parent == nullptr ; //
                    if ( !parentRoot ) {
                        
                        // notify grandparent
                        notifyParent( parent->parent, parent, ignore_sampled ) ;
                        
                    }
                    else {
                        
                        // parent is also root: remove from root list
                        roots.erase( parent ) ;
                        
                    }
                    
                    delete parent ; // free memory
                    --nnodes ;
            
                }
                else if ( nChildren == 1 ) {
                    // parent is a redundant mid node: do merge move
                    mergeParentChild( parent ) ;
                    
                }
                // else do nothing

            }
            
        }
            
    } ;
    
    
    /*
     
     Merges the only child of 'midNode' with its parent, then
     removes 'midNode'.
     
     This function should not be called on nodes that:
     
         1. Have not exactly one child
         2. Are EXTANT.
         3. Have been SAMPLED.
     
     See also 'notifyParent'.
     
     */
    
    void mergeParentChild( LineageTreeNode<T,U>* midNode )  {
        
        assert( midNode->children.size() == 1 ) ;
        assert( !midNode->extant ) ;
        
        if ( midNode->parent != nullptr ) { // midNode is an intermediate node O->X->O
            
            midNode->children[0]->parent = midNode->parent ;
            midNode->parent->eraseChild( midNode ) ;
            midNode->parent->children.push_back( midNode->children[0] ) ;
            //midNode->parent->children_branching_times[ midNode->children[0]->lng ] = midNode->parent->children_branching_times[ midNode->lng ] ;
            //midNode->parent->children_branching_times.erase( midNode->lng ) ;

            midNode->children[0]->tBranchParent = midNode->tBranchParent ;
            
        }
        else { // midNode is a root with a single child @->X->O, hence child becomes root
            
            midNode->children[0]->parent = nullptr ;
            roots.erase( midNode ) ;
            roots.insert( midNode->children[0] ) ;
            //midNode->children_branching_times.erase( midNode->lng ) ;
            midNode->children[0]->tBranchParent = midNode->children[0]->t ; // This should be OK because branching time is irrelevant for roots
        
        }
        
        delete midNode ;
        --nnodes ;
        
    } ;
    
    
    /*
     
     // ??? Revise the description
     
     This is a recursive function that when called on a root node of a transmission tree,
     yields the subtree containing only nodes marked as 'needed' (i.e. nodes that had been
     sampled or that are instrumental to reconstruct ancestral relationships between sampled
     lineage).
     
     The function is called on the 'node' of the original transmission tree, while 'parent' is
     a copy of node->parent in the subtree.
     
     Extracts subtree using the 'needed' attribute.
     
     Starting from a current tree node ('node') and attach node to 'parent', which is a newly created node
     that builds the output, the subsampled tree
     
     */
    
    LineageTreeNode<T,U>* extractSubTree( LineageTreeNode<T,U>* node, LineageTreeNode<T,U>* parent ) {
        
        assert( node != nullptr ) ;
        
        // copy original node into new node (output)
        LineageTreeNode<T,U>* newNode = new LineageTreeNode<T,U>( node->lng, node->data, node->t, node->extant, parent ) ;
        
        // Add rest of sampling information
        newNode->sampled   = node->sampled ;
        newNode->tSample   = node->tSample ;
        newNode->locSample = node->locSample ;
        newNode->tBranchParent = node->tBranchParent ;
        //newNode->children_branching_times = {} ;

        for ( auto& child : node->children ) {
            
            if ( child->needed ) {
                
                LineageTreeNode<T,U>* newChild = extractSubTree( child, newNode ) ;
                ( newNode->children ).push_back( newChild ) ;
                //newNode->children_branching_times[ child->lng ] = node->children_branching_times[ child->lng ] ;
                
            }
            
        }
        
        return newNode ;
        
    } ;
    
    /*
     
     Recursively fills 'lngs' vector with SAMPLED lineages.
     
     If 'node' is SAMPLED, store it back. Then move to its children.
     
     */
    void getSampledLineagesRecursive( LineageTreeNode<T,U>* node, std::vector<T>& lngs ) {
        
        if ( node->sampled ) // if sampled, add node to vector
            lngs.push_back( node->lng ) ;
        
        for ( auto child : node->children )
            getSampledLineagesRecursive( child, lngs ) ;
        
    } ;

    
    /*
     
     Set 'node' to 'needed' if it should be included in the reduced transmission tree.
     
     This function is used recursively: returns whether 'node' is 'needed' or not.
     At the same time it fills the vector 'neededLngs' with SAMPLED nodes and some of their
     ancestors.
     
     */
    
    bool markNodeNeeded( LineageTreeNode<T,U>* node, const std::vector<T>& neededLngs ) {
        
        if ( node->extant ) { // is extant
            
            bool needed = false ;
            for ( const auto& ll : neededLngs ) {
                
                if ( ll == node->lng ) {
                    
                    needed = true ;
                    break ;
                    
                }
            
            }
            
            if ( node->sampled ) // if sampled previously
                needed = true ;
            
            if ( needed ) {
                node->needed = true ;
                for ( auto& child : node->children )
                    markNodeNeeded( child, neededLngs ) ;
            }
            else { // if not sampled directly
                
                for ( auto& child : node->children ) {
                    
                    bool isChildNeeded = markNodeNeeded( child, neededLngs ) ;
                    node->needed = ( node->needed or isChildNeeded ) ;
                    
                }
                
            }
            
            return node->needed ;
        
        }
        else { // is not extant, check if needed due to children or if sampled
            
            node->needed = false ;
            
            if ( node->sampled )
                node->needed = true ;
            
            for ( auto& child : node->children ) {
                
                bool isChildNeeded = markNodeNeeded( child, neededLngs ) ;
                node->needed = ( node->needed or isChildNeeded ) ;
                
            }
            
        }
        
        return node->needed ;
        
    } ;

    /*
     
     Helper function to print a transmission tree, edge by edge.
     
     */
    /*
    void printParentChildEdge( LineageTreeNode<T,U>* child, bool recursive = true )  {
        
        std::string stringParent = "@";
        if ( child->parent != nullptr )
            stringParent = utils::to_string( child->parent->lng ) ;
            
        std::cout << stringParent << "->" << utils::to_string( child->lng )<< "\n";
        
        if ( recursive and child->parent != nullptr )
            printParentChildEdge( child->parent, true ) ;
        
    } ;*/

    
    
} ;

//====== LineageTreeNode helper functions ======//


/*
 
 Prunes a LineageTreeNode<T,U> transmission tree with 'root' as root node,
 using a vector 'sampledLngs' of lineages.
 
 Returns the reduced transmission tree.
 
 */

template <typename T, typename U>
LineageTreeNode<T,U>* eliminateRedundantNodes( LineageTreeNode<T,U>* root, const std::vector<T>& sampledLngs ) {
    
    // find leaves
    std::unordered_set<LineageTreeNode<T,U>*> leaves = {} ;
    leaves.reserve( sampledLngs.size() ) ;

    findLeaves( leaves, root ) ; // find leaves
    
    // starting from each leaf, recurse back to root and remove intermediate leaves
    for ( auto leaf : leaves )
        removeRedundantNodeMerge( leaf->parent, sampledLngs );
        
    return findRoot( *( leaves.begin() ) ) ; // return tree as the root of the first sampled lineage.
    
} ;

/*
 
 Returns the root of the tree starting from 'lngNode'.
 
 N.B. This could be 'lngNode' itself if it is a root node.
 
 */

template <typename T, typename U>
LineageTreeNode<T,U>* findRoot( LineageTreeNode<T,U>* lngNode ) {
    
    if ( lngNode->parent == nullptr )
        return lngNode ;
    else
        return findRoot( lngNode->parent ) ;
    
}

/*
 
 Fills the vector 'leaves' with EXTREMAL lineages descending from 'node'.
 
 If 'node' is a root node, it fills 'leaves' with all EXTREMAL lineages.
 
 */

template <typename T, typename U>
void findLeaves( std::unordered_set<LineageTreeNode<T,U>*>& leaves, LineageTreeNode<T,U>* node ) {
    
    if ( node->children.size() == 0 ) {
        
        //if ( node->parent != nullptr ) // root node can not be leaf
        leaves.insert( node ) ; // root node can be a leaf if it is the only node
        
    }
    else {
        
        for ( auto child : node->children )
            findLeaves( leaves, child ) ;
        
    }
    
}

/*
 
 Helper function that removes 'midNode' in a transmission tree if it is
 found to be redundant given a vector 'sampledLngs' with all sampled lineages.
 
 See also 'eliminateRedundantNodes'.
 
 */

template <typename T, typename U>
void removeRedundantNodeMerge( LineageTreeNode<T,U>* midNode, const std::vector<T>& sampledLngs ) {
    
    if ( !midNode )
        return ;
    
    if ( midNode->children.size() == 1 ) { // may merge if has one child only
    
        // check if node is sampled
        bool isSampled = false ;
        for ( auto& lng : sampledLngs ) {
            
            if ( lng == midNode->lng ) {
            
                isSampled = true ;
                break ;
            
            }
        
        }
        if ( !isSampled ) {  // if not sampled, remove
            
            if ( midNode->parent == nullptr ) { // if parent is root
                
                // promote only child to root
                midNode->children[0]->parent = nullptr ;
                midNode->children[0]->tBranchParent = midNode->children[0]->t ; // this should be OK given that tBranchParent is irrelevant for roots
                
            }
            else {
                
                midNode->children[0]->parent = midNode->parent ;
                midNode->parent->eraseChild( midNode ) ;
                midNode->parent->children.push_back( midNode->children[0] ) ;

                //midNode->parent->children_branching_times[ midNode->children[0]->lng ] = midNode->parent->children_branching_times[ midNode->lng ] ;
                //midNode->parent->children_branching_times.erase( midNode->lng ) ;
                midNode->children[0]->tBranchParent = midNode->tBranchParent ;
                removeRedundantNodeMerge( midNode->parent, sampledLngs ) ;
                
                
            }
            
            delete midNode ;
            
        }
        else { // recurse up
            
            if ( midNode->parent != nullptr ) // recurse up until root
                removeRedundantNodeMerge( midNode->parent, sampledLngs ) ;
            
        }
        
    }
    else { // recurse up
        
        if ( midNode->parent != nullptr ) // recurse up until root
            removeRedundantNodeMerge( midNode->parent, sampledLngs ) ;
        
    }
    
}


//====== PhyloNode ======//

/*
 
 Phylogenetic tree building block
 
 Holds info about:
 - internal nodes
 - leaf nodes
 - node metadata
 - node timing
 - branch length
 
*/

template <typename T, typename U>
struct PhyloNode {
    
    PhyloNode( const T& lng, PhyloNode* parent = nullptr  ): lng( lng ), leftChild( nullptr ), rightChild( nullptr ), parent( parent ), depth( 0 ), depthChild( 0 ), depthAttachSampledNode( -1 ), t( 0 ), dt( 0 ), locSample( "NA" ) {} ;
    PhyloNode* leftChild ;
    PhyloNode* rightChild ;
    PhyloNode* parent ;
    uint depth ; // initialised to 0, used to track depth of internal nodes
    uint depthChild ; // initialised to 0, used to track child index
    uint depthAttachSampledNode ; // initialised to -1, used to track where sampled ancestors must be placed
    double t ; // node time (infection time if internal, sampling time if leaf)
    double dt ; // branch length (wrt parent)
    std::string locSample ; // sampling location (if any; default is NA)
    T lng ;  // lineage identity
    U data ; // extra data
    friend std::ostream& operator<<(std::ostream& os, const PhyloNode<T,U>*& dt);
        
} ;



/*
 
 Frees memory initially allocated to a tree composed of PhyloNode<T,U>.
 
 */

template <typename T, typename U>
void deletePhyloNodeTree( PhyloNode<T,U>* root ) {
    
    if ( root->leftChild != nullptr )
        deletePhyloNodeTree( root->leftChild ) ;
    
    if ( root->rightChild != nullptr )
        deletePhyloNodeTree( root->rightChild ) ;
    
    delete root ;
    
}





/*
 
 Returns a phylogenetic tree from a reduced transmission tree.
 
 'node' can be any node in the reduced transmission tree, but
 the complete phylogenetic tree is returned by calling the
 function on the root node of the reduced transmission tree.
 
 'getAncestralTree' will be called recursively on downstream
 nodes.
 
 In the resulting tree, sampled lineages appear as leaf nodes,
 while internal nodes correspond to past infection events.
 However, it may also happen that some ancestral lineages are also leaves.
 This might happen if we manage to sample the parent of a lineage.
 
 */


template <typename T, typename U>
PhyloNode<T,U>* getAncestralTree( LineageTreeNode<T,U>* node, PhyloNode<T,U>* phyloParent = nullptr ) {

    if ( node == nullptr )
        return nullptr ;
        
    bool isChild   = ( phyloParent == nullptr ) ? false : true ; // true : is child of some other node; false : is root node
    bool isSampled = node->sampled ;
    
    // create phylo node (internal or tip)
    PhyloNode<T,U>* newNode = new PhyloNode<T,U>( node->lng, phyloParent ) ;
    
    // add data
    newNode->data = node->data ;

    // calculate node depth and set depthChild
    newNode->depth      = 0 ;
    newNode->depthChild = 0 ;

    if ( isChild and ( newNode->lng == phyloParent->lng ) ) { // this child is part of a chain of transmission events from the same source
        newNode->depth = phyloParent->depth + 1 ; // is incremented
        newNode->depthChild = phyloParent->depthChild ; // must be incremented only when used
        newNode->depthAttachSampledNode = phyloParent->depthAttachSampledNode ;

    }
    
    // sort children chronologically (according to tBranchParent)
    auto& children_sorted = node->children ;
    uint nChildren = static_cast<uint>( children_sorted.size() ) ;
    
    if ( ( nChildren > 1 ) and ( newNode->depth == 0 ) ) {
        
        std::sort( children_sorted.begin(), children_sorted.end(), [&]( LineageTreeNode<T,U>*& node1, LineageTreeNode<T,U>* node2 ) {
                return node1->tBranchParent < node2->tBranchParent ;  // sort by time from the map
            } ) ;
    
    }
    
    if ( nChildren == 0 ) { // sampled leaf with no children
        
        assert( isSampled ) ; // must be a leaf
        newNode->t = node->tSample ;
        if ( phyloParent == nullptr ) newNode->dt = 0 ;
        else {
            newNode->dt = newNode->t - phyloParent->t ;
            assert( newNode->t >= phyloParent->t ) ;
        }
        newNode->locSample = node->locSample ; // only for sampled nodes
        
    }
    else { // has some children: this leaf is also an ancestor for some other nodes (sampled ancestor)
        
        if ( isSampled ) { // node was sampled
            
            double tSample = node->tSample ;
            
            // calculate (ONLY ONCE when depth = 0) where node should be placed
            
            if ( newNode->depth == 0 ) { // calculate only once
                
                newNode->depthAttachSampledNode = 0 ; // position where sampled lineage should be placed
                for ( auto& child : children_sorted ) {
                    
                    if ( tSample < child->tBranchParent ) // stop before sampling time exceeds branching time
                        break ;
                    else
                        newNode->depthAttachSampledNode++ ;

                }
                
                //if ( newNode->depthAttachSampledNode == nChildren ) // node has been sampled after all children (notice that nChildren > 0 )
                //    newNode->depthAttachSampledNode -= 1 ; // Explanation : it does not matter if the ancestor was sampled after the last or before-to-last sampling time; the attachment order is always nChildren - 1 (last intenr
                
            }
                        
            if ( newNode->depthAttachSampledNode < nChildren ) { // node is sampled before some children are created
                
                if ( newNode->depth == newNode->depthAttachSampledNode ) { // attach sampled node
                    
                    newNode->t = tSample ;
                    if ( phyloParent == nullptr ) {
                        newNode->dt = 0. ;
                    }
                    else {
                        newNode->dt = newNode->t - phyloParent->t ;
                        assert( newNode->t >= phyloParent->t ) ;
                    }
                    
                    PhyloNode<T,U>* sampledNode = new PhyloNode<T,U>( node->lng, newNode ) ;
                    sampledNode->t = node->tSample ;
                    sampledNode->dt = 0. ; // Sampled ancestor has 0 branch length..
                    sampledNode->depth = newNode->depth + 1 ;
                    sampledNode->data = node->data ;
                    sampledNode->locSample = node->locSample ; // only for sampled nodes
                    sampledNode->leftChild = nullptr ;
                    sampledNode->rightChild = nullptr ;
                    newNode->rightChild = sampledNode ;
                    
                    if ( newNode->depthChild == nChildren - 1  ) {
                        auto& child = children_sorted[ newNode->depthChild ] ;
                        newNode->depthChild += 1 ;
                        newNode->leftChild = getAncestralTree( child, newNode ) ;
                    }
                    else {
                        newNode->leftChild = getAncestralTree( node, newNode ) ;
                    }

                
                }
                else { // attach child + internal node (or two nodes)
                    
                    auto& child = children_sorted[ newNode->depthChild ] ;
                    newNode->t = child->tBranchParent ;
                    newNode->depthChild += 1 ;
                    
                    if ( phyloParent == nullptr ) {
                        newNode->dt = 0 ;
                    }
                    else {
                        newNode->dt = newNode->t - phyloParent->t ;
                        assert( newNode->t >= phyloParent->t ) ;
                    }
                    
                    newNode->leftChild = getAncestralTree( child, newNode ) ;
                    
                    
                    if (  newNode->depth == nChildren - 1 ) {
                        
                        auto& child2 = children_sorted[ newNode->depthChild ] ;
                        newNode->rightChild = getAncestralTree( child2,  newNode ) ;
                        
                    }
                    else {
                        
                        newNode->rightChild = getAncestralTree( node,  newNode ) ;
                        
                    }
                    
                }
                
            }
            else { // node is sampled after all children are created
                
                newNode->t = children_sorted[ newNode->depth ]->tBranchParent ;
                if ( phyloParent == nullptr ) {
                    newNode->dt = 0 ;
                }
                else {
                    newNode->dt = newNode->t - phyloParent->t ;
                    assert( newNode->t >= phyloParent->t ) ;
                }
                
                if ( newNode->depth < nChildren - 1 ) { // attach child + internal node
                    
                    auto& child = children_sorted[ newNode->depth ] ;
                    newNode->leftChild  = getAncestralTree( child, newNode ) ;
                    newNode->rightChild = getAncestralTree( node, newNode ) ;
                       
                }
                else { // add last child and sampled node
 
                    auto& child = children_sorted[ newNode->depth ] ;
                    
                    // create sampled ancestor node (a leaf)
                    PhyloNode<T,U>* sampledNode = new PhyloNode<T,U>( node->lng, newNode ) ;
                    sampledNode->t = node->tSample ;
                    sampledNode->dt = node->tSample - children_sorted[ newNode->depth ]->tBranchParent ;
                    sampledNode->depth = newNode->depth + 1 ;
                    sampledNode->data = node->data ;
                    sampledNode->locSample = node->locSample ; // only for sampled nodes
                    sampledNode->leftChild = nullptr ;
                    sampledNode->rightChild = nullptr ;

                    newNode->leftChild  = getAncestralTree( child, newNode ) ;
                    newNode->rightChild  = sampledNode ;
                    
                }
                
            }
            
        }
        else { // node is not sampled (easy case)
            
            assert( nChildren >= 2 ) ; // node must have at least two children (otherwise it would have been removed already)

            bool isLastEvent = ( newNode->depth == nChildren ) ? true : false ;
            assert( !isLastEvent ) ; // we should never get to this point
            
            newNode->t = children_sorted[ newNode->depth ]->tBranchParent ;
            if ( phyloParent == nullptr )
                newNode->dt = 0 ;
            else {
                newNode->dt = newNode->t - phyloParent->t ;
                assert( newNode->t >= phyloParent->t ) ;
            }

            if ( newNode->depth < nChildren - 2 ) { // recurse down
                
                auto& child = children_sorted[ newNode->depth ] ;
                //newNode->depthChild += 1 ;
                newNode->leftChild  = getAncestralTree( child, newNode ) ;
                newNode->rightChild = getAncestralTree( node,  newNode ) ;
                
            }
            else { // stop recursion: last cherry in the tree
                
                auto& child1 = children_sorted[ newNode->depth ] ;
                auto& child2 = children_sorted[ newNode->depth + 1 ] ;

                newNode->leftChild  = getAncestralTree( child1, newNode ) ;
                newNode->rightChild = getAncestralTree( child2, newNode ) ;
                
            }
            
            
        }
        
        
    }
    
    return newNode ;
    
}

/*
 
 Converts 'lng' to string.
 N.B. MAY REQUIRE OVERLOADING '<<' OPERATOR IF 'T' IS NOT A BASE TYPE
 
 */
template <typename T>
std::string lng2string( T lng ) {
    
    std::ostringstream strm;
    strm << lng ;
    return strm.str() ;
    
}

/*
 
 Converts data from 'node' to string.
 N.B. REQUIREs OVERLOADING '<<' OPERATOR
 
 */
template <typename T, typename U>
std::string data2string( PhyloNode<T,U>* node ) {
    
    std::ostringstream strm;
    strm << node->data ;
    return strm.str() ;
    
}


/*
 Yields a phylogenetic tree in NHX format
 */

template <typename T, typename U>
std::string getNHX( PhyloNode<T,U>* root ) {
    
    std::string nhx ; // holds result
    PhyloNode2NHX( nhx, root ) ; // begin recursion
    nhx += ";" ; // closing character
    
    return nhx ;
    
}

/*
Recursive function to build a phylogenetic tree in NHX format
 */

template <typename T, typename U>
void PhyloNode2NHX( std::string& nhx, PhyloNode<T,U>* node ) {
    
    bool isLeaf = ( node->leftChild == nullptr ) ? true : false ;
    
    if ( isLeaf ) {
        
        nhx += lng2string( node->lng ) + ":" + std::to_string( node->dt ) ;
        nhx += "[&&NHX:" + data2string( node ) + ":" + std::to_string( node->t ) + "]";
        
    }
    else { // manage branching event
        
        nhx += "(" ;
        PhyloNode2NHX( nhx, node->leftChild ) ;
        nhx += "," ;
        PhyloNode2NHX( nhx, node->rightChild ) ;
        nhx += ")" ;
        nhx += lng2string( node->lng ) + "-" + std::to_string( node->depth ) ;
        nhx += ":" + std::to_string( node->dt ) ;
        nhx += "[&&NHX:" + data2string( node ) + ":" + std::to_string( node->t ) + "]";


    }
    
}

/*
 Yields a phylogenetic tree in Newick format
 */

template <typename T, typename U>
std::string getSimpleNewick( PhyloNode<T,U>* root ) {
    
    std::string nhx ; // holds result
    PhyloNode2Newick( nhx, root ) ; // begin recursion
    nhx += ";" ; // closing character
    
    return nhx ;
    
}

/*
Recursive function to build a phylogenetic tree in NHX format
 */

template <typename T, typename U>
void PhyloNode2Newick( std::string& nhx, PhyloNode<T,U>* node ) {
    
    bool isLeaf = ( node->leftChild == nullptr ) ? true : false ;
    
    if ( isLeaf ) {
        
        nhx += lng2string( node->lng ) + ":" + std::to_string( node->dt ) ;
        
    }
    else { // manage branching event
        
        nhx += "(" ;
        PhyloNode2Newick( nhx, node->leftChild ) ;
        nhx += "," ;
        PhyloNode2Newick( nhx, node->rightChild ) ;
        nhx += ")" ;
        nhx += lng2string( node->lng ) + "-" + std::to_string( node->depth ) ;
        nhx += ":" + std::to_string( node->dt ) ;


    }
    
}

#endif /* tree_h */
