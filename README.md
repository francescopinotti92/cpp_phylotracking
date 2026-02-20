## Intro

This repository is meant to illustrate how to generate phylogenetic trees during epidemic simulations in C++.

The header-only library `src/tree.hpp` provides data structures and functions to track transmission chains within simulations and output phylogenetic trees where tips represent sampled infections and branch lengths are expressed in time units.

The tracking algorithm uses a dynamical pruning strategy to reduce memory requirements. Briefly, stored transmission chains are pruned continually following the extinction of lineages (i.e. infections) that can not be sampled or that will never appear in phylogenetic trees as internal nodes. Similarly, paths along transmission chains are shortened as much as possible to retain only the informations necessary to draw ancestral relationships between extant and sampled lineages.

`tree.hpp` currently supports newick and New Hampshire X (NHX) tree formats. Both newick and NHX trees appears as strings of text containing ancestral relationships, tips/internal node names and branch lengths. The NHX format adds the possibility to include additional node metadata.

## How to use?

We exemplify the use of `tree.hpp` in the case of Birth & Death (BD) model. The latter simulates a branching process where individual lineages reproduce with rate $\beta$ and die with rate $\alpha$. Dying lineages are sampled with probability $\rho$ and will therefore appear in the final tree.

The BD model is a very simple transmission model that we only use to illustrate how to use `tree.hpp`; if you want to use this library alongside your own simulator you will need to follow the same steps.

First, `tree.hpp` introduces the `LineageTreeNode<T,U,H>` class as the basic transmission chain building block. Transmission chains are organised in a `LineageTree<T,U,H>` object

 The template argument `T,U,H` denote:
- The type associated with a lineage identifier (`T`).
- The type associated with a lineage metadate (`U`).
- The hash function associated with `T` (`H`). There is no need to specify `H` if `T` is a basic type like `int`. In that case `H` simply defaults to `std::Hash<T>`.

In the BD example every lineage has a unique integer identifier, hence `T=int`. We are not interested in metadata either, so we simply set `U=int` and ignore it. We then endow our `Simulator` class with a `LineageTree<int,int>` instance (`tree_mngr`).

Now, `Simulator` is in charge of notifying `tree_mngr` whenever a new lineage is created, and whenever extant lineage die or are sampled.

Regarding lineage creation we must distinguish between introduction and transmission events. Introductions do not involve a parent-child relationship and add a new transmission chain in `tree_mngr`, with the introduced lineage as its root.

To record an introduction, you need to call the following function:

```cpp
tree_mngr->addExtantLineageExternal( time, lng_ID, event_data );
```

where `time` denotes the time of this event `lng_ID` is the lineage identifier and `event_data` is the eventual metadata.

To record a transmission event where `lng_child_ID` was generated from `lng_parent_ID`, call:

```cpp
tree_mngr->addExtantLineage( time, lng_child_ID, event_data, lng_parent_ID );
```

Whenever a lineage dies, call:

```cpp
tree_mngr->removeExtantLineage( lng_ID );
```

If you sample a lineage at time `time` you need to mark it as sampled by calling:

```cpp
tree_mngr->sampleExtantLineage( lng_ID, time );
```

### Collecting the tree

The next instructions show how to get a phylogenetic tree from the transmission chains. Importantly, the tips of the tree correspond to sampled lineages.

First, we get `tree_mngr` from `simulator`:

```cpp
LineageTree<int,int>* tree_mngr = simulator.get_tree();
```

Second, we prune the transmission tree from unnecessary transmission events and lineages (unnecessary with regard to the phylogenetic tree of sampled lineages):    

```cpp
LineageTreeNode<int,int>* rtree = tree_mngr->subSampleTree()[0];
```

Please note that `tree_mngr->subSampleTree()` yields a vector of reduced transmission trees, whose size corresponds to the number of independent transmission chains sampled lineages belong to. This is not an issue in the BD example since we know that the epidemic starts from a single seed and hence a single transmission chain. In general, however, users may end up with more than transmission chain and hence multiple independent trees.
Given a reduced transmission chain `rtree`, we proceed to extract the phylogenetic tree (`atree`):

```cpp
PhyloNode<int,int>* atree = getAncestralTree( rtree );
```
`atree` is a `PhyloNode<int,int>*` object (`PhyloNode<T,U>*` in general), which is the building block for phylogenetic trees. `atree` is actually the root node of the tree!

Finally, we generate a newick/NHX string from the tree:

```cpp
std::string nwk = getSimpleNewick( atree );
```

## Using custom lineage identifiers

What if lineage identifiers are not basic C++ types?
Let us consider an identifier specified by two `int` as in the `CumstomID` class:

```cpp
struct CumstomID {
    CustomID( int a, int b ): A(a), B(b) {};
    int A;
    int B;
} ;
```
We need to be able to tell whether two lineages are identical. To do so, we endow `CumstomID` with its own `==` operator:

```cpp
struct CumstomID {
    CustomID( int a, int b ): A(a), B(b) {};
    int A;
    int B;
    bool operator==(const CumstomID& other) const { return A == other.A && B == other.B ; }
} ;
```
Then we need to specify a custom hash function. The reason is that C++ maps require hash functions to manage their keys.

```cpp
namespace std {
    template<>
    struct hash<CustomID> {
        std::size_t operator()(const CustomID& ID) const noexcept {
            std::size_t h1 = std::hash<int>()(ID.A);
            std::size_t h2 = std::hash<int>()(ID.B);
            return h1 ^ (h2 << 1);
        }
    };
}
```

Finally, we need to tell C++ how to turn `CumstomID` into a `string`. This can be done by overloading the `<<` operator:

```cpp
std::ostream& operator<<( std::ostream& os, CustomID& ID ) {
   os << ID.A << "-" << ID.B ;
   return os ;
}
```

Alternatively, we can overload the `lng2string` function:

```cpp
std::string lng2string( const CustomID& ID ) {
    return std::to_string( ID.A ) + "-" + std::to_string( ID.B ) ;
}
```

The corresponding string will look like `A-B`.

## Using custom metadata

What if lineage metadata are not basic C++ types?

Let us consider some custom metadata class `CumstomData` :

```cpp
struct CumstomData {
    CumstomData( int x ): X(x) {};
    int X;
} ;
```
We can then overload the `<<` operator or the `data2string` function:

```cpp
std::string data2string( PhyloNode<CustomID,CustomData>* node ) {
    
    CustomData& data = node->data ;
    return "X=" + std::to_string( CustomData.X ) ;
    
}
```

## A note on sampled ancestors

In some applications it might happen that we record sampled ancestors. These are internal nodes that are also leaves, and might appear if, e.g., we sample an infection but do not remove it, and subsequently sample some of its descendants (This can not happen in the BD example since lineages are sampled only upon removal).

In the resulting phylogenetic tree, sampled ancestors appear as leaves with 0 branch length. This structure preserves the binary tree structure. However, note that other software may represent sampled ancestor differently, e.g. as additional internal nodes with a single child.

# Running the code

We wrapped the code to generate BD trees in a python module (`pysimBD`). To compile the code into a module, open the terminal and move to this folder, then type `make`. Please make sure to install `pybind11` before that and modify the `makefile` variables `CXX`, `CXXFLAGS`, `INC` and `EXT` to match the specifics of your system.
