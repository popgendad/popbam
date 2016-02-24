/** \file pop_tree.h
 *  \brief Header for the pop_tree.cpp file
 *  \author Daniel Garrigan
 *  \version 0.4
 * Much of the code for the NJ algorithm heavily borrows from the PHYLIP
 * package (v 3.6) written by Mary Kuhner, Jon Yamato, Joseph Felsenstein,
 * Akiko Fuseki, Sean Lamont, and Andrew Keefe at the University of Washington
*/

#ifndef POP_TREE_H
#define POP_TREE_H

//
// Define data structures
//

/*!
 * \struct node
 * \brief A structure to represent an individual node in a bifurcating tree
 */

typedef struct _node
{
    struct _node *next;               //!< Pointer to the next node in the list
    struct _node *back;               //!< Pointer to the previous node in the list
    bool tip;                         //!< Is the node a tip?
    int index;                        //!< Index number of the node
    double v;                         //!< Branch length above the node
} node;

/*!
 * \var typedef node **ptarray
 * \brief A pointer to a node pointer
 */

typedef node **ptarray;

/*!
 * \struct tree
 * \brief A structure to represent a bifurcating tree
 */

typedef struct _tree
{
    node *start;                      //!< Pointer to the start node of the tree
    ptarray nodep;                    //!< Pointer to the linked list of nodes in the tree
} tree;

/*!
 * \class treeData
 * \brief A derived class for passing parameters and data to the tree function
 */

class treeData: public popbamData
{
public:
    // constructor
    treeData(const popbamOptions&);

    // destructor
    ~treeData(void);

    // member public variables
    hData_t hap;                   //!< Structure to hold haplotype data (public)
    uint64_t *pop_sample_mask;     //!< Bit mask for samples covered from a specific population
    char *refid;                   //!< Pointer to string that holds the reference sequence name
    std::string dist;              //!< Pointer to the name of the desired distance metric  (-d switch)
    uint16_t **diff_matrix;        //!< Array of pairwise sequence differences
    double **dist_matrix;          //!< Array of divergence calculations
    int npops;                     //!< Number of populations represented in BAM
    int minSites;                  //!< User-specified minimum number of aligned sites to perform analysis
    int ntaxa;                     //!< Total number of tips in the tree
    int *enterorder;               //!< Array containing the input order of OTUs for the NJ algorithm

    // member public functions
    int makeNJ(const std::string);
    int calcDistMatrix(void);
    int allocTree(void);
    void joinTree(tree, node**);
    void printTree(node*, node*);
    void hookup(node*, node*);
    void setupTree(tree*);
    void initTree(ptarray *);
    void freeTree(ptarray *);
};

///
/// Function prototypes
///

extern int mainTree(int argc, char *argv[]);

#endif
