/** \file pop_tree.h
 *  \brief Header for the pop_tree.cpp file
 *  \author Daniel Garrigan
 *  \version 0.4
 * Much of the code for the NJ algorithm heavily borrows from the PHYLIP
 * package (v 3.6) written by Mary Kuhner, Jon Yamato, Joseph Felsenstein,
 * Akiko Fuseki, Sean Lamont, and Andrew Keefe at the University of Washington
*/

#include "pop_base.h"

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
		hData_t hap;                            //!< Structure to hold haplotype data (public)
		unsigned long long *pop_sample_mask;    //!< Bit mask for samples covered from a specific population
		char *refid;                            //!< Pointer to string that holds the reference sequence name
		std::string dist;                       //!< Pointer to the name of the desired distance metric	(-d switch)
		unsigned short **diff_matrix;           //!< Array of pairwise sequence differences
		double **dist_matrix;                   //!< Array of divergence calculations
		int minSites;                           //!< User-specified minimum number of aligned sites to perform analysis
		int ntaxa;                              //!< Total number of tips in the tree
		int *enterorder;                        //!< Array containing the input order of OTUs for the NJ algorithm

		// member public functions
		int makeNJ(const std::string);
		int calc_dist_matrix(void);
		int allocTree(void);
		void printUsage(const std::string);
		void join_tree(tree, node**);
		void print_tree(node*, node*);
		void hookup(node*, node*);
		void setup_tree(tree*);
		void tree_init(ptarray *);
		void free_tree(ptarray *);
};

///
/// Function prototypes
///

/*!
* \fn unsigned long long *callBase(bam_sample_t *sm, errmod_t *em, int n, const bam_pileup1_t *pl)
* \brief Calls the base from the pileup at each position
* \param sm     Pointer to the sample data structure
* \param em     Pointer to the error model structure
* \param n      The number of reads in the pileup
* \param pl     Pointer to the pileup
* \return       Pointer to the consensus base call information for the individuals
*/
template unsigned long long* callBase<treeData>(treeData *t, int n, const bam_pileup1_t *pl);

/*!
 * \fn int make_tree(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
 * \brief Runs the neighbor-joining tree construction procedure
 * \param tid Chromosome identifier
 * \param pos Genomic position
 * \param n The read depth
 * \param pl A pointer to the alignment covering a single position
 * \param data A pointer to the user-passed data
 */
int makeTree(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data);

void calc_diff_matrix(treeData*);
