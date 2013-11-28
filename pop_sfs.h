/** \file pop_sfs.h
 *  \brief Header for the pop_sfs.cpp file
 *  \author Daniel Garrigan
 *  \version 0.4
*/

#include "pop_base.h"

///
/// Additional include headers
///
#include <limits>

///
/// Definitions
///

//
// Define data structures
//
/*!
 * \class sfsData
 * \brief A derived class for passing parameters and data to the sfs function
 */
class sfsData: public popbamData
{
	public:
		// constructor
		sfsData(const popbamOptions&);

		// destructor
		~sfsData(void);

		// member variables
		double minSites;                        //!< User-specified minimum proportion of aligned sites to perform analysis
		unsigned long *ns;                      //!< Number of aligned sites within each population
		unsigned int **ncov;                    //!< Sample size per population per segregating site
		int *num_snps;                          //!< Number of SNPs in a given window
		unsigned int *pop_cov;                  //!< Boolean for population coverage
		double minPop;                          //!< Minimum proportion of samples present
		std::string outgroup;                   //!< Sample name of outgroup to use
		int outidx;                             //!< Index of outgroup sequence
		int outpop;                             //!< Population of outgroup sequence
		double **dw;                            //!< Matrix of weights for Tajima's D calculation
		double **hw;                            //!< Matrix of weights for Fay and Wu's H calculation
		double *a1;                             //!< Constants for Tajima's D calculation
		double *a2;                             //!< Constants for Tajima's D calculation
		double *e1;                             //!< Constants for Tajima's D calculation
		double *e2;                             //!< Constants for Tajima's D calculation
		double *td;                             //!< Pointer to the array of Tajima's D statistics
		double *fwh;                            //!< Pointer to the array of standardized Fay and Wu's H statistics

		// member functions
		int allocSFS(void);
		int printSFS(const std::string);
		int assignOutpop(void);
		int calc_dw(void);
		int calc_hw(void);
		int calcSFS(void);
		int calc_a1(void);
		int calc_a2(void);
		int calc_e1(void);
		int calc_e2(void);
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
template unsigned long long* callBase<sfsData>(sfsData *t, int n, const bam_pileup1_t *pl);

/*!
 * \fn int makeSFS(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
 * \brief Runs the site frequency spectrum analysis
 * \param tid Chromosome identifier
 * \param pos Genomic position
 * \param n The read depth
 * \param pl A pointer to the alignment covering a single position
 * \param data A pointer to the user-passed data
 */
int makeSFS(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data);

void usageSFS(const std::string);
