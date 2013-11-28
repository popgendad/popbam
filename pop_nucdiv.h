/** \file pop_nucdiv.h
 *  \brief Header for the pop_nucdiv.cpp file
 *  \author Daniel Garrigan
 *  \version 0.4
*/

#include "pop_base.h"

///
/// Definitions
///

/*! \def EPSILON
 *  \brief Value of epsilon for testing whether double is zero
 */
#define EPSILON 1e-08

//
// Define data structures
//

/*!
 * \class nucdivData
 * \brief A derived class for passing parameters and data to the nucdiv function
 */
class nucdivData: public popbamData
{
	public:
		// constructor
		nucdivData(const popbamOptions&);

		// destructor
		~nucdivData(void);

		// member public variables
		unsigned int *pop_cov;                  //!< Boolean for population coverage
		unsigned int **ncov;                    //!< Sample size per population per segregating site
		unsigned long *ns_within;               //!< Number of aligned sites with each population
		unsigned long *ns_between;              //!< Number of aligned sites between each pair of populations
		int *num_snps;                          //!< Number of SNPs in a given window
		double minPop;                          //!< Minimum proportion of samples present

		// member public functions
		int calcNucdiv(void);
		int allocNucdiv(void);
		int printNucdiv(const std::string);

	private:
		// member private variables
		double minSites;                        //!< User-specified minimum proportion of aligned sites to perform analysis
		double *piw;                            //!< Array of within-population nucleotide diversity
		double *pib;                            //!< Array of between-population Dxy values
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
template unsigned long long* callBase<nucdivData>(nucdivData *t, int n, const bam_pileup1_t *pl);

/*!
 * \fn int makeNucdiv(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
 * \brief Runs the nucleotide diversity calculations
 * \param tid Chromosome identifier
 * \param pos Genomic position
 * \param n The read depth
 * \param pl A pointer to the alignment covering a single position
 * \param data A pointer to the user-passed data
 */
int makeNucdiv(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data);

void usageNucdiv(const std::string);
