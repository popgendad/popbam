/** \file pop_diverge.h
 *  \brief Header for the pop_diverge.cpp file
 *  \author Daniel Garrigan
 *  \version 0.4
*/

#include "pop_base.h"
//
// Define data structures
//

/*!
 * \class divergeData
 * \brief A derived class for passing parameters and data to the diverge function
 */
class divergeData: public popbamData
{
	public:
		// constructor
		divergeData(const popbamOptions&);

		// destructor
		~divergeData(void);

		// member public variables
		int output;                             //!< Analysis output option
		std::string outgroup;                   //!< Sample name of outgroup to use
		int outidx;                             //!< Index of outgroup sequence
		hData_t hap;                            //!< Structure to hold haplotype data
		unsigned long long *pop_sample_mask;    //!< Bit mask for samples covered from a specific population
		std::string dist;                       //!< Pointer to the name of the desired distance metric	(-d switch)
		unsigned short *min_pop_n;              //!< Minimum sample size per population
		int *num_snps;                          //!< Number of SNPs in a given window

		// member public functions
		int calcDiverge(void);
		int allocDiverge(void);
		int setMinPop_n(void);
		int printDiverge(const std::string);

	private:
		// member private variables
		int minSites;                           //!< User-specified minimum number of aligned sites to perform analysis
		unsigned short *pop_div;                //!< Array of mean population divergence calculations
		unsigned short *ind_div;                //!< Array of individual divergence calculations
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
template unsigned long long* callBase<divergeData>(divergeData *t, int n, const bam_pileup1_t *pl);

/*!
 * \fn int makeDiverge(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
 * \brief Calculate divergence with reference genome sequence
 * \param tid Chromosome identifier
 * \param pos Genomic position
 * \param n The read depth
 * \param pl A pointer to the alignment covering a single position
 * \param data A pointer to the user-passed data
 */
int makeDiverge(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data);

void usageDiverge(const std::string);
