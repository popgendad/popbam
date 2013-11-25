/** \file pop_snp.h
 *  \brief Header for the pop_snp.cpp file
 *  \author Daniel Garrigan
 *  \version 0.4
*/

#include "pop_base.h"

char bam_nt16_rev_table[16] = {'=','A','C','M','G','R','S','V','T','W','Y','H','K','D','B','N'};

/*!
 * \class snpData
 * \brief A derived class for passing parameters and data to the snp function
 */
class snpData: public popbamData
{
	public:
		// constructor
		snpData();

		// destructor
		~snpData() {}

		// member public variables
		hData_t hap;                            //!< Structure to hold haplotype data
		unsigned int win_size;                  //!< User-specified window size in kilobases
		unsigned int *pop_cov;                  //!< Boolean for population coverage
		unsigned int **ncov;                    //!< Sample size per population per segregating site
		unsigned long long **pop_sample_mask;   //!< Bit mask for samples covered from a specific population
		int output;                             //!< User-specified output mode
		std::string outgroup;                   //!< Sample name of outgroup to use
		int outidx;                             //!< Index of outgroup sequence
		double min_pop;                         //!< Minimum proportion of samples present

		// member public functions
		std::string parseCommandLine(int, char**);
		void init_snp(void);
		int printMSHeader(long);
		void print_snp(int);
		void destroy_snp(void);

	private:
		// member private functions
		int printSNP(int);
		int printSweep(int);
		int printMS(int);
		void printUsage(std::string);
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
template unsigned long long* callBase<snpData>(snpData *t, int n, const bam_pileup1_t *pl);

/*!
 * \fn int make_snp(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
 * \brief Runs the SNP analysis
 * \param tid Chromosome identifier
 * \param pos Genomic position
 * \param n The read depth
 * \param pl A pointer to the alignment covering a single position
 * \param data A pointer to the user-passed data
 */
int make_snp(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data);

typedef int(snpData::*snp_func)(int);
