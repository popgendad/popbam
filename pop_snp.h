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
        snpData(const popbamOptions&);

        // destructor
        ~snpData(void);

        // member public variables
        hData_t hap;                            //!< Structure to hold haplotype data
        unsigned int *pop_cov;                  //!< Boolean for population coverage
        unsigned int **ncov;                    //!< Sample size per population per segregating site
        unsigned long long **pop_sample_mask;   //!< Bit mask for samples covered from a specific population
        int output;                             //!< User-specified output mode
        std::string outgroup;                   //!< Sample name of outgroup to use
        int outidx;                             //!< Index of outgroup sequence
        double minPop;                          //!< Minimum proportion of samples present

        // member public functions
        int allocSNP(void);
        int printMSHeader(long);
        int print_SNP(const std::string);

    private:
        // member private functions
        int printSNP(const std::string);
        int printSweep(const std::string);
        int printMS(const std::string);
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
 * \fn int makeSNP(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
 * \brief Runs the SNP analysis
 * \param tid Chromosome identifier
 * \param pos Genomic position
 * \param n The read depth
 * \param pl A pointer to the alignment covering a single position
 * \param data A pointer to the user-passed data
 */
int makeSNP(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data);

void usageSNP(const std::string);

typedef int(snpData::*snp_func)(const std::string);
