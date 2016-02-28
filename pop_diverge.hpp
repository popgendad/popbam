/** \file pop_diverge.hpp
 *  \brief Header for the pop_diverge.cpp file
 *  \author Daniel Garrigan
 *  \version 0.5
*/

#ifndef POP_DIVERGE_HPP
#define POP_DIVERGE_HPP

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
    divergeData(void);

    // destructor
    ~divergeData(void);

    // member public variables
    int npops;                      //!< Number of populations represented in BAM
    int outidx;                     //!< Index of outgroup sequence
    int *num_snps;                  //!< Number of SNPs in a given window
    uint16_t *min_pop_n;            //!< Minimum sample size per population
    uint64_t *pop_sample_mask;      //!< Bit mask for samples covered from a specific population
    hData_t hap;                    //!< Structure to hold haplotype data


    // member public functions
    int calc_diverge(const pop_diverge_parser *param);
    int alloc_diverge(void);
    int set_min_pop_n(void);
    int print_diverge(const pop_diverge_parser *param, const char *scaffold, int win_size);

private:
    // member private variables
    int minSites;              //!< User-specified minimum number of aligned sites to perform analysis
    uint16_t *pop_div;         //!< Array of mean population divergence calculations
    uint16_t *ind_div;         //!< Array of individual divergence calculations
};

///
/// Function prototypes
///

extern int main_diverge(int argc, char *argv[]);
extern int main_diverge_bam(pop_diverge_parser *param);
extern int main_diverge_vcf(pop_diverge_parser *param);
extern int make_diverge(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data);

#endif

