/** \file pop_diverge.h
 *  \brief Header for the pop_diverge.cpp file
 *  \author Daniel Garrigan
 *  \version 0.5
*/

#ifndef POP_DIVERGE_H
#define POP_DIVERGE_H

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
    int output;                     //!< Analysis output option
    std::string outgroup;           //!< Sample name of outgroup to use
    int outidx;                     //!< Index of outgroup sequence
    hData_t hap;                    //!< Structure to hold haplotype data
    uint64_t *pop_sample_mask;      //!< Bit mask for samples covered from a specific population
    std::string dist;               //!< Pointer to the name of the desired distance metric	(-d switch)
    uint16_t *min_pop_n;            //!< Minimum sample size per population
    int *num_snps;                  //!< Number of SNPs in a given window
    int npops;                      //!< Number of populations represented in BAM

    // member public functions
    int calcDiverge(void);
    int allocDiverge(void);
    int setMinPop_n(void);
    int printDiverge(const std::string);

private:
    // member private variables
    int minSites;              //!< User-specified minimum number of aligned sites to perform analysis
    uint16_t *pop_div;         //!< Array of mean population divergence calculations
    uint16_t *ind_div;         //!< Array of individual divergence calculations
};

///
/// Function prototypes
///

extern int mainDiverge(int argc, char *argv[]);

#endif
