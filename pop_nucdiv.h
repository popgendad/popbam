/** \file pop_nucdiv.h
 *  \brief Header for the pop_nucdiv.cpp file
 *  \author Daniel Garrigan
 *  \version 0.5
*/

#ifndef POP_NUCDIV_H
#define POP_NUCDIV_H

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
    uint32_t *pop_cov;         //!< Boolean for population coverage
    uint32_t **ncov;           //!< Sample size per population per segregating site
    uint64_t *ns_within;       //!< Number of aligned sites with each population
    uint64_t *ns_between;      //!< Number of aligned sites between each pair of populations
    int npops;                 //!< Number of populations represented in BAM
    int *num_snps;             //!< Number of SNPs in a given window
    double minPop;             //!< Minimum proportion of samples present

    // member public functions
    int calcNucdiv(void);
    int allocNucdiv(void);
    int printNucdiv(const std::string);

private:
    // member private variables
    double minSites;           //!< User-specified minimum proportion of aligned sites to perform analysis
    double *piw;               //!< Array of within-population nucleotide diversity
    double *pib;               //!< Array of between-population Dxy values
};

///
/// Function prototypes
///

extern int mainNucdiv(int argc, char *argv[]);

#endif
