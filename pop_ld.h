/** \file pop_ld.h
 *  \brief Header for the pop_ld.cpp file
 *  \author Daniel Garrigan
 *  \version 0.4
*/

#ifndef POP_LD_H
#define POP_LD_H

///
/// Definitions
///

//
// Define data structures
//

/*!
 * \class ldData
 * \brief A derived class for passing parameters and data to the ld function
 */
class ldData: public popbamData
{
public:
    // constructor
    ldData(const popbamOptions&);

    // destructor
    ~ldData(void);

    // member variables
    int output;                  //!< Analysis output option
    uint32_t *pop_cov;           //!< Boolean for population coverage
    int minSNPs;                 //!< Minimum number of snps for a window to be considered
    uint16_t minFreq;            //!< Minimum allele count in LD calculation
    int *num_snps;               //!< Number of SNPs in a given window
    int npops;
    double minSites;             //!< Minimum proportion of aligned sites for a window to be considered
    double *omegamax;            //!< Pointer to array of omega_max values
    double *wallb;               //!< Pointer to array of Wall's B statistic
    double *wallq;               //!< Pointer to array of Wall's Q statistic
    double *zns;                 //!< Pointer to array of ZnS values

    // member functions
    int calcZns(void);
    int calcOmegamax(void);
    int calcWall(void);
    int allocLD(void);
    int printLD(const std::string);
};

///
/// Function prototypes
///

extern int mainLD(int argc, char *argv[]);

#endif

