/** \file pop_sfs.h
 *  \brief Header for the pop_sfs.cpp file
 *  \author Daniel Garrigan
 *  \version 0.4
*/

#ifndef POP_SFS_H
#define POP_SFS_H

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
    double minSites;          //!< User-specified minimum proportion of aligned sites to perform analysis
    uint32_t *ns;             //!< Number of aligned sites within each population
    uint32_t **ncov;          //!< Sample size per population per segregating site
    int *num_snps;            //!< Number of SNPs in a given window
    uint32_t *pop_cov;        //!< Boolean for population coverage
    double minPop;            //!< Minimum proportion of samples present
    std::string outgroup;     //!< Sample name of outgroup to use
    int npops;                //!< Number of populations represented in BAM
    int outidx;               //!< Index of outgroup sequence
    int outpop;               //!< Population of outgroup sequence
    double **dw;              //!< Matrix of weights for Tajima's D calculation
    double **hw;              //!< Matrix of weights for Fay and Wu's H calculation
    double *a1;               //!< Constants for Tajima's D calculation
    double *a2;               //!< Constants for Tajima's D calculation
    double *e1;               //!< Constants for Tajima's D calculation
    double *e2;               //!< Constants for Tajima's D calculation
    double *td;               //!< Pointer to the array of Tajima's D statistics
    double *fwh;              //!< Pointer to the array of standardized Fay and Wu's H statistics

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

extern int mainSFS(int argc, char *argv[]);

#endif
