/** \file pop_haplo.h
 *  \brief Header for the pop_haplo.cpp file
 *  \author Daniel Garrigan
 *  \version 0.5
*/

#ifndef POP_HAPLO_H
#define POP_HAPLO_H

//
// Define data structures
//

/*!
 * \class haploData
 * \brief A derived class for passing parameters and data to the haplo function
 */

class haploData: public popbamData
{
public:
    // constructor
    haploData(const popbamOptions&);

    // destructor
    ~haploData(void);

    // member public variables
    double minSites;                            //!< User-specified minimum number of aligned sites to perform analysis
    double minPop;                              //!< Minimum proportion of samples present
    int npops;                                  //!< Number of populations represented in BAM
    uint32_t *nsite_matrix;                     //!< Matrix of pairwise aligned sites
    uint32_t *diff_matrix;                      //!< Matrix of pairwise sequence difference
    uint32_t *minDxy;                           //!< Array of minimum between-population Dxy
    std::vector<std::vector<std::string>> hap;  //!< Vector of strings to hold the haplotypes

    // member public functions
    int allocHaplo(void);
    int calcHaplo(void);
    int printHaplo(const std::string);

private:
    // member private variables
    int output;                             //!< Analysis output option
    int *nhaps;                             //!< Array of number of haplotypes within populations
    double *hdiv;                           //!< Array of haplotype diversities within populations
    double *piw;                            //!< Array of within-population heterozygosities
    double *pib;                            //!< Array of between-population heterozygosities
    double *ehhs;                           //!< Array of site-specfic haplotype homozygosity statistics within populations

    // member private function
    int calcNhaps(void);
    int calcEHHS(void);
    int calcGmin(void);
};

///
/// Function prototypes
///

extern int mainHaplo(int argc, char *argv[]);

#endif
