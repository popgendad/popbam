/** \file pop_snp.h
 *  \brief Header for the pop_snp.cpp file
 *  \author Daniel Garrigan
 *  \version 0.4
*/

#ifndef POP_SNP_H
#define POP_SNP_H

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
        uint32_t *pop_cov;                      //!< Boolean for population coverage
        uint32_t **ncov;                        //!< Sample size per population per segregating site
        int npops;
        uint64_t **pop_sample_mask;             //!< Bit mask for samples covered from a specific population
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

extern int mainSNP(int argc, char *argv[]);

#endif
