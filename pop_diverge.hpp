/** \file pop_diverge.hpp
 *  \brief Header for the pop_diverge.cpp file
 *  \author Daniel Garrigan
 *  \version 0.5
*/

#ifndef POP_DIVERGE_HPP
#define POP_DIVERGE_HPP

#include <cstdint>

#include "sam.h"
#include "bam.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/faidx.h"

#include "pop_error.hpp"
#include "pop_sample.hpp"
#include "pop_utils.hpp"
#include "pop_diverge_parser.hpp"

//
// Define data structures
//

typedef struct _diverge_data_bam
{
    int tid;                       //!< Reference chromosome/scaffold identifier
    int beg;                       //!< Reference coordinate of the beginning of the current region
    int end;                       //!< Reference coordinate of the end of current region
    int len;                       //!< Length of the reference sequence for current region
    int num_sites;                 //!< Total number of aligned sites
    int segsites;                  //!< Total number of segregating sites in entire sample
    int npops;                     //!< Number of populations represented in BAM
    int outidx;                    //!< Index of outgroup sequence
    int *num_snps;                 //!< Number of SNPs in a given window
    uint16_t *min_pop_n;           //!< Minimum sample size per population
    uint64_t *pop_sample_mask;     //!< Bit mask for samples covered from a specific population
    uint8_t *pop_nsmpl;            //!< Sample size per population
    uint64_t *types;               //!< The site type for each aligned site
    uint64_t *pop_mask;            //!< Bit mask for which individuals are in which population
    char *ref_base;                //!< Reference sequence string for specified region
    hData_t hap;                   //!< Structure to hold haplotype data
    int minSites;                  //!< User-specified minimum number of aligned sites to perform analysis
    uint16_t *pop_div;             //!< Array of mean population divergence calculations
    uint16_t *ind_div;             //!< Array of individual divergence calculations
    samfile_t *bam_in;             //! handle for input BAM file
    bam_sample_t *sm = NULL;       //! pointer to the sample information for the input BAM file
    bam_header_t *h;               //! pointer to BAM header
    bam_index_t *idx;              //! pointer to BAM index
    faidx_t *fai_file;             //! handle for indexed reference fastA file
    errmod_t *em;                  //! pointer to error model structure
} diverge_data_bam;

typedef struct _diverge_data_vcf
{
    int tid;                       //!< Reference chromosome/scaffold identifier
    int beg;                       //!< Reference coordinate of the beginning of the current region
    int end;                       //!< Reference coordinate of the end of current region
    int len;                       //!< Length of the reference sequence for current region
    int num_sites;                 //!< Total number of aligned sites
    int segsites;                  //!< Total number of segregating sites in entire sample
    int npops;                     //!< Number of populations represented in BAM
    int outidx;                    //!< Index of outgroup sequence
    int *num_snps;                 //!< Number of SNPs in a given window
    uint16_t *min_pop_n;           //!< Minimum sample size per population
    uint64_t *pop_sample_mask;     //!< Bit mask for samples covered from a specific population
    uint8_t *pop_nsmpl;            //!< Sample size per population
    uint64_t *types;               //!< The site type for each aligned site
    uint64_t *pop_mask;            //!< Bit mask for which individuals are in which population
    char *ref_base;                //!< Reference sequence string for specified region
    hData_t hap;                   //!< Structure to hold haplotype data
    int min_sites;                 //!< User-specified minimum number of aligned sites to perform analysis
    uint16_t *pop_div;             //!< Array of mean population divergence calculations
    uint16_t *ind_div;             //!< Array of individual divergence calculations
} diverge_data_vcf;

///
/// Function prototypes
///

extern int main_diverge (int argc, char *argv[]);
extern int main_diverge_bam (pop_diverge_parser *param);
extern int main_diverge_vcf (pop_diverge_parser *param);
extern int init_diverge_bam (diverge_data_bam *ddb, const pop_diverge_parser *param);
extern int init_diverge_vcf (diverge_data_vcf *ddv, const pop_diverge_parser *param);
extern int make_diverge (uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *par, void *data);
extern int calc_diverge (diverge_data_bam *ddb, const pop_diverge_parser *param);
extern int alloc_diverge_bam (diverge_data_bam *ddb, const pop_diverge_parser *param);
extern void dealloc_diverge_bam (diverge_data_bam *ddb, const pop_diverge_parser *param);
extern int alloc_diverge_vcf (diverge_data_vcf *ddv, const pop_diverge_parser *param);
extern void dealloc_diverge_vcf (diverge_data_vcf *ddv, const pop_diverge_parser *param);
extern int set_min_pop_n (diverge_data_bam *ddb, const pop_diverge_parser *param);
extern int print_diverge_bam (diverge_data_bam *ddb, const pop_diverge_parser *param, const char *scaffold, int win_size);
extern int print_diverge_vcf (diverge_data_vcf *ddv, const pop_diverge_parser *param, const char *scaffold, int win_size);
extern uint64_t *call_base_diverge (diverge_data_bam *ddb, const pop_diverge_parser *param, int n, const bam_pileup1_t *pl);
extern int check_BAM (diverge_data_bam *ddb, const pop_diverge_parser *param);
extern int assign_pops (diverge_data_bam *ddb, const pop_diverge_parser *param);

#endif

