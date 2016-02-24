/** \file pop_global.h
 *  \brief Global declarations for the popbam program
 *  \author Daniel Garrigan
 *  \version 0.4
 */

#ifndef POP_GLOBAL_H
#define POP_GLOBAL_H

#include "tables.h"
#include "ksort.h"
#include "khash.h"
#include "kstring.h"
#include "bam.h"
#include "sam.h"
#include "faidx.h"

///
/// Definitions
///

/*! \def BAM_VARIANT
 *  \brief Flag for the -v command line switch-- only output variable sites
 */
#define BAM_VARIANT 0x01

/*! \def BAM_ILLUMINA
 *  \brief Flag for the -i command line switch-- qualities are in Illumina 1.3+ format
 */
#define BAM_ILLUMINA 0x02

/*! \def BAM_WINDOW
 *  \brief Flag for the -w command line switch-- do a sliding window analysis
 */
#define BAM_WINDOW 0x04

/*! \def BAM_MINPOPSAMPLE
 *  \brief Flag for user designated minimum sample sizes per population
 */
// TODO: This is deprecated
#define BAM_MINPOPSAMPLE 0x08

/*! \def BAM_SUBSTITUTE
 *  \brief Flag for only counting fixed substitutions in diverge function
 */
#define BAM_SUBSTITUTE 0x10

/*! \def BAM_HETEROZYGOTE
 *  \brief Flag for outputting heterozygous site calls in snp function
 */
#define BAM_HETEROZYGOTE 0x20

/*! \def BAM_OUTGROUP
 *  \brief Flag for changing the outgroup from the reference
 */
#define BAM_OUTGROUP 0x40

/*! \def BAM_HEADERIN
 *  \brief Flag for presence of user BAM header file
 */
#define BAM_HEADERIN 0x80

/*! \def BAM_NOSINGLETONS
 *  \brief Flag to exclude singleton polymorphisms from the analysis
 */
#define BAM_NOSINGLETONS 0x100

/*! \def POPBAM_RELEASE
 *  \brief Version number of popbam program
 */
#define POPBAM_RELEASE "0.4b"

/*! \def NBASES
 *  \brief The number of possible bases
 */
#define NBASES 4

/*! \def IUPAC_N
 *  \brief The integer representation of the IUPAC ambiguity symbol 'N'
 */
#define IUPAC_N 0xf

/*! \def KB
 *  \brief Integer for length of a kilobase
 */
#define KB 1000

/*! \def CHECK_BIT(var,pos)
 *  \brief A macro to check if a bit is set at pos in the unsigned long long var
 */
#define CHECK_BIT(var,pos) ((var) & (0x1ULL << (pos)))

/*! \def SEG_IDX(segsite)
 *  \brief A macro access index of a segregating site
 */
#define SEG_IDX(seg) (((seg) - 1) / 64)

/*! \def UTIDX(nrows,row,col)
 *  \brief A macro to access index of 1D array from 2D data structure
 */
#define UTIDX(n,i,j) ((2*(n)-((i)+1))*(((i)+1)-1)/2-((i)+1)+((j)+1)-1)

/*! \def SQ(x)
 *  \brief A macro to calculate the square of its argument
 */
#define SQ(x) ((x) * (x))

/*! \def BINOM(x)
 *  \brief A macro to calculate binomial coefficient
 */
#define BINOM(x) ((x) * ((x) - 1) / 2)

//
// Define data structures
//

/*!
 * struct hData_t
 * \brief A structure to represent a haplotype data set
 */
typedef struct
{
    uint64_t **seq;             //!< binary encoding of haplotype data
    uint32_t *pos;              //!< reference coordinate for each position
    uint32_t *idx;              //!< position index of each segregating site
    uint8_t *ref;               //!< reference allele at each position
    uint8_t **base;             //!< consensus base at each position in each individual
    uint16_t **rms;             //!< root mean square mapping score at each position
    uint16_t **snpq;            //!< SNP quality score at each position
    uint16_t **num_reads;       //!< number of reads at each position in each individual
} hData_t;

//
// Define some global variables
//

/*! \def popbam_func_t
 *  \brief A enum data type that holds the popbam function identifier
 */
enum popbam_func_t {SNP, FASTA, DIVERGE, HAPLO, TREE, NUCDIV, LD, SFS};

#endif
