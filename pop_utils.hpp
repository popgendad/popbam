/** \file pop_utils.hpp
 *  \brief Header for the pop_utils.cpp file
 *  \author Daniel Garrigan
 *  \version 0.5
 */
#ifndef POP_UTILS_HPP
#define POP_UTILS_HPP

#include <cstdint>
#include <cmath>

#include "sam.h"

///
/// Definitions
///

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

/*! \def NBASES
 *  \brief The number of possible bases
 */
#define NBASES 4

/*! \def popbam_func_t
 *  \brief A enum data type that holds the popbam function identifier
 */
enum popbam_func_t {SNP, FASTA, DIVERGE, HAPLO, TREE, NUCDIV, LD, SFS};

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

///
/// Function prototypes
///

/*!
 * \fn uint64_t gl2cns(float q[16], uint16_t k)
 * \brief Calculates a consensus base call from genotype likelihoods
 * \param q  Probabilites associated with each base
 * \param k  Number of reads mapping to a position in an individual
 */

extern uint64_t gl2cns (float q[16], uint16_t k);

/*!
 * \fn uint64_t  qualFilter(int num_samples, uint64_t *cb, int min_rmsQ, int min_depth, int max_depth)
 * \brief Filters data based on quality threshholds
 * \param num_samples  The number of samples in the pileup
 * \param cb  The consensus base call information for the individual
 * \param min_rmsQ  Minimum root-mean square of mapping quality for site to be considered
 * \param min_depth  Minimum read depth per individual for site to be considered
 * \param max_depth  Maximum read depth per individual for site to be considered
 */

extern uint64_t qual_filter (int num_samples, uint64_t *cb, int min_rmsQ,
                             int min_depth, int max_depth);

/*!
 * \fn int seg_base(int num_samples, uint64_t *cb, char ref, int min_snpq)
 * \brief Determines whether a base position is segregating or not
 * \param num_samples  The number of samples in the pileup
 * \param cb  The consensus base call information for the individual
 * \param ref  The reference base
 * \param min_snpq  The minimum acceptable SNP score to consider a site a variant
 */

extern int seg_base (int num_samples, uint64_t *cb, char ref, int min_snpq);

/*!
 * \fn void clean_hets(int num_samples, uint64_t *cb, int ref, int min_snpq)
 * \brief Reconfigures heterozygous base calls
 * \param num_samples  The number of samples in the pileup
 * \param cb  The consensus base call information for the individual
 * \param ref  The reference base
 * \param min_snpq  The minimum acceptable SNP score to consider a site a variant
 */

extern void clean_hets (int num_samples, uint64_t *cb, int ref, int min_snpq);

/*!
 * \fn static double* logbinomial_table(const int n_size)
 * \brief Construct logbinomial table
 * \param n_size Size of the table
 */

extern double *logbinomial_table (const int n_size);

/*!
 * \fn char *get_refid(char *htext)
 * \brief Function to extract reference identifier from BAM header
 * \param htext Pointer to unformatted BAM header text
 */

extern char *get_refid (char *htext);

/*!
 * \fn int fetch_func (const bam1_t *b, void *data)
 * \brief Assigns functions to the pileup push
 * \param b Pointer to the alignment structure
 * \param data User defined data structure
 */
extern int fetch_func (const bam1_t *b, void *data);

/*!
 * \fn inline uint16_t bitcount64 (uint64_t x)
 * \brief Function count the number of bits set in a 64-bit integer
 * \param x the 64-bit integer
 * \return unsigned integer
 * Returns the number of bits set
 */

inline uint16_t
bitcount64 (uint64_t x)
{
    x = (x & 0x5555555555555555ULL) + ((x >> 1) & 0x5555555555555555ULL);
    x = (x & 0x3333333333333333ULL) + ((x >> 2) & 0x3333333333333333ULL);
    x = (x & 0x0F0F0F0F0F0F0F0FULL) + ((x >> 4) & 0x0F0F0F0F0F0F0F0FULL);
    return (x * 0x0101010101010101ULL) >> 56;
}

/*!
 * \fn inline uint32_t hamming_distance (uint64_t x, uint64_t y)
 * \brief Function to compute the hamming distance between two 64-bit integers
 * \param x the first 64-bit integer
 * \param y the second 64-bit integer
 * \return unsigned int of hamming distance
 * Returns the number of bits set
 */

inline uint32_t
hamming_distance (uint64_t x, uint64_t y)
{
    uint32_t dist = 0;
    uint64_t val = x^y;

    while (val)
        {
            ++dist;
            val &= val - 1;
        }
    return dist;
}

/*!
 * \fn inline uint64_t calculate_sitetype (int n, uint64_t *cb)
 * \brief Function to calculate site type representation
 * \param n
 * \param cb Consensus base data structure
 * \return encoded site type
 */

inline uint64_t
calculate_sitetype (int n, uint64_t *cb)
{
    uint64_t site_type = 0;

    for (int i = 0; i < n; i++)
        {
            if ((cb[i] & 0x3ULL) == 0x3ULL)
                {
                    site_type |= 0x1ULL << i;
                }
        }
    return site_type;
}

#endif
