/** \file pop_sample.hpp
 *  \brief Header for the pop_sample.cpp file
 *  \author Daniel Garrigan
 *  \version 0.5
 */

#ifndef POP_SAMPLE_HPP
#define POP_SAMPLE_HPP

#include "bam.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"

/*!
 * \struct bam_sample_t
 * \brief A structure to represent a sample in the BAM file
 */

typedef struct __bam_sample_t
{
    int npops;          //!< Number of populations in the BAM file
    int b;              //!< Counter for population configuration
    int n;              //!< Number of samples in the BAM file
    int m;              //!< Counter for sample configuration
    char **smpl;        //!< Pointer to array of sample names
    char **popul;       //!< Pointer to array of population names
    void *rg2smid;      //!< Pointer to hash for read group to sample id lookup
    void *sm2popid;     //!< Pointer to hash for sample to population id lookup
    void *sm2id;        //!< Pointer to hash for sample to identifier lookup
    void *pop2sm;       //!< Pointer to hash for population to sample lookup
} bam_sample_t;

KHASH_MAP_INIT_STR(sm, int)

///
/// Function prototypes
///

/*!
 * \fn bam_sample_t *bam_smpl_init(void)
 * \brief Initialize the sample data structure
 */

extern bam_sample_t *bam_smpl_init (void);

/*!
 * \fn int bam_smpl_add(bam_sample_t *sm, const char *bamfile)
 * \brief Add a sample data structure
 * \param sm Pointer to sample data structure
 * \param abs Pointer to name of input BAM file
 * \param txt Pointer to unformatted BAM header txt
 */

extern int bam_smpl_add (bam_sample_t *sm, bam_header_t *h, const char *bamfile);

/*!
 * \fn int bam_smpl_rg2smid(const bam_sample_t *sm, const char *fn, const char *rg, kstring_t *str)
 * \brief Get the sample id of a read group
 * \param sm Pointer to sample data structure
 * \param fn Pointer to the name of the input BAM file
 * \param rg Pointer to the name of the read group
 * \param str Pointer to the name of the sample
 */

extern int bam_smpl_rg2smid (const bam_sample_t *sm, const char *fn,
                             const char *rg, kstring_t *str);

/*!
 * \fn int bam_smpl_sm2popid(const bam_sample_t *sm, const char *fn, const char *smpl, kstring_t *str)
 * \brief Get the population id of a sample
 * \param sm Pointer to sample data structure
 * \param fn Pointer to the name of the input BAM file
 * \param smpl Pointer to the name of the sample
 * \param str Pointer to the name of the population
 */

extern int bam_smpl_sm2popid (const bam_sample_t *sm, const char *fn,
                             const char *smpl, kstring_t *str);

/*!
 * \fn void bam_smpl_destroy(bam_sample_t *sm)
 * \brief Free a sample data structure from memory
 * \param sm Pointer to sample data structure
 */
extern void bam_smpl_destroy (bam_sample_t *sm);

/*!
 * \fn static void add_sample_pair(bam_sample_t *sm, khash_t(sm) *pop2sm, const char *key, const char *val)
 * \brief
 * \param sm Pointer to sample data structure
 * \param pop2sm
 * \param key
 * \param val
 */

extern void add_sample_pair (bam_sample_t *sm, khash_t(sm) *pop2sm,
                             const char *key, const char *val);

/*!
 * \fn static void add_pop_pair(bam_sample_t *sm, khash_t(sm) *pop2sm, const char *key, const char *val)
 * \brief
 * \param sm Pointer to sample data structure
 * \param pop2sm
 * \param key
 * \param val
 */

extern void add_pop_pair (bam_sample_t *sm, khash_t(sm) *pop2sm, const char *key,
                          const char *val);

#endif
