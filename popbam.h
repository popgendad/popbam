/** \file popbam.h
 *  \brief Header for the popbam program
 *  \author Daniel Garrigan
 *  \version 0.4
 */
#ifndef POPBAM_H
#define POPBAM_H

///
/// Define classes
///

/*!
 * \class popbamData
 * \brief The abstract base class for passing parameters and data
 */
class popbamData
{
	public:
		// default constructor
		popbamData();

		// destructor
		~popbamData() {}

		// member functions
		int assignPops(const popbamOptions *p);

		// member variables
		std::string bamfile;           //!< Name of bamfile used for indexing purposes
		bam_sample_t *sm;              //!< Pointer to the sample information for the input BAM file
		char *ref_base;                //!< Reference sequence string for specified region
		int tid;                       //!< Reference chromosome/scaffold identifier
		int beg;                       //!< Reference coordinate of the beginning of the current region
		int end;                       //!< Reference coordinate of the end of current region
		int len;                       //!< Length of the reference sequence for current region
		uint16_t flag;                 //!< Bit flag to hold user options
		int num_sites;                 //!< Total number of aligned sites
		int segsites;                  //!< Total number of segregating sites in entire sample
		uint8_t *pop_nsmpl;            //!< Sample size per population
		uint64_t *types;               //!< The site type for each aligned site
		uint64_t *pop_mask;            //!< Bit mask for which individuals are in which population
		int minDepth;                  //!< User-specified minimumm read depth
		int maxDepth;                  //!< User-specified maximum read depth
		int minRMSQ;                   //!< User-specified minimum rms mapping quality
		int minSNPQ;                   //!< User-specified minimum SNP quality score
		uint8_t minMapQ;               //!< User-specified minimum individual read mapping quality
		uint8_t minBaseQ;              //!< User-specified minimum inidividual base quality
		double hetPrior;               //!< Prior probability of heterozygous genotype
		errmod_t *em;                  //!< Error model data structure
		popbam_func_t derived_type;    //!< Type of the derived class
};

#endif
