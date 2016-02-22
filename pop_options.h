/** \file pop_options.h
 *  \brief Header for pop_options.cpp
 *  \author Daniel Garrigan
 *  \version 0.4
 */
 
#ifndef POP_OPTIONS_H
#define POP_OPTIONS_H

///
/// Define classes
///

class popbamOptions
{
public:
	// constructor
	popbamOptions(int, char**);

	// destructor
	~popbamOptions() {}

	// member variables
	samfile_t *bam_in;          //!< BAM input file stream
	faidx_t *fai_file;          //!< Fasta reference file index
	bam_index_t *idx;           //!< Pointer to the BAM input file index
	bam_header_t *h;            //!< Pointer to the header text for the input BAM file
	uint16_t flag;              //!< Bit flag to hold user options
	int output;                 //!< Analysis output option
	int errorCount;             //!< Flag to indicate error in reading user options
	int minDepth;               //!< User-specified minimumm read depth
	int maxDepth;               //!< User-specified maximum read depth
	int minRMSQ;                //!< User-specified minimum rms mapping quality
	int minSNPQ;                //!< User-specified minimum SNP quality score
	uint32_t winSize;           //!< User-specified window size in kilobases
	uint8_t minMapQ;            //!< User-specified minimum individual read mapping quality
	uint8_t minBaseQ;           //!< User-specified minimum inidividual base quality
	double minSites;            //!< User-specified minimum number of aligned sites to perform analysis
	double minPop;              //!< Minimum proportion of samples present
	double hetPrior;            //!< Prior probability for calling heterozygous genotypes
	std::string dist;           //!< Pointer to the name of the desired distance metric	(-d switch)
	std::string bamfile;        //!< File name for the input BAM file
	std::string reffile;        //!< File name for the input reference Fasta file
	std::string headfile;       //!< File name for optional BAM header input file
	std::string region;         //!< Region on which to perform the analysis
	std::string errorMsg;       //!< String to hold any error messages
	std::string popFunc;        //!< The popbam function being invoked

	// member functions
	int checkBAM(void);
};

#endif
