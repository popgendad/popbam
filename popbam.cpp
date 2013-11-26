/** \file popbam.cpp
 *  \brief Main entry point for evolutionary analysis of BAM files
 *  \author Daniel Garrigan
 *  \version 0.4
*/
#include "popbam.h"
#include "tables.h"

const int bam_nt16_nt4_table[16] = { 4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4 };

const char iupac[16] = {'A','M','R','W','N','C','S','Y','N','N','G','K','N','N','N','T'};

const unsigned char bam_nt16_table[256] =
{
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	 1, 2, 4, 8, 15,15,15,15, 15,15,15,15, 15, 0,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15, 1,14, 2, 13,15,15, 4, 11,15,15,12, 15, 3,15,15,
	15,15, 5, 6,  8,15, 7, 9, 15,10,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15,
	15,15,15,15, 15,15,15,15, 15,15,15,15, 15,15,15,15
};

const unsigned char iupac_rev[256] =
{
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14, 0,14, 1, 14,14,14, 2, 14,14,14,14, 14,14,14,14,
	14,14,14,14,  3,14,14,14, 14,14,14,14, 14,14,14,14,
	14, 0,14, 1, 14,14,14, 2, 14,14,14,14, 14,14,14,14,
	14,14,14,14,  3,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14,
	14,14,14,14, 14,14,14,14, 14,14,14,14, 14,14,14,14
};

int main(int argc, char *argv[])
{
	if (argc < 2)
		return popbam_usage();
	if (!strcmp(argv[1], "snp"))
		return mainSNP(argc-1, argv+1);
	else if (!strcmp(argv[1], "haplo"))
		return mainHaplo(argc-1, argv+1);
	else if (!strcmp(argv[1], "diverge"))
		return mainDiverge(argc-1, argv+1);
	else if (!strcmp(argv[1], "tree"))
		return mainTree(argc-1, argv+1);
	else if (!strcmp(argv[1], "nucdiv"))
		return mainNucdiv(argc-1, argv+1);
	else if (!strcmp(argv[1], "ld"))
		return mainLD(argc-1, argv+1);
	else if (!strcmp(argv[1], "sfs"))
		return mainSFS(argc - 1, argv + 1);
	else if (!strcmp(argv[1], "fasta"))
		return 0;
	else
	{
		std::cerr << "Error: unrecognized command: " << argv[1] << std::endl;
		return 1;
	}
	return 0;
}

popbamData::popbamData(void)
{
	flag = 0x0;
	num_sites = 0;
	tid = -1;
	beg = 0;
	end = 0x7fffffff;
	min_rmsQ = 25;
	min_snpQ = 25;
	min_depth = 3;
	max_depth = 255;
	min_mapQ = 13;
	min_baseQ = 13;
	het_prior = 0.0001;
}

void popbamData::checkBAM(void)
{
	std::string msg;

	bam_in = samopen(bamfile.c_str(), "rb", 0);

	// check if BAM file is readable
	if (!bam_in)
	{
		msg = "Cannot read BAM file " + bamfile;
		fatalError(msg);
	}

	// check if BAM header is returned
	if (!bam_in->header)
	{
		msg = "Cannot read BAM header from file " + bamfile;
		fatalError(msg);
	}
	else
		h = bam_in->header;

	//read in new header text
	if (flag & BAM_HEADERIN)
	{
		std::ifstream headin(headfile);
		headin.seekg(0, std::ios::end);
		h->l_text = headin.tellg();
		headin.seekg(0, std::ios::beg);
		h->text = (char*)realloc(h->text, (size_t)h->l_text);
		headin.read(h->text, h->l_text);
		headin.close();
	}

	// check for bam index file
	if (!(idx = bam_index_load(bamfile.c_str())))
	{
		msg = "Index file not available for BAM file " + bamfile;
		fatalError(msg);
	}

	// check if fastA reference index is available
	fai_file = fai_load(reffile.c_str());
	if (!fai_file)
	{
		msg = "Failed to load index for fastA reference file: " + reffile;
		fatalError(msg);
	}
}

void popbamData::assign_pops(void)
{
	int si = -1;
	kstring_t buf;

	memset(&buf, 0, sizeof(kstring_t));

	for (int i = 0; i < sm->n; i++)
	{
		if (sm->smpl[i])
			si = bam_smpl_sm2popid(sm, bamfile.c_str(), sm->smpl[i], &buf);

		if (si < 0)
			si = bam_smpl_sm2popid(sm, bamfile.c_str(), 0, &buf);

		if (si < 0)
		{
			std::string msg;
			std::string missing_sample(sm->smpl[i]);
			msg = "Sample " + missing_sample + " not assigned to a population.\nPlease check BAM header file definitions";
			fatalError (msg);
		}

		pop_mask[si] |= 0x1ULL << i;
		pop_nsmpl[si]++;
	}
}

int popbam_usage(void)
{
	std::cerr << std::endl;
	std::cerr << "Program: popbam " << std::endl;
	std::cerr << "(Tools to perform evolutionary analysis from BAM files)" << std::endl;
	std::cerr << "Version: " << POPBAM_RELEASE << std::endl;
	std::cerr << "Usage: popbam <command> [options] <in.bam> [region]"  << std::endl << std::endl;
	std::cerr << "Commands:  snp       output consensus base calls" << std::endl;
	std::cerr << "           fasta     output alignment as multi-fasta file" << std::endl;
	std::cerr << "           haplo     output haplotype-based analyses" << std::endl;
	std::cerr << "           diverge   output divergence from reference" << std::endl;
	std::cerr << "           tree      output neighbor-joining trees" << std::endl;
	std::cerr << "           nucdiv    output nucleotide diversity statistics" << std::endl;
	std::cerr << "           ld        output linkage disequilibrium analysis" << std::endl;
	std::cerr << "           sfs       output site frequency spectrum analysis" << std::endl;
	std::cerr << std::endl;
	return 1;
}
