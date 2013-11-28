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

	std::string userFunc(argv[1]);

	if (userFunc.compare(std::string("snp")) == 0)
		return mainSNP(argc, argv);
	else if (userFunc.compare(std::string("haplo")) == 0)
		return mainHaplo(argc, argv);
	else if (userFunc.compare(std::string("diverge")) == 0)
		return mainDiverge(argc, argv);
	else if (userFunc.compare(std::string("tree")) == 0)
		return mainTree(argc, argv);
	else if (userFunc.compare(std::string("nucdiv")) == 0)
		return mainNucdiv(argc, argv);
	else if (userFunc.compare(std::string("ld")) == 0)
		return mainLD(argc, argv);
	else if (userFunc.compare(std::string("sfs")) == 0)
		return mainSFS(argc, argv);
	else if (userFunc.compare(std::string("fasta")) == 0)
		return 0;
	else
	{
		std::cerr << "Error: unrecognized command: " << userFunc << std::endl;
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
	minRMSQ = 25;
	minSNPQ = 25;
	minDepth = 3;
	maxDepth = 255;
	minMapQ = 13;
	minBaseQ = 13;
	hetPrior = 0.0001;
}

int popbamData::assignPops(const popbamOptions *p)
{
	int si = -1;
	kstring_t buf;

	memset(&buf, 0, sizeof(kstring_t));

	for (int i = 0; i < sm->n; i++)
	{
		if (sm->smpl[i])
			si = bam_smpl_sm2popid(sm, p->bamfile.c_str(), sm->smpl[i], &buf);

		if (si < 0)
			si = bam_smpl_sm2popid(sm, p->bamfile.c_str(), 0, &buf);

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

	return 0;
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
