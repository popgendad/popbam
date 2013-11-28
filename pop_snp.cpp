/** \file pop_snp.cpp
 *  \brief Functions for extracting SNP calls from BAM files
 *  \author Daniel Garrigan
 *  \version 0.4
*/
#include "pop_snp.h"

int mainSNP(int argc, char *argv[])
{
	bool found = false;           //! is the outgroup found?
	int chr = 0;                  //! chromosome identifier
	int beg = 0;                  //! beginning coordinate for analysis
	int end = 0;                  //! end coordinate for analysis
	int ref = 0;                  //! reference allele
	long nWindows = 0;            //! total number of windows
	std::string msg;              //! string for error message
	bam_sample_t *sm = nullptr;   //!< Pointer to the sample information for the input BAM file
	bam_plbuf_t *buf = nullptr;   //! pileup buffer

	// initialize user command line options
	popbamOptions p(argc, argv);

	// check input BAM file for errors
	p.checkBAM();

	// initialize the sample data structure
	sm = bam_smpl_init();

	// add samples
	bam_smpl_add(sm, &p);

	// initialize the snp data structre
	snpData t(p);
	t.sm = sm;

	// initialize error model
	t.em = errmod_init(0.17);

	// if outgroup option is used check to make sure it exists
	if (p.flag & BAM_OUTGROUP)
	{
		for (int i = 0; i < t.sm->n; ++i)
		{
			if (strcmp(t.sm->smpl[i], t.outgroup.c_str()) == 0)
			{
				t.outidx = i;
				found = true;
			}
		}

		if (!found)
		{
			msg = "Specified outgroup " + t.outgroup + " not found";
			fatalError(msg);
		}
	}

	// parse genomic region
	int k = bam_parse_region(p.h, p.region, &chr, &beg, &end);
	if (k < 0)
	{
		msg = "Bad genome coordinates: " + p.region;
		fatalError(msg);
	}

	// fetch reference sequence
	t.ref_base = faidx_fetch_seq(p.fai_file, p.h->target_name[chr], 0, 0x7fffffff, &(t.len));

	// calculate the number of windows
	if (p.flag & BAM_WINDOW)
		nWindows = ((end - beg) - 1) / p.winSize;
	else
	{
		p.winSize = end - beg;
		nWindows = 1;
	}

	// iterate through all windows along specified genomic region
	for (long j = 0; j < nWindows; ++j)
	{
		// construct genome coordinate string
		std::string scaffold_name(p.h->target_name[chr]);
		std::ostringstream winc(scaffold_name);

		winc.seekp(0, std::ios::end);
		winc << ':' << beg + (j * p.winSize) + 1 << '-' << ((j + 1) * p.winSize) + (beg - 1);

		std::string winCoord = winc.str();

		// initialize number of sites to zero
		t.num_sites = 0;

		// parse the BAM file and check if region is retrieved from the reference
		if (p.flag & BAM_WINDOW)
		{
			k = bam_parse_region(p.h, winCoord, &ref, &(t.beg), &(t.end));

			if (k < 0)
			{
				msg = "Bad window coordinates " + winCoord;
				fatalError(msg);
			}
		}
		else
		{
			ref = chr;
			t.beg = beg;
			t.end = end;

			if (ref < 0)
			{
				msg = "Bad scaffold name: " + p.region;
				fatalError(msg);
			}
		}

		// initialize diverge specific variables
		t.allocSNP();

		// create population assignments
		t.assignPops(&p);

		// print ms header if first window iteration
		if ((t.output == 2) && (j == 0))
			t.printMSHeader(nWindows);

		// initialize pileup
		buf = bam_plbuf_init(makeSNP, &t);

		// fetch region from bam file
		if ((bam_fetch(p.bam_in->x.bam, p.idx, ref, t.beg, t.end, buf, fetch_func)) < 0)
		{
			msg = "Failed to retrieve region " + p.region + " due to corrupted BAM index file";
			fatalError(msg);
		}

		// finalize pileup
		bam_plbuf_push(0, buf);

		// print results to stdout
		t.print_SNP(std::string(p.h->target_name[chr]));

		// take out the garbage
		bam_plbuf_destroy(buf);
	}
	// end of window interation

	errmod_destroy(t.em);
	samclose(p.bam_in);
	bam_index_destroy(p.idx);
	bam_smpl_destroy(sm);
	free(t.ref_base);

	return 0;
}

int makeSNP(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
{
	int i = 0;
	int fq = 0;
	unsigned long long sample_cov = 0;
	unsigned long long *cb = nullptr;
	snpData *t = nullptr;

	// get control data structure
	t = (snpData*)data;

	// only consider sites located in designated region
	if ((t->beg <= (int)pos) && (t->end > (int)pos))
	{
		// call bases
		cb = callBase(t, n, pl);

		// resolve heterozygous sites
		if (!(t->flag & BAM_HETEROZYGOTE))
			cleanHeterozygotes(t->sm->n, cb, (int)t->ref_base[pos], t->minSNPQ);

		// determine if site is segregating
		fq = segBase(t->sm->n, cb, t->ref_base[pos], t->minSNPQ);

		// determine how many samples pass the quality filters
		sample_cov = qualFilter(t->sm->n, cb, t->minRMSQ, t->minDepth, t->maxDepth);
		
		unsigned int *ncov = nullptr;
		ncov = new unsigned int [t->sm->npops]();

		// determine population coverage
		for (i = 0; i < t->sm->npops; ++i)
		{
			unsigned long long pc = 0;

			pc = sample_cov & t->pop_mask[i];
			ncov[i] = bitcount64(pc);

			unsigned int req = (unsigned int)((t->minPop * t->pop_nsmpl[i]) + 0.4999);

			if (ncov[i] >= req)
				t->pop_cov[t->num_sites] |= 0x1U << i;
		}	

		if (t->pop_cov[t->num_sites] > 0)
		{
			if (fq > 0)
			{
				// calculate the site type
				for (i = 0; i < t->sm->npops; ++i)
					t->ncov[i][t->segsites] = ncov[i];
				t->types[t->segsites] = calculateSiteType(t->sm->n, cb);

				// add to the haplotype matrix
				t->hap.pos[t->segsites] = pos;
				t->hap.ref[t->segsites] = bam_nt16_table[(int)t->ref_base[pos]];

				for (i = 0; i < t->sm->n; i++)
				{
					t->hap.rms[i][t->segsites] = (cb[i] >> (CHAR_BIT * 6)) & 0xffff;
					t->hap.snpq[i][t->segsites] = (cb[i] >> (CHAR_BIT * 4)) & 0xffff;
					t->hap.num_reads[i][t->segsites] = (cb[i]>>(CHAR_BIT * 2)) & 0xffff;
					t->hap.base[i][t->segsites] = bam_nt16_table[(int)iupac[(cb[i] >> CHAR_BIT) & 0xff]];

					if (cb[i] & 0x2ULL)
						t->hap.seq[i][t->segsites/64] |= 0x1ULL << t->segsites % 64;
				}
				t->hap.idx[t->segsites] = t->num_sites;
				t->segsites++;
			}
			t->num_sites++;
		}

		// take out the garbage
		delete [] cb;
		delete [] ncov;
	}

	return 0;
}

int snpData::print_SNP(const std::string scaffold)
{
	snp_func fp[3] = {&snpData::printSNP, &snpData::printSweep, &snpData::printMS};
	(this->*fp[output])(scaffold);
}

int snpData::printSNP(const std::string scaffold)
{
	int i = 0;
	int j = 0;

	for (i = 0; i < segsites; i++)
	{
		std::stringstream out;
		out << scaffold << '\t' << hap.pos[i] + 1 << '\t';
		out << bam_nt16_rev_table[hap.ref[i]];

		for (j = 0; j < sm->n; j++)
		{
			out << '\t' << bam_nt16_rev_table[hap.base[j][i]];
			out << '\t' << hap.snpq[j][i];
			out << '\t' << hap.rms[j][i];
			out << '\t' << hap.num_reads[j][i];
		}

		std::cout << out.str() << std::endl;
	}

	return 0;
}

int snpData::printSweep(const std::string scaffold)
{
	int i = 0;
	int j = 0;
	unsigned short freq = 0;
	unsigned short pop_n = 0;
	unsigned long long pop_type = 0;

	for (i = 0; i < segsites; i++)
	{
		std::stringstream out;
		out << scaffold << '\t' << hap.pos[i] + 1;

		for (j = 0; j < sm->npops; j++)
		{
			// population-specific site type
			pop_type = types[i] & pop_mask[j];

			// polarize the mutation at the site
			if ((flag & BAM_OUTGROUP) && CHECK_BIT(types[i], outidx))
				freq = ncov[j][i] - bitcount64(pop_type);
			else
				freq = bitcount64(pop_type);

			out << '\t' << freq << '\t' << ncov[j][i];
		}

		std::cout << out.str() << std::endl;
	}

	return 0;
}

int snpData::printMS(const std::string scaffold)
{
	int i = 0;
	int j = 0;
	std::stringstream out;

	out << "//\n" << "segsites: " << segsites << "\npositions: ";

	for (i = 0; i < segsites; i++)
		out << std::setprecision(8) << (double)(hap.pos[i] - beg) / (end - beg) << ' ';
	out << '\n';

	for (i = 0; i < sm->n; i++)
	{
		for (j = 0; j < segsites; j++)
		{
			if ((flag & BAM_OUTGROUP) && CHECK_BIT(types[hap.idx[j]], outidx))
			{
				if (CHECK_BIT(hap.seq[i][j/64], j % 64))
					out << '0';
				else
					out << '1';
			}
			else
			{
				if (CHECK_BIT(hap.seq[i][j/64], j % 64))
					out << '1';
				else
					out << '0';
			}
		}
		out << '\n';
	}
	std::cout << out.str() << std::endl;
	return 0;
}

int snpData::printMSHeader(long nwindows)
{
	int i = 0;
	std::stringstream out;

	if (sm->npops > 1)
	{
		out << "ms " << sm->n << ' ' << nwindows;
		out << " -t 5.0 -I " << sm->npops << ' ';

		for (i = 0; i < sm->npops; i++)
			out << (int)(pop_nsmpl[i]) << ' ';
	}
	else
	{
		out << "ms " << sm->n << ' ' << nwindows << " -t 5.0 ";
	}

	out << "\n1350154902";
	std::cout << out.str() << std::endl << std::endl;
}


snpData::snpData(const popbamOptions &p)
{
	// inherit values from popbamOptions
	flag = p.flag;
	minDepth = p.minDepth;
	maxDepth = p.maxDepth;
	minRMSQ = p.minRMSQ;
	minSNPQ = p.minSNPQ;
	minMapQ = p.minMapQ;
	minBaseQ = p.minBaseQ;
	hetPrior = p.hetPrior;
	minPop = p.minPop;
	output = p.output;

	// initialize native variables
	derived_type = SNP;
	outidx = 0;
}

int snpData::allocSNP(void)
{
	int i = 0;
	const int length = end - beg;
	const int n = sm->n;
	const int npops = sm->npops;

	segsites = 0;

	try
	{
		types = new unsigned long long [length]();
		pop_mask = new unsigned long long [npops]();
		pop_nsmpl = new unsigned char [npops]();
		pop_cov = new unsigned int [length]();
		hap.pos = new unsigned int [length]();
		hap.idx = new unsigned int [length]();
		hap.ref = new unsigned char [length]();
		hap.seq = new unsigned long long* [n];
		hap.base = new unsigned char* [n];
		hap.rms = new unsigned short* [n];
		hap.snpq = new unsigned short* [n];
		hap.num_reads = new unsigned short* [n];
		ncov = new unsigned int* [npops];

		for (i = 0; i < n; i++)
		{
			hap.seq[i] = new unsigned long long [length]();
			hap.base[i] = new unsigned char [length]();
			hap.rms[i] = new unsigned short [length]();
			hap.snpq[i] = new unsigned short [length]();
			hap.num_reads[i] = new unsigned short [length]();
		}

		for (i = 0; i < npops; i++)
			ncov[i] = new unsigned int [length]();
	}
	catch (std::bad_alloc& ba)
	{
		std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
	}
}

snpData::~snpData(void)
{
	int i = 0;
	int npops = sm->npops;

	delete [] pop_mask;
	delete [] pop_nsmpl;
	delete [] pop_cov;
	delete [] types;
	delete [] hap.pos;
	delete [] hap.idx;
	delete [] hap.ref;
	for (i = 0; i < sm->n; i++)
	{
		delete [] hap.seq[i];
		delete [] hap.base[i];
		delete [] hap.num_reads[i];
		delete [] hap.snpq[i];
		delete [] hap.rms[i];
	}
	for (i = 0; i < npops; i++)
		delete [] ncov[i];
	delete [] ncov;
	delete [] hap.seq;
	delete [] hap.base;
	delete [] hap.snpq;
	delete [] hap.rms;
	delete [] hap.num_reads;
}

void snpData::printUsage(const std::string msg)
{
	std::cerr << msg << std::endl << std::endl;
	std::cerr << "Usage:   popbam snp [options] <in.bam> [region]" << std::endl << std::endl;
	std::cerr << "Options: -i          base qualities are Illumina 1.3+               [ default: Sanger ]" << std::endl;
	std::cerr << "         -h  FILE    Input header file                              [ default: none ]" << std::endl;
	std::cerr << "         -v          output variant sites only                      [ default: all sites ]" << std::endl;
	std::cerr << "         -z  FLT     output heterozygous base calls                 [ default: consensus ]" << std::endl;
	std::cerr << "         -w  INT     use sliding window of size (kb)" << std::endl;
	std::cerr << "         -p  STR     sample name of outgroup                        [ default: reference ]" << std::endl;
	std::cerr << "         -o  INT     output format                                  [ default: 0 ]" << std::endl;
	std::cerr << "                     0 : popbam snp format" << std::endl;
	std::cerr << "                     1 : SweepFinder snp format" << std::endl;
	std::cerr << "                     2 : MS format" << std::endl;
	std::cerr << "         -n  FLT     minimum proportion of population covered       [ default: 1.0 ]" << std::endl;
	std::cerr << "         -f  FILE    Reference fastA file" << std::endl;
	std::cerr << "         -m  INT     minimum read coverage                          [ default: 3 ]" << std::endl;
	std::cerr << "         -x  INT     maximum read coverage                          [ default: 255 ]" << std::endl;
	std::cerr << "         -q  INT     minimum rms mapping quality                    [ default: 25 ]" << std::endl;
	std::cerr << "         -s  INT     minimum snp quality                            [ default: 25 ]" << std::endl;
	std::cerr << "         -a  INT     minimum map quality                            [ default: 13 ]" << std::endl;
	std::cerr << "         -b  INT     minimum base quality                           [ default: 13 ]" << std::endl << std::endl;
	exit(EXIT_FAILURE);
}
