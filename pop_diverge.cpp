/** \file pop_diverge.cpp
 *  \brief Functions for calculating divergence from reference genome
 *  \author Daniel Garrigan
 *  \version 0.4
*/

#include "pop_diverge.h"

int mainDiverge(int argc, char *argv[])
{
	bool found = false;           //! is the outgroup sequence found?
	int chr = 0;                  //! chromosome identifier
	int beg = 0;                  //! beginning coordinate for analysis
	int end = 0;                  //! end coordinate for analysis
	int ref = 0;                  //! ref
	long nWindows = 0;            //! number of windows
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

	// initialize the diverge data structre
	divergeData t(p);
	t.sm = sm;

	// initialize error model
	t.em = errmod_init(0.17);

	// if outgroup option is used check to make sure it exists
	if (p.flag & BAM_OUTGROUP)
	{
		for (int i = 0; i < t.sm->n; i++)
			if (strcmp(t.sm->smpl[i], t.outgroup.c_str()) == 0)
			{
				t.outidx = i;
				found = true;
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
		t.allocDiverge();

		// create population assignments
		t.assignPops(&p);

		// set default minimum sample size as
		// the number of samples in the population
		t.setMinPop_n();

		// initialize pileup
		buf = bam_plbuf_init(makeDiverge, &t);

		// fetch region from bam file
		if ((bam_fetch(p.bam_in->x.bam, p.idx, ref, t.beg, t.end, buf, fetch_func)) < 0)
		{
			msg = "Failed to retrieve region " + p.region + " due to corrupted BAM index file";
			fatalError(msg);
		}

		// finalize pileup
		bam_plbuf_push(0, buf);

		// print results to stdout
		t.calcDiverge();
		t.printDiverge(std::string(p.h->target_name[chr]));

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

int makeDiverge(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
{
	int i = 0;
	int fq = 0;
	unsigned long long sample_cov = 0;
	unsigned long long *cb = nullptr;
	divergeData *t = nullptr;

	// get control data structure
	t = (divergeData*)data;

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

		for (i = 0; i < t->sm->npops; i++)
			t->pop_sample_mask[i] = sample_cov & t->pop_mask[i];

		if (bitcount64(sample_cov) == t->sm->n)
		{
			// calculate the site type
			t->types[t->num_sites] = calculateSiteType(t->sm->n, cb);

			if (fq > 0)
			{
				t->hap.pos[t->segsites] = pos;
				t->hap.ref[t->segsites] = (unsigned char)bam_nt16_table[(int)t->ref_base[pos]];
				for (i = 0; i < t->sm->n; i++)
				{
					t->hap.rms[i][t->segsites] = (cb[i] >> (CHAR_BIT * 6)) & 0xffff;
					t->hap.snpq[i][t->segsites] = (cb[i] >> (CHAR_BIT * 4)) & 0xffff;
					t->hap.num_reads[i][t->segsites] = (cb[i] >> (CHAR_BIT * 2)) & 0xffff;
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
	}
	return 0;
}

int divergeData::calcDiverge(void)
{
	int i = 0;
	int j = 0;
	unsigned short freq = 0;
	unsigned long long pop_type = 0;

	// calculate number of differences with reference sequence
	switch (output)
	{
		case 0:
			for (i = 0; i < sm->n; i++)
				for (j = 0; j <= SEG_IDX(segsites); j++)
					ind_div[i] += bitcount64(hap.seq[i][j]);
			break;
		case 1:
			for (i = 0; i < sm->npops; i++)
			{
				num_snps[i] = 0;
				for (j = 0; j < segsites; j++)
				{
					pop_type = types[hap.idx[j]] & pop_mask[i];

					// check if outgroup is different from reference
					if ((flag & BAM_OUTGROUP) && CHECK_BIT(types[hap.idx[j]], outidx))
						freq = pop_nsmpl[i] - bitcount64(pop_type);
					else
						freq = bitcount64(pop_type);
					if ((freq > 0) && (freq < pop_nsmpl[i]) && !(flag & BAM_NOSINGLETONS))
						++num_snps[i];
					else if ((freq > 1) && (freq < pop_nsmpl[i]) && (flag & BAM_NOSINGLETONS))
						++num_snps[i];
					else if (freq == pop_nsmpl[i])
						++pop_div[i];
				}
			}
			break;
		default:
			break;
	}

	return 0;
}

divergeData::divergeData(const popbamOptions &p)
{
	// inherit values from popbamOptions
	bamfile = p.bamfile;
	flag = p.flag;
	minDepth = p.minDepth;
	maxDepth = p.maxDepth;
	minRMSQ = p.minRMSQ;
	minSNPQ = p.minSNPQ;
	minMapQ = p.minMapQ;
	minBaseQ = p.minBaseQ;
	hetPrior = p.hetPrior;
	dist = p.dist;
	minSites = p.minSites;
	output = p.output;

	// initialize native variables
	derived_type = DIVERGE;
}

int divergeData::allocDiverge(void)
{
	int i = 0;
	int length = end - beg;
	int n = sm->n;
	int npops = sm->npops;

	segsites = 0;

	try
	{
		types = new unsigned long long [length]();
		pop_mask = new unsigned long long [npops]();
		pop_nsmpl = new unsigned char [npops]();
		pop_sample_mask = new unsigned long long [npops]();
		min_pop_n = new unsigned short [npops]();
		num_snps = new int [npops]();
		hap.pos = new unsigned int [length]();
		hap.idx = new unsigned int [length]();
		hap.ref = new unsigned char [length]();
		hap.seq = new unsigned long long* [n];
		hap.base = new unsigned char* [n];
		hap.rms = new unsigned short* [n];
		hap.snpq = new unsigned short* [n];
		hap.num_reads = new unsigned short* [n];
		switch (output)
		{
			case 0:
				ind_div = new unsigned short [n]();
				break;
			case 1:
				pop_div = new unsigned short [npops]();
			default:
				break;
		}
		for (i = 0; i < n; i++)
		{
			hap.seq[i] = new unsigned long long [length]();
			hap.base[i] = new unsigned char [length]();
			hap.rms[i] = new unsigned short [length]();
			hap.snpq[i] = new unsigned short [length]();
			hap.num_reads[i] = new unsigned short [length]();
		}
	}
	catch (std::bad_alloc& ba)
	{
		std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
	}

	return 0;
}

divergeData::~divergeData(void)
{
	int i = 0;

	delete [] pop_mask;
	delete [] types;
	delete [] pop_nsmpl;
	delete [] pop_sample_mask;
	delete [] num_snps;
	delete [] min_pop_n;
	delete [] hap.pos;
	delete [] hap.idx;
	delete [] hap.ref;
	switch (output)
	{
		case 0:
			delete [] ind_div;
			break;
		case 1:
			delete [] pop_div;
			break;
		default:
			break;
	}
	for (i = 0; i < sm->n; i++)
	{
		delete [] hap.seq[i];
		delete [] hap.base[i];
		delete [] hap.num_reads[i];
		delete [] hap.snpq[i];
		delete [] hap.rms[i];
	}
	delete [] hap.seq;
	delete [] hap.base;
	delete [] hap.snpq;
	delete [] hap.rms;
	delete [] hap.num_reads;
}

int divergeData::printDiverge(const std::string scaffold)
{
	int i = 0;
	double pdist = 0.0;
	double jc = 0.0;
	std::stringstream out;

	out << scaffold << '\t' << beg + 1 << '\t' << end + 1 << '\t' << num_sites;

	switch (output)
	{
		case 0:
			for (i = 0; i < sm->n; i++)
			{
				if (num_sites >= minSites)
				{
					if (dist == "pdist")
					{
						out << "\td[" << sm->smpl[i] << "]:";
						out << '\t' << std::fixed << std::setprecision(5) << (double)(ind_div[i]) / num_sites;
					}
					else if (dist == "jc")
					{
						pdist = (double)(ind_div[i]) / num_sites;
						jc = -0.75 * log(1.0 - pdist * (4.0 / 3.0));
						out << "\td[" << sm->smpl[i] << "]:";
						out << '\t' << std::fixed << std::setprecision(5) << jc;
					}
					else
						out << "\td[" << sm->smpl[i] << "]:\t" << std::setw(7) << "NA";
				}
				else
					out << "\td[" << sm->smpl[i] << "]:\t" << std::setw(7) << "NA";
			}
			break;
		case 1:
			for (i = 0; i < sm->npops; i++)
			{
				if (num_sites >= minSites)
				{
					out << "\tFixed[" << sm->popul[i] << "]:\t" << pop_div[i];
					out << "\tSeg[" << sm->popul[i] << "]:\t" << num_snps[i];
					out << "\td[" << sm->popul[i] << "]:";
					if (dist == "pdist")
					{
						if (flag & BAM_SUBSTITUTE)
							out << '\t' << std::fixed << std::setprecision(5) << (double)(pop_div[i]) / num_sites;
						else
							out << '\t' << std::fixed << std::setprecision(5) << (double)(pop_div[i] + num_snps[i]) / num_sites;
					}
					else if (dist == "jc")
					{
						if (flag & BAM_SUBSTITUTE)
							pdist = (double)(pop_div[i]) / num_sites;
						else
							pdist = (double)(pop_div[i] + num_snps[i]) / num_sites;
						jc = -0.75 * log(1.0 - pdist * (4.0 / 3.0));
						out << '\t' << std::fixed << std::setprecision(5) << jc;
					}
					else
					{
						out << "\tFixed[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
						out << "\tSeg[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
						out << "\td[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
					}
				}
				else
				{
					out << "\tFixed[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
					out << "\tSeg[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
					out << "\td[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
				}
			}
			break;
		default:
			break;
	}
	std::cout << out.str() << std::endl;

	return 0;
}

int divergeData::setMinPop_n(void)
{
	int j = 0;

	for (j = 0; j < sm->npops; j++)
		min_pop_n[j] = (unsigned short)pop_nsmpl[j];

	return 0;
}

void divergeData::printUsage(const std::string msg)
{
	std::cerr << msg << std::endl << std::endl;
	std::cerr << "Usage:   popbam diverge [options] <in.bam> [region]" << std::endl;
	std::cerr << std::endl;
	std::cerr << "Options: -i          base qualities are Illumina 1.3+     [ default: Sanger ]" << std::endl;
	std::cerr << "         -h  FILE    Input header file                    [ default: none ]" << std::endl;
	std::cerr << "         -d  STR     distance metric (pdist or jc)        [ default: pdist ]" << std::endl;
	std::cerr << "         -o  INT     analysis option                      [ default: 0 ]" << std::endl;
	std::cerr << "                     0 : output individual divergence" << std::endl;
	std::cerr << "                     1 : population divergence statistics" << std::endl;
	std::cerr << "         -p  STR     sample name of outgroup              [ default: reference ]" << std::endl;
	std::cerr << "         -w  INT     use sliding window of size (kb)" << std::endl;
	std::cerr << "         -k  INT     minimum number of sites in window    [ default: 10 ]" << std::endl;
	std::cerr << "         -n  INT     minimum sample size per population   [ default: all samples present ]" << std::endl;
	std::cerr << "         -t          only count substitutions" << std::endl;
	std::cerr << "         -e          exclude singleton polymorphisms" << std::endl;
	std::cerr << "         -f  FILE    Reference fastA file" << std::endl;
	std::cerr << "         -m  INT     minimum read coverage                [ default: 3 ]" << std::endl;
	std::cerr << "         -x  INT     maximum read coverage                [ default: 255 ]" << std::endl;
	std::cerr << "         -q  INT     minimum rms mapping quality          [ default: 25 ]" << std::endl;
	std::cerr << "         -s  INT     minimum snp quality                  [ default: 25 ]" << std::endl;
	std::cerr << "         -a  INT     minimum map quality                  [ default: 13 ]" << std::endl;
	std::cerr << "         -b  INT     minimum base quality                 [ default: 13 ]" << std::endl;
	std::cerr << std::endl;
	exit(EXIT_FAILURE);
}
