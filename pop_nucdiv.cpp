/** \file pop_nucdiv.cpp
 *  \brief Functions for calculating nucleotide diversity statistics
 *  \author Daniel Garrigan
 *  \version 0.4
*/

#include "pop_nucdiv.h"

int mainNucdiv(int argc, char *argv[])
{
	int chr = 0;                  //! chromosome identifier
	int beg = 0;                  //! beginning coordinate for analysis
	int end = 0;                  //! end coordinate for analysis
	int ref = 0;                  //! ref
	long nWindows = 0;            //! number of windows
	std::string msg;              //! string for error message
	bam_sample_t *sm = nullptr;   //! Pointer to the sample information for the input BAM file
	bam_plbuf_t *buf;             //! pileup buffer

	// initialize user command line options
	popbamOptions p(argc, argv);

	if (p.errorCount > 0)
		usageNucdiv(p.errorMsg);

	// check input BAM file for errors
	p.checkBAM();

	// initialize the sample data structure
	sm = bam_smpl_init();

	// add samples
	bam_smpl_add(sm, &p);

	// initialize the nucdiv data structre
	nucdivData t(p);
	t.sm = sm;

	// initialize error model
	t.em = errmod_init(0.17);

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

		// initialize nucdiv variables
		t.allocNucdiv();

		// create population assignments
		t.assignPops(&p);

		// initialize pileup
		buf = bam_plbuf_init(makeNucdiv, &t);

		// fetch region from bam file
		if ((bam_fetch(p.bam_in->x.bam, p.idx, ref, t.beg, t.end, buf, fetch_func)) < 0)
		{
			msg = "Failed to retrieve region " + p.region + " due to corrupted BAM index file";
			fatalError(msg);
		}

		// finalize pileup
		bam_plbuf_push(0, buf);

		// calculate nucleotide diversity in window
		t.calcNucdiv();

		// print results to stdout
		t.printNucdiv(p.h->target_name[chr]);

		// take out the garbage
		bam_plbuf_destroy(buf);
	}
	// end of window iteration

	errmod_destroy(t.em);
	samclose(p.bam_in);
	bam_index_destroy(p.idx);
	bam_smpl_destroy(sm);
	free(t.ref_base);

	return 0;
}

int makeNucdiv(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
{
	int i = 0;
	int j = 0;
	int fq = 0;
	unsigned long long sample_cov = 0;
	unsigned long long *cb = nullptr;
	nucdivData *t = nullptr;

	// get control data structure
	t = (nucdivData*)data;

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

		// record site type if the site is variable
		if (t->pop_cov[t->num_sites] > 0)
		{
			t->num_sites++;
			if (fq > 0)
			{
				for (j = 0; j < t->sm->npops; ++j)
					t->ncov[j][t->segsites] = ncov[j];
				t->types[t->segsites++] = calculateSiteType(t->sm->n, cb);
			}
		}

		// take out the garbage
		delete [] cb;
		delete [] ncov;
	}

	return 0;
}

int nucdivData::calcNucdiv(void)
{
	int i = 0;
	int j = 0;
	int k = 0;
	unsigned short **freq = nullptr;
	unsigned long long pop_type = 0;
	double sum = 0.0;

	freq = new unsigned short* [sm->npops];
	for (i = 0; i < sm->npops; i++)
		freq[i] = new unsigned short [segsites];

	// calculate the number of aligned sites within and between populations
	for (i = 0; i < num_sites; ++i)
	{
		for (j = 0; j < sm->npops; ++j)
			if (CHECK_BIT(pop_cov[i],j))
				++ns_within[j];
		for (j = 0; j < sm->npops - 1; ++j)
			for (k = j + 1; k < sm->npops; ++k)
				if (CHECK_BIT(pop_cov[i],j) && CHECK_BIT(pop_cov[i],k))
					++ns_between[UTIDX(sm->npops,j,k)];
	}

	// calculate within population heterozygosity
	for (i = 0; i < sm->npops; i++)
	{
		num_snps[i] = 0;
		sum = 0.0;
		for (j = 0; j < segsites; j++)
		{
			pop_type = types[j] & pop_mask[i];
			freq[i][j] = bitcount64(pop_type);

			if (ncov[i][j] > 1)
				if (((flag & BAM_NOSINGLETONS) && (freq[i][j] > 1)) || !(flag & BAM_NOSINGLETONS))
					sum += (2.0 * freq[i][j] * (ncov[i][j] - freq[i][j])) / SQ(ncov[i][j]-1);
		}
		piw[i] = sum / ns_within[i];
	}

	// calculate between population heterozygosity
	// this still will always include singletons
	for (i = 0; i < sm->npops - 1; i++)
	{
		for (j = i + 1; j < sm->npops; j++)
		{
			sum = 0.0;
			for (k = 0; k < segsites; k++)
				if ((ncov[i][k] > 0) && (ncov[j][k] > 0))
					sum += (double)(freq[i][k] * (ncov[j][k] - freq[j][k]) + freq[j][k] * (ncov[i][k] - freq[i][k])) / (ncov[i][k] * ncov[j][k]);
			pib[UTIDX(sm->npops,i,j)] = sum / ns_between[UTIDX(sm->npops,i,j)];
		}
	}

	// take out the garbage
	for (i = 0; i < sm->npops; i++)
		delete [] freq[i];
	delete [] freq;

	return 0;
}


int nucdivData::printNucdiv(const std::string scaffold)
{
	int i = 0;
	int j = 0;
	std::stringstream out;

	out << scaffold << '\t' << beg + 1 << '\t' << end + 1;

	for (i = 0; i < sm->npops; i++)
	{
		out << "\tns[" << sm->popul[i] << "]:";
		out << '\t' << ns_within[i];
		if (ns_within[i] >= minSites)
		{
			out << "\tpi[" << sm->popul[i] << "]:";
			out << '\t' << std::fixed << std::setprecision(5) << piw[i];
		}
		else
			out << "\tpi[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
	}

	for (i = 0; i < sm->npops - 1; i++)
	{
		for (j = i + 1; j < sm->npops; j++)
		{
			out << "\tns[" <<  sm->popul[i] << "-" << sm->popul[j] << "]:";
			out << '\t' << ns_between[UTIDX(sm->npops,i,j)];
			if (ns_between[UTIDX(sm->npops,i,j)] >= (unsigned long int)((end - beg) * minSites))
			{
				out << "\tdxy[" << sm->popul[i] << "-" << sm->popul[j] << "]:";
				out << '\t' << std::fixed << std::setprecision(5) << pib[UTIDX(sm->npops,i,j)];
			}
			else
				out << "\tdxy[" << sm->popul[i] << "-" << sm->popul[j] << "]:\t" << std::setw(7) << "NA";
		}
	}

	std::cout << out.str() << std::endl;

	return 0;
}


nucdivData::nucdivData(const popbamOptions &p)
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
	minSites = p.minSites;
	minPop = p.minPop;

	// initialize native variables
	derived_type = NUCDIV;
}

int nucdivData::allocNucdiv(void)
{
	int i = 0;
	int length = end - beg;
	int npops = sm->npops;

	segsites = 0;

	try
	{
		types = new unsigned long long [length]();
		ns_within = new unsigned long [npops] ();
		ns_between = new unsigned long [npops*(npops-1)] ();
		pop_mask = new unsigned long long [npops]();
		pop_cov = new unsigned int [length]();
		ncov = new unsigned int* [npops];
		pop_nsmpl = new unsigned char [npops]();
		piw = new double [npops]();
		pib = new double [npops*(npops-1)]();
		num_snps = new int [npops]();
		for (i=0; i < npops; ++i)
			ncov[i] = new unsigned int [length]();
	}
	catch (std::bad_alloc& ba)
	{
		std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
	}

	return 0;
}

nucdivData::~nucdivData(void)
{
	int npops = sm->npops;

	delete [] pop_mask;
	delete [] types;
	delete [] pop_cov;
	delete [] pop_nsmpl;
	delete [] ns_within;
	delete [] ns_between;
	delete [] piw;
	delete [] pib;
	delete [] num_snps;
	for (int i = 0; i < npops; ++i)
		delete [] ncov[i];
	delete [] ncov;
}

void usageNucdiv(const std::string msg)
{
	std::cerr << msg << std::endl << std::endl;
	std::cerr << "Usage:   popbam nucdiv [options] <in.bam> [region]" << std::endl;
	std::cerr << std::endl;
	std::cerr << "Options: -i          base qualities are Illumina 1.3+               [ default: Sanger ]" << std::endl;
	std::cerr << "         -h  FILE    Input header file                              [ default: none ]" << std::endl;
	std::cerr << "         -w  INT     use sliding window of size (kb)" << std::endl;
	std::cerr << "         -k  FLT     minimum proportion of sites covered in window  [ default: 0.5 ]" << std::endl;
	std::cerr << "         -n  FLT     minimum proportion of population covered       [ default: 1.0 ]" << std::endl;
	std::cerr << "         -e          exclude singleton polymorphisms" << std::endl;
	std::cerr << "         -f  FILE    Reference fastA file" << std::endl;
	std::cerr << "         -m  INT     minimum read coverage                          [ default: 3 ]" << std::endl;
	std::cerr << "         -x  INT     maximum read coverage                          [ default: 255 ]" << std::endl;
	std::cerr << "         -q  INT     minimum rms mapping quality                    [ default: 25 ]" << std::endl;
	std::cerr << "         -s  INT     minimum snp quality                            [ default: 25 ]" << std::endl;
	std::cerr << "         -a  INT     minimum map quality                            [ default: 13 ]" << std::endl;
	std::cerr << "         -b  INT     minimum base quality                           [ default: 13 ]" << std::endl;
	std::cerr << std::endl;
	exit(EXIT_FAILURE);
}
