/** \file pop_nucdiv.cpp
 *  \brief Functions for calculating nucleotide diversity statistics
 *  \author Daniel Garrigan
 *  \version 0.4
*/

#include "pop_nucdiv.h"
#include "tables.h"

int main_nucdiv(int argc, char *argv[])
{
	int k;
	int chr;                  //! chromosome identifier
	int beg;                  //! beginning coordinate for analysis
	int end;                  //! end coordinate for analysis
	int ref;                  //! ref
	long num_windows;         //! number of windows
	long cw;                  //! counter for windows
	std::string msg;          //! string for error message
	bam_plbuf_t *buf;         //! pileup buffer
	nucdivData *t;            //! nucdiv data structure

	// allocate memory for nucdiv data structre
	t = new nucdivData;

	// parse the command line options
	std::string region = t->parseCommandLine(argc, argv);

	// check input BAM file for errors
	t->checkBAM();

	// initialize the sample data structure
	t->bam_smpl_init();

	// add samples
	t->bam_smpl_add();

	// initialize error model
	t->em = errmod_init(1.-0.83);

	// parse genomic region
	k = bam_parse_region(t->h, region, &chr, &beg, &end);
	if (k < 0)
	{
		msg = "Bad genome coordinates: " + region;
		fatal_error(msg, __FILE__, __LINE__, 0);
	}

	// fetch reference sequence
	t->ref_base = faidx_fetch_seq(t->fai_file, t->h->target_name[chr], 0, 0x7fffffff, &(t->len));

	// calculate the number of windows
	if (t->flag & BAM_WINDOW)
		num_windows = ((end - beg) - 1) / t->win_size;
	else
	{
		t->win_size = end - beg;
		num_windows = 1;
	}

	// iterate through all windows along specified genomic region
	for (cw=0; cw < num_windows; cw++)
	{
		// construct genome coordinate string
		std::string scaffold_name(t->h->target_name[chr]);
		std::ostringstream winc(scaffold_name);
		winc.seekp(0, std::ios::end);
		winc << ":" << beg + (cw * t->win_size) + 1 << "-" << ((cw + 1) * t->win_size) + (beg - 1);
		std::string winCoord = winc.str();

		// initialize number of sites to zero
		t->num_sites = 0;

		// parse the BAM file and check if region is retrieved from the reference
		if (t->flag & BAM_WINDOW)
		{
			k = bam_parse_region(t->h, winCoord, &ref, &(t->beg), &(t->end));
			if (k < 0)
			{
				msg = "Bad window coordinates " + winCoord;
				fatal_error(msg, __FILE__, __LINE__, 0);
			}
		}
		else
		{
			ref = chr;
			t->beg = beg;
			t->end = end;
			if (ref < 0)
			{
				msg = "Bad scaffold name: " + region;
				fatal_error(msg, __FILE__, __LINE__, 0);
			}
		}

		// initialize nucdiv variables
		t->init_nucdiv();

		// create population assignments
		t->assign_pops();

		// initialize pileup
		buf = bam_plbuf_init(make_nucdiv, t);

		// fetch region from bam file
		if ((bam_fetch(t->bam_in->x.bam, t->idx, ref, t->beg, t->end, buf, fetch_func)) < 0)
		{
			msg = "Failed to retrieve region " + region + " due to corrupted BAM index file";
			fatal_error(msg, __FILE__, __LINE__, 0);
		}

		// finalize pileup
		bam_plbuf_push(0, buf);

		// calculate nucleotide diversity in window
		t->calc_nucdiv();

		// print results to stdout
		t->print_nucdiv(chr);

		// take out the garbage
		t->destroy_nucdiv();
		bam_plbuf_destroy(buf);
	}
	// end of window iteration

	errmod_destroy(t->em);
	samclose(t->bam_in);
	bam_index_destroy(t->idx);
	t->bam_smpl_destroy();
	free(t->ref_base);
	delete t;

	return 0;
}

int make_nucdiv(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
{
	int i;
	int fq;
	unsigned long long sample_cov;
	unsigned long long *cb = NULL;
	nucdivData *t = NULL;

	// get control data structure
	t = (nucdivData*)data;

	// only consider sites located in designated region
	if ((t->beg <= (int)pos) && (t->end > (int)pos))
	{
		// allocate memory pileup data
		try
		{
			cb = new unsigned long long [t->sm->n]();
		}
		catch (std::bad_alloc& ba)
		{
			std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
		}

		// call bases
		t->call_base(n, pl, cb);

		// resolve heterozygous sites
		if (!(t->flag & BAM_HETEROZYGOTE))
			clean_heterozygotes(t->sm->n, cb, (int)t->ref_base[pos], t->min_snpQ);

		// determine if site is segregating
		fq = segbase(t->sm->n, cb, t->ref_base[pos], t->min_snpQ);

		// determine how many samples pass the quality filters
		sample_cov = qfilter(t->sm->n, cb, t->min_rmsQ, t->min_depth, t->max_depth);

		// determine population coverage
		for (i=0; i < t->sm->npops; ++i)
		{
			unsigned long long pc = 0;
			pc = sample_cov & t->pop_mask[i];
			unsigned int ncov = bitcount64(pc);
			unsigned int req = (unsigned int)((t->min_pop * t->pop_nsmpl[i]) + 0.4999);
			if (ncov >= req)
				t->pop_cov[t->num_sites] |= 0x1U << i;
		}

		// record site type if the site is variable
		if (t->pop_cov[t->num_sites] > 0)
		{
			t->num_sites++;
			if (fq > 0)
				t->types[t->segsites++] = t->cal_site_type(cb);
		}

		// take out the garbage
		delete [] cb;
	}

	return 0;
}

void nucdivData::calc_nucdiv(void)
{
	int i, j, k;
	unsigned int sum;
	unsigned short **freq;

	freq = new unsigned short* [sm->npops];
	for (i=0; i < sm->npops; i++)
		freq[i] = new unsigned short [segsites];

	// calculate the number of aligned sites within and between populations
	for (i=0; i < num_sites; ++i)
	{
		for (j=0; j < sm->npops; ++j)
			if (CHECK_BIT(pop_cov[i],j))
				++ns_within[j];
		for (j=0; j < sm->npops-1; ++j)
			for (k=j+1; k < sm->npops; ++k)
				if (CHECK_BIT(pop_cov[i],j) && CHECK_BIT(pop_cov[i],k))
					++ns_between[UTIDX(sm->npops,j,k)];
	}

	// calculate within population heterozygosity
	for (i=0; i < sm->npops; i++)
	{
		num_snps[i] = 0;
		sum = 0;
		for (j=0; j < segsites; j++)
		{
			unsigned long long pop_type = types[j] & pop_mask[i];
			freq[i][j] = bitcount64(pop_type);

			// calculate within population heterozygosity
			if (((flag & BAM_NOSINGLETONS) && (freq[i][j] > 1)) || !(flag & BAM_NOSINGLETONS))
				sum += 2 * freq[i][j] * (pop_nsmpl[i] - freq[i][j]);
		}
		//TODO: This uses a constant sample size for the bias correction, instead of per site sample size
		if (pop_nsmpl[i] > 1)
			piw[i] = (double)(sum) / (ns_within[i] * SQ(pop_nsmpl[i]-1));
		else
			piw[i] = 0.0;
	}

	// calculate between population heterozygosity
	// this still will always include singletons
	for (i=0; i < sm->npops-1; i++)
	{
		for (j=i+1; j < sm->npops; j++)
		{
			sum = 0;
			for (k=0; k < segsites; k++)
				sum += freq[i][k] * (pop_nsmpl[j] - freq[j][k]) + freq[j][k] * (pop_nsmpl[i] - freq[i][k]);
			pib[UTIDX(sm->npops,i,j)] = (double)(sum) / (ns_between[UTIDX(sm->npops,i,j)] * pop_nsmpl[i] * pop_nsmpl[j]);
		}
	}

	for (i=0; i < sm->npops; i++)
		delete [] freq[i];
	delete [] freq;
}


void nucdivData::print_nucdiv(int chr)
{
	int i, j;

	std::cout << h->target_name[chr] << "\t" << beg+1 << "\t" << end+1;

	for (i=0; i < sm->npops; i++)
	{
		std::cout << "\tns[" << sm->popul[i] << "]:";
		std::cout << "\t" << ns_within[i];
		if (ns_within[i] >= min_sites)
		{
			std::cout << "\tpi[" << sm->popul[i] << "]:";
			std::cout << "\t" << std::fixed << std::setprecision(5) << piw[i];
		}
		else
			std::cout << "\tpi[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
	}

	for (i=0; i < sm->npops-1; i++)
	{
		for (j=i+1; j < sm->npops; j++)
		{
			std::cout << "\tns[" <<  sm->popul[i] << "-" << sm->popul[j] << "]:";
			std::cout << "\t" << ns_between[UTIDX(sm->npops,i,j)];
			if (ns_between[UTIDX(sm->npops, i, j)] >= (unsigned long int)(win_size * KB * min_sites))
			{
				std::cout << "\tdxy[" << sm->popul[i] << "-" << sm->popul[j] << "]:";
				std::cout << "\t" << std::fixed << std::setprecision(5) << pib[UTIDX(sm->npops,i,j)];
			}
			else
				std::cout << "\tdxy[" << sm->popul[i] << "-" << sm->popul[j] << "]:\t" << std::setw(7) << "NA";
		}
	}
	std::cout << std::endl;
}

std::string nucdivData::parseCommandLine(int argc, char *argv[])
{
	std::vector<std::string> glob_opts;
	std::string msg;
#ifdef _MSC_VER
	struct _stat finfo;
#else
	struct stat finfo;
#endif

	//read command line options
	GetOpt::GetOpt_pp args(argc, argv);
	args >> GetOpt::Option('f', reffile);
	args >> GetOpt::Option('h', headfile);
	args >> GetOpt::Option('m', min_depth);
	args >> GetOpt::Option('x', max_depth);
	args >> GetOpt::Option('q', min_rmsQ);
	args >> GetOpt::Option('s', min_snpQ);
	args >> GetOpt::Option('a', min_mapQ);
	args >> GetOpt::Option('b', min_baseQ);
	args >> GetOpt::Option('k', min_sites);
	args >> GetOpt::Option('n', min_pop);
	args >> GetOpt::Option('w', win_size);
	if (args >> GetOpt::OptionPresent('w'))
	{
		win_size *= 1000;
		flag |= BAM_WINDOW;
	}
	if (args >> GetOpt::OptionPresent('h'))
		flag |= BAM_HEADERIN;
	if (args >> GetOpt::OptionPresent('p'))
		flag |= BAM_OUTGROUP;
	if (args >> GetOpt::OptionPresent('i'))
		flag |= BAM_ILLUMINA;
	if (args >> GetOpt::OptionPresent('e'))
		flag |= BAM_NOSINGLETONS;
	args >> GetOpt::GlobalOption(glob_opts);

	// run some checks on the command line

	// if no input BAM file is specified -- print usage and exit
	if (glob_opts.size() < 2)
	{
		msg = "Need to specify BAM file name";
		fatal_error(msg, __FILE__, __LINE__, &nucdivUsage);
	}
	else
		bamfile = glob_opts[0];

	// check if specified BAM file exists on disk
	if ((stat(bamfile.c_str(), &finfo)) != 0)
	{
		msg = "Specified input file: " + bamfile + " does not exist";
		switch(errno)
		{
			case ENOENT:
				std::cerr << "File not found" << std::endl;
				break;
			case EINVAL:
				std::cerr << "Invalid parameter to stat" << std::endl;
				break;
			default:
				std::cerr << "Unexpected error in stat" << std::endl;
				break;
		}
		fatal_error(msg, __FILE__, __LINE__, 0);
	}

	// check if fastA reference file is specified
	if (reffile.empty())
	{
		msg = "Need to specify fasta reference file";
		fatal_error(msg, __FILE__, __LINE__, &nucdivUsage);
	}

	// check is fastA reference file exists on disk
	if ((stat(reffile.c_str(), &finfo)) != 0)
	{
		switch(errno)
		{
			case ENOENT:
				std::cerr << "File not found" << std::endl;
				break;
			case EINVAL:
				std::cerr << "Invalid parameter to stat" << std::endl;
				break;
			default:
				std::cerr << "Unexpected error in stat" << std::endl;
				break;
		}
		msg = "Specified reference file: " + reffile + " does not exist";
		fatal_error(msg, __FILE__, __LINE__, 0);
	}

	//check if BAM header input file exists on disk
	if (flag & BAM_HEADERIN)
	{
		if ((stat(headfile.c_str(), &finfo)) != 0)
		{
			switch(errno)
			{
				case ENOENT:
					std::cerr << "File not found" << std::endl;
					break;
				case EINVAL:
					std::cerr << "Invalid parameter to stat" << std::endl;
					break;
				default:
					std::cerr << "Unexpected error in stat" << std::endl;
					break;
			}
			msg = "Specified header file: " + headfile + " does not exist";
			fatal_error(msg, __FILE__, __LINE__, 0);
		}
	}

	// return the index of first non-optioned argument
	return glob_opts[1];
}

nucdivData::nucdivData(void)
{
	derived_type = NUCDIV;
	min_sites = 0.5;
	win_size = 0;
	min_pop = 1.0;
}

void nucdivData::init_nucdiv(void)
{
	int length = end-beg;
	int npops = sm->npops;

	segsites = 0;

	try
	{
		types = new unsigned long long [length]();
		ns_within = new unsigned long [npops] ();
		ns_between = new unsigned long [npops*(npops-1)] ();
		pop_mask = new unsigned long long [npops]();
		pop_cov = new unsigned int [length]();
		pop_nsmpl = new unsigned char [npops]();
		piw = new double [npops]();
		pib = new double [npops*(npops-1)]();
		num_snps = new int [npops]();
	}
	catch (std::bad_alloc& ba)
	{
		std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
	}
}

void nucdivData::destroy_nucdiv(void)
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
}

void nucdivData::nucdivUsage(void)
{
	std::cerr << std::endl;
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
