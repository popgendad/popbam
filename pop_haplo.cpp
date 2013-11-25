/** \file pop_haplo.cpp
 *  \brief Functions for calculating haplotype-based statistics
 *  \author Daniel Garrigan
 *  \version 0.4
*/

#include "pop_haplo.h"

int main_haplo(int argc, char *argv[])
{
	int k = 0;
	int chr = 0;                  //! chromosome identifier
	int beg = 0;                  //! beginning coordinate for analysis
	int end = 0;                  //! end coordinate for analysis
	int ref = 0;                  //! ref
	long num_windows = 0;         //! number of windows
	long cw = 0;                  //! window counter
	std::string msg;              //! string for error message
	std::string region;           //! the scaffold/chromosome region string
	bam_plbuf_t *buf = nullptr;   //! pileup buffer
	haploData *t = nullptr;       //! pointer to the function data structure

	// allocate memory for nucdiv data structre
	t = new haploData;

	// parse the command line options
	region = t->parseCommandLine(argc, argv);

	// check input BAM file for errors
	t->checkBAM();

	// initialize the sample data structure
	t->bam_smpl_init();

	// add samples
	t->bam_smpl_add();

	// initialize error model
	t->em = errmod_init(0.17);

	// parse genomic region
	k = bam_parse_region(t->h, region, &chr, &beg, &end);
	if (k < 0)
	{
		msg = "Bad genome coordinates: " + region;
		fatalError(msg);
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
				fatalError(msg);
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
				fatalError(msg);
			}
		}

		// initialize diverge specific variables
		t->init_haplo();

		// create population assignments
		t->assign_pops();

		// initialize pileup
		buf = bam_plbuf_init(make_haplo, t);

		// fetch region from bam file
		if ((bam_fetch(t->bam_in->x.bam, t->idx, ref, t->beg, t->end, buf, fetch_func)) < 0)
		{
			msg = "Failed to retrieve region " + region + " due to corrupted BAM index file";
			fatalError(msg);
		}

		// finalize pileup
		bam_plbuf_push(0, buf);

		// calculate haplotype-based statistics
		t->calc_haplo();

		// print results to stdout
		t->print_haplo(chr);

		// take out the garbage
		t->destroy_haplo();
		bam_plbuf_destroy(buf);
	}
	// end of window interation

	errmod_destroy(t->em);
	samclose(t->bam_in);
	bam_index_destroy(t->idx);
	t->bam_smpl_destroy();
	free(t->ref_base);
	delete t;

	return 0;
}

int make_haplo(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
{
	int i = 0;
	int j = 0;
	int fq = 0;
	unsigned long long sample_cov = 0;
	unsigned long long *cb = nullptr;
	haploData *t = nullptr;

	// get control data structure
	t = (haploData*)data;

	// only consider sites located in designated region
	if ((t->beg <= (int)pos) && (t->end > (int)pos))
	{
		// call bases
		cb = callBase(t, n, pl);

		// resolve heterozygous sites
		if (!(t->flag & BAM_HETEROZYGOTE))
			cleanHeterozygotes(t->sm->n, cb, (int)t->ref_base[pos], t->min_snpQ);

		// determine if site is segregating
		fq = segBase(t->sm->n, cb, t->ref_base[pos], t->min_snpQ);

		// determine how many samples pass the quality filters
		sample_cov = qualFilter(t->sm->n, cb, t->min_rmsQ, t->min_depth, t->max_depth);
		
		// calculate the site type
		t->types[t->segsites] = calculateSiteType(t->sm->n, cb);

		// Update the aligned sites and difference matrices
		for (i=0; i < t->sm->n - 1; i++)
		{
			for (j=i+1; j < t->sm->n; j++)
			{
				if (CHECK_BIT(sample_cov,i) && CHECK_BIT(sample_cov,j))
				{
					t->nsite_matrix[UTIDX(t->sm->n,i,j)]++;
					if ((CHECK_BIT(t->types[t->segsites],i) && !(CHECK_BIT(t->types[t->segsites],j))) || (!(CHECK_BIT(t->types[t->segsites],i)) && CHECK_BIT(t->types[t->segsites],j)))
						t->diff_matrix[UTIDX(t->sm->n,i,j)]++;
				}
			}
		}
		t->segsites++;

		// take out the garbage
		delete [] cb;
	}
	return 0;
}

void haploData::calc_haplo(void)
{
	haplo_func do_haplo[3] = {&haploData::calc_nhaps, &haploData::calc_EHHS, &haploData::calc_Gmin};
	(this->*do_haplo[output])();
}

int haploData::calc_nhaps(void)
{
	int i = 0;
	int j = 0;
	int k = 0;
	int f = 0;
	int nelem = 0;

	for (i=0; i < sm->npops; i++)
	{
		nelem = pop_nsmpl[i];
		if (nelem > 1)
		{
			std::vector<int> b;
			std::vector<int>::iterator it1;
			std::vector<int>::iterator it2;

			// initialize the haplotype identity vector
			for (j=0; j < sm->n; j++)
				if (CHECK_BIT(pop_mask[i], j))
					b.push_back(j);

			// assign haplotype identifiers
			for (j=0, it1=b.begin(); j < nelem-1; j++, it1++)
			{
				it2 = b.begin();
				for (k=j+1, std::advance(it2, j+1); k < nelem; k++, it2++)
					if ((diff_matrix[UTIDX(sm->n,j,k)] == 0) && (*it2 > *it1))
						b.at(k)=j;
			}

			// count number of haplotypes and calculate haplotype diversity
			int ff = 0;
			double sh = 0.0;

			for (j=0; j < (int)b.size(); j++)
			{
				if ((f = count(b.begin(), b.end(), j)) > 0)
					++nhaps[i];
				ff += SQ(f);
			}
			sh = (double)(ff) / SQ(nelem);
			hdiv[i] = 1.0 - ((1.0 - sh) * (double)(nelem / (nelem - 1)));
		}
		else
		{
			nhaps[i] = 1;
			hdiv[i] = 1.0;
		}
	}

	return 0;
}

int haploData::calc_EHHS(void)
{
	int i = 0;
	int j = 0;

	calc_nhaps();

	for (i=0; i < sm->npops; i++)
	{
		if (pop_nsmpl[i] < 4)
			ehhs[i] = std::numeric_limits<double>::quiet_NaN();
		else
		{
			unsigned short popf;
			unsigned long long pop_type;
			std::list<unsigned long long> pop_site;

			// make list container of all non-singleton partitions present in population i
			for (j=0; j < segsites; j++)
			{
				pop_type = types[j] & pop_mask[i];
				popf = bitcount64(pop_type);
				if ((popf > 1) && (popf < (pop_nsmpl[i] - 1)))
					pop_site.push_back(pop_type);
			}

			// count unique partitions
			int before = 0;
			int after = 0;
			int part_count = 0;
			int part_max_count = 0;
			unsigned long long part_type = 0;
			unsigned long long part_type_comp = 0;
			unsigned long long max_site = 0;
			double sh = 0.0;
			std::list<unsigned long long> uniq(pop_site);
			std::list<unsigned long long>::iterator it;
			uniq.sort();
			uniq.unique();

			for (it=uniq.begin(); it != uniq.end(); it++)
			{
				part_type = *it;

				// find the complement of part_type
				for (j=0; j < sm->n; j++)
					if (~CHECK_BIT(part_type,j) && CHECK_BIT(pop_mask[i],j))
						part_type_comp |= 0x1ULL << j;
				before = static_cast<int>(pop_site.size());
				pop_site.remove(part_type);
				pop_site.remove(part_type_comp);
				after = static_cast<int>(pop_site.size());
				part_count = (before - after) + 1;
				if (part_count > part_max_count)
				{
					part_max_count = part_count;
					max_site = part_type;
				}
			}

			// calculate site heterozygosity
			popf = bitcount64(max_site);
			sh = (1.0 - ((double)(SQ(popf) + ((pop_nsmpl[i] - popf) * (pop_nsmpl[i] - popf))) / SQ(pop_nsmpl[i]))) * (double)(pop_nsmpl[i] / (pop_nsmpl[i] - 1));

			// calculate site-specific extended haplotype homozygosity
			ehhs[i] = hdiv[i] / (1.0 - sh);
		}
	}

	return 0;
}

int haploData::calc_Gmin(void)
{
	int i = 0;
	int j = 0;
	int v = 0;
	int w = 0;
	const int npops = sm->npops;
	const int n = sm->n;

	for (i=0; i < npops; i++)
	{
		for (j=i+1; j < npops; j++)
		{
			minDxy[UTIDX(npops,i,j)] = std::numeric_limits<unsigned int>::max();
			for (v=0; v < n-1; v++)
			{
				for (w=v+1; w < n; w++)
				{
					if (CHECK_BIT(pop_mask[i],v) && CHECK_BIT(pop_mask[j],w))
					{
						pib[UTIDX(npops,i,j)] += (double)(diff_matrix[UTIDX(n,v,w)]);
						minDxy[UTIDX(npops,i,j)] = minDxy[UTIDX(npops,i,j)] < diff_matrix[UTIDX(n,v,w)] ? minDxy[UTIDX(npops,i,j)] : diff_matrix[UTIDX(n,v,w)];
					}
				}
			}
			pib[UTIDX(npops,i,j)] *= 1.0 / (double)(pop_nsmpl[i] * pop_nsmpl[j]);
		}
	}

	return 0;
}

void haploData::print_haplo(int chr)
{
	int i = 0;
	int j = 0;

	//print coordinate information and number of aligned sites
	std::cout << h->target_name[chr] << "\t" << beg + 1 << "\t" << end+1 << "\t" << num_sites;

	switch(output)
	{
	case 0:
		for (i=0; i < sm->npops; i++)
		{
			if (num_sites >= min_sites)
			{
				std::cout << "\tK[" << sm->popul[i] << "]:\t" << nhaps[i];
				std::cout << "\tKdiv[" << sm->popul[i] << "]:";
				std::cout << "\t" << std::fixed << std::setprecision(5) << 1.0 - hdiv[i];
			}
			else
			{
				std::cout << "\tK[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
				std::cout << "\tKdiv[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
			}
		}
		break;
	case 1:
		for (i=0; i < sm->npops; i++)
		{
			if (num_sites >= min_sites)
			{
				if (isnan(ehhs[i]))
					std::cout << "\tEHHS[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
				else
				{
					std::cout << "\tEHHS[" << sm->popul[i] << "]:";
					std::cout << "\t" << std::fixed << std::setprecision(5) << ehhs[i];
				}
			}
			else
				std::cout << "\tEHHS[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
		}
		break;
	case 2:
		for (i=0; i < sm->npops; i++)
		{
			if (num_sites >= min_sites)
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
				if (num_sites >= min_sites)
				{
					std::cout << "\tdxy[" << sm->popul[i] << "-" << sm->popul[j] << "]:";
					std::cout << "\t" << std::fixed << std::setprecision(5) << pib[UTIDX(sm->npops,i,j)];
					std::cout << "\tmin[" << sm->popul[i] << "-" << sm->popul[j] << "]:";
					std::cout << "\t" << minDxy[UTIDX(sm->npops,i,j)];
				}
				else
				{
					std::cout << "\tdxy[" << sm->popul[i] << "-" << sm->popul[j] << "]:\t" << std::setw(7) << "NA";
					std::cout << "\tmin[" << sm->popul[i] << "-" << sm->popul[j] << "]:\t" << std::setw(7) << "NA";
				}
			}
		}
		break;
	default:
		break;
	}
	std::cout << std::endl;
}

std::string haploData::parseCommandLine(int argc, char *argv[])
{
#ifdef _MSC_VER
	struct _stat finfo;
#else
	struct stat finfo;
#endif
	std::vector<std::string> glob_opts;
	std::string msg;

	GetOpt::GetOpt_pp args(argc, argv);
	args >> GetOpt::Option('f', reffile);
	args >> GetOpt::Option('h', headfile);
	args >> GetOpt::Option('o', output);
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
		win_size *= KB;
		flag |= BAM_WINDOW;
	}
	if (args >> GetOpt::OptionPresent('h'))
		flag |= BAM_HEADERIN;
	if (args >> GetOpt::OptionPresent('i'))
		flag |= BAM_ILLUMINA;
	args >> GetOpt::GlobalOption(glob_opts);

	// run some checks on the command line

	// check if output option is valid
	if ((output < 0) || (output > 2))
		printUsage("Not a valid output option");

	// if no input BAM file is specified -- print usage and exit
	if (glob_opts.size() < 2)
		printUsage("Need to specify BAM file name");
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
		fatalError(msg);
	}

	// check if fastA reference file is specified
	if (reffile.empty())
		printUsage("Need to specify fastA reference file");

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
		fatalError(msg);
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
			fatalError(msg);
		}
	}

	// return the index of first non-optioned argument
	return glob_opts[1];
}

haploData::haploData(void)
{
	derived_type = HAPLO;
	output = 0;
	min_sites = 0.5;
	min_pop = 1.0;
	win_size = 0;
}

void haploData::init_haplo(void)
{
	const int length = end - beg;
	const int npairs = BINOM(sm->n);
	const int npops = sm->npops;

	segsites = 0;

	try
	{
		types = new unsigned long long [length]();
		pop_mask = new unsigned long long [npops]();
		pop_nsmpl = new unsigned char [npops]();
		nhaps = new int [npops]();
		hdiv = new double [npops]();
		piw = new double [npops]();
		pib = new double [npops*(npops-1)]();
		ehhs = new double [npops]();
		minDxy = new unsigned int [npops*(npops-1)]();
		diff_matrix = new unsigned int [npairs];
		nsite_matrix = new unsigned int [npairs];
	}
	catch (std::bad_alloc& ba)
	{
		std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
	}
}

void haploData::destroy_haplo(void)
{
	delete [] pop_mask;
	delete [] types;
	delete [] pop_nsmpl;
	delete [] hdiv;
	delete [] nhaps;
	delete [] piw;
	delete [] pib;
	delete [] ehhs;
	delete [] minDxy;
	delete [] diff_matrix;
	delete [] nsite_matrix;
}

void haploData::printUsage(std::string msg)
{
	std::cerr << msg << std::endl << std::endl;
	std::cerr << "Usage:   popbam haplo [options] <in.bam> [region]" << std::endl;
	std::cerr << std::endl;
	std::cerr << "Options: -i          base qualities are Illumina 1.3+               [ default: Sanger ]" << std::endl;
	std::cerr << "         -h  FILE    Input header file                              [ default: none ]" << std::endl;
	std::cerr << "         -w  INT     use sliding window of size (kb)" << std::endl;
	std::cerr << "         -k  FLT     minimum proportion of sites covered in window  [ default: 0.5 ]" << std::endl;
	std::cerr << "         -n  FLT     minimum proportion of population covered       [ default: 1.0 ]" << std::endl;
	std::cerr << "         -o  INT     analysis to output                             [ default: 0 ]" << std::endl;
	std::cerr << "                     0 : number of haplotypes" << std::endl;
	std::cerr << "                     1 : extended haplotype homozygosity statistic" << std::endl;
	std::cerr << "                     2 : Gmin statistic" << std::endl;
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
