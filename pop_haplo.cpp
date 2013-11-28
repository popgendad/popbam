/** \file pop_haplo.cpp
 *  \brief Functions for calculating haplotype-based statistics
 *  \author Daniel Garrigan
 *  \version 0.4
*/

#include "pop_haplo.h"

int mainHaplo(int argc, char *argv[])
{
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

	if (p.errorCount > 0)
		usageHaplo(p.errorMsg);

	// check input BAM file for errors
	p.checkBAM();

	// initialize the sample data structure
	sm = bam_smpl_init();

	// add samples
	bam_smpl_add(sm, &p);

	// initialize the haplo data structure
	haploData t(p);
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
	for (long i = 0; i < nWindows; ++i)
	{
		// construct genome coordinate string
		std::string scaffold_name(p.h->target_name[chr]);
		std::ostringstream winc(scaffold_name);

		winc.seekp(0, std::ios::end);
		winc << ':' << beg + (i * p.winSize) + 1 << '-' << ((i + 1) * p.winSize) + (beg - 1);

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
		t.allocHaplo();

		// create population assignments
		t.assignPops(&p);

		// initialize pileup
		buf = bam_plbuf_init(makeHaplo, &t);

		// fetch region from bam file
		if ((bam_fetch(p.bam_in->x.bam, p.idx, ref, t.beg, t.end, buf, fetch_func)) < 0)
		{
			msg = "Failed to retrieve region " + p.region + " due to corrupted BAM index file";
			fatalError(msg);
		}

		// finalize pileup
		bam_plbuf_push(0, buf);

		// calculate haplotype-based statistics
		t.calcHaplo();

		// print results to stdout
		t.printHaplo(std::string(p.h->target_name[chr]));

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

int makeHaplo(unsigned int tid, unsigned int pos, int n, const bam_pileup1_t *pl, void *data)
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
			cleanHeterozygotes(t->sm->n, cb, (int)t->ref_base[pos], t->minSNPQ);

		// determine if site is segregating
		fq = segBase(t->sm->n, cb, t->ref_base[pos], t->minSNPQ);

		// determine how many samples pass the quality filters
		sample_cov = qualFilter(t->sm->n, cb, t->minRMSQ, t->minDepth, t->maxDepth);

		// determine population coverage
		for (i = 0; i < t->sm->npops; ++i)
		{
			unsigned long long pc = 0;
			pc = sample_cov & t->pop_mask[i];
			unsigned int ncov = bitcount64(pc);
			unsigned int req = (unsigned int)((t->minPop * t->pop_nsmpl[i]) + 0.4999);
			unsigned long long type = calculateSiteType(t->sm->n, cb);
			unsigned short segi = bitcount64(type & t->pop_mask[i]);
			int k = 0;
			if ((ncov == t->pop_nsmpl[i]) && (segi > 0) && (segi < t->pop_nsmpl[i]))
			{
				for (j = 0, k = 0; j < t->sm->n; j++)
				{
					if ((CHECK_BIT(type,j)) && (CHECK_BIT(t->pop_mask[i],j)))
					{
						t->hap[i][k].push_back('1');
						++k;
					}
					else if (!(CHECK_BIT(type,j)) && (CHECK_BIT(t->pop_mask[i],j)))
					{
						t->hap[i][k].push_back('0');
						++k;
					}
				}
			}
		}

		// calculate the site type
		t->types[t->segsites] = calculateSiteType(t->sm->n, cb);

		// Update the aligned sites and difference matrices
		for (i = 0; i < t->sm->n - 1; i++)
		{
			for (j = i + 1; j < t->sm->n; j++)
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

int haploData::calcHaplo(void)
{
	haplo_func do_haplo[3] = {&haploData::calcNhaps, &haploData::calcEHHS, &haploData::calcGmin};
	(this->*do_haplo[output])();

	return 0;
}

int haploData::calcNhaps(void)
{
	// iterate over populations
	for (int i = 0; i < sm->npops; i++)
	{
		std::set<std::string> hapcount;
		std::multiset<std::string> hapfreq;
		double hom = 0.0;
		for (unsigned int j = 0; j < pop_nsmpl[i]; j++)
		{
			hapcount.insert(hap[i][j]);
			hapfreq.insert(hap[i][j]);
		}
		nhaps[i] = hapcount.size();
		std::set<std::string>::iterator it;
		for (it = hapcount.begin(); it != hapcount.end(); ++it)
		{
			int k = hapfreq.count(*it);
			hom += (double)(SQ(k)) / SQ(pop_nsmpl[i]);
		}
		if ((nhaps[i] > 1) || (pop_nsmpl[i] > 1))
			hdiv[i] = 1.0 - ((1.0 - hom) * (double)(pop_nsmpl[i] / (pop_nsmpl[i] - 1)));
		else
			hdiv[i] = 0.0;
	}

	return 0;
}

int haploData::calcEHHS(void)
{
	int i = 0;
	int j = 0;
	unsigned short popf = 0;
	unsigned long long pop_type = 0;
	int before = 0;
	int after = 0;
	int part_count = 0;
	int part_max_count = 0;
	unsigned long long part_type = 0;
	unsigned long long part_type_comp = 0;
	unsigned long long max_site = 0;
	double sh = 0.0;

	calcNhaps();

	for (i = 0; i < sm->npops; i++)
	{
		if (pop_nsmpl[i] < 4)
			ehhs[i] = std::numeric_limits<double>::quiet_NaN();
		else
		{
			std::list<unsigned long long> pop_site;

			// make list container of all non-singleton partitions present in population i
			for (j = 0; j < segsites; j++)
			{
				pop_type = types[j] & pop_mask[i];
				popf = bitcount64(pop_type);
				if ((popf > 1) && (popf < (pop_nsmpl[i] - 1)))
					pop_site.push_back(pop_type);
			}

			// count unique partitions
			std::list<unsigned long long> uniq(pop_site);
			std::list<unsigned long long>::iterator it;
			uniq.sort();
			uniq.unique();

			for (it = uniq.begin(); it != uniq.end(); it++)
			{
				part_type = *it;

				// find the complement of part_type
				for (j = 0; j < sm->n; j++)
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

int haploData::calcGmin(void)
{
	int i = 0;
	int j = 0;
	int v = 0;
	int w = 0;
	const int npops = sm->npops;
	const int n = sm->n;

	for (i = 0; i < npops; i++)
	{
		for (j = i + 1; j < npops; j++)
		{
			minDxy[UTIDX(npops,i,j)] = std::numeric_limits<unsigned int>::max();
			for (v = 0; v < n - 1; v++)
			{
				for (w = v + 1; w < n; w++)
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

int haploData::printHaplo(const std::string scaffold)
{
	int i = 0;
	int j = 0;
	std::stringstream out;

	//print coordinate information and number of aligned sites
	out << scaffold << '\t' << beg + 1 << '\t' << end+1 << '\t' << num_sites;

	switch(output)
	{
	case 0:
		for (i = 0; i < sm->npops; i++)
		{
			if (num_sites >= minSites)
			{
				out << "\tK[" << sm->popul[i] << "]:\t" << nhaps[i];
				out << "\tKdiv[" << sm->popul[i] << "]:";
				out << '\t' << std::fixed << std::setprecision(5) << 1.0 - hdiv[i];
			}
			else
			{
				out << "\tK[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
				out << "\tKdiv[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
			}
		}
		break;
	case 1:
		for (i = 0; i < sm->npops; i++)
		{
			if (num_sites >= minSites)
			{
				if (isnan(ehhs[i]))
					out << "\tEHHS[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
				else
				{
					out << "\tEHHS[" << sm->popul[i] << "]:";
					out << '\t' << std::fixed << std::setprecision(5) << ehhs[i];
				}
			}
			else
				out << "\tEHHS[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
		}
		break;
	case 2:
		for (i = 0; i < sm->npops; i++)
		{
			if (num_sites >= minSites)
			{
				out << "\tpi[" << sm->popul[i] << "]:";
				out << '\t' << std::fixed << std::setprecision(5) << piw[i];
			}
			else
				out << "\tpi[" << sm->popul[i] << "]:\t" << std::setw(7) << "NA";
		}

		for (i = 0; i < sm->npops-1; i++)
		{
			for (j = i + 1; j < sm->npops; j++)
			{
				if (num_sites >= minSites)
				{
					out << "\tdxy[" << sm->popul[i] << "-" << sm->popul[j] << "]:";
					out << '\t' << std::fixed << std::setprecision(5) << pib[UTIDX(sm->npops,i,j)];
					out << "\tmin[" << sm->popul[i] << "-" << sm->popul[j] << "]:";
					out << '\t' << minDxy[UTIDX(sm->npops,i,j)];
				}
				else
				{
					out << "\tdxy[" << sm->popul[i] << "-" << sm->popul[j] << "]:\t" << std::setw(7) << "NA";
					out << "\tmin[" << sm->popul[i] << "-" << sm->popul[j] << "]:\t" << std::setw(7) << "NA";
				}
			}
		}
		break;
	default:
		break;
	}

	// print final output stream
	std::cout << out.str() << std::endl;

	return 0;
}

haploData::haploData(const popbamOptions &p)
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
	output = p.output;

	// initialize native variables
	derived_type = HAPLO;
}

int haploData::allocHaplo(void)
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
		hap.resize(npops);
		for (int i = 0; i < npops; ++i)
		{
			hap[i].resize(pop_nsmpl[i]);
			for (int j = 0; j < pop_nsmpl[i]; ++j)
				hap[i][j] = "";
		}
	}
	catch (std::bad_alloc& ba)
	{
		std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
	}

	return 0;
}

haploData::~haploData(void)
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

void usageHaplo(const std::string msg)
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
