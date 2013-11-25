
#include "popbam.h"
#include "tables.h"

template <class T> unsigned long long* callBase(T *t, int n, const bam_pileup1_t *pl)
{
	int i = 0;
	int j = 0;
	int qq = 0;
	int baseQ = 0;
	int tmp_baseQ = 0;
	int b = 0;
	int si = -1;
	int rmsq = 0;
	int n_smpl = t->sm->n;
	unsigned short k = 0;
	unsigned short *bases = nullptr;
	unsigned long long rms = 0;
	unsigned long long *cb;
	std::vector<int> depth(n_smpl, 0);
	unsigned char *s = nullptr;
	float q[16];
	std::string msg;
	bam_pileup1_t ***p = nullptr;
	kstring_t buf;

	// allocate memory pileup data
	try
	{
		cb = new unsigned long long [n_smpl]();
		p = new bam_pileup1_t** [n_smpl];
		for (i = 0; i < n_smpl; i++)
			p[i] = new bam_pileup1_t* [t->max_depth];
	}
	catch (std::bad_alloc& ba)
	{
		std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
	}

	memset(&buf, 0, sizeof(kstring_t));

	// partition pileup according to sample
	for (i = 0; i < n; i++)
	{
		if ((pl+i)->is_del || (pl+i)->is_refskip || ((pl+i)->b->core.flag & BAM_FUNMAP))
			continue;
		s = bam_aux_get((pl+i)->b, "RG");

		// skip reads with no read group tag
		if (!s)
			continue;
		else
			si = bam_smpl_rg2smid(t->sm, t->bamfile.c_str(), (char*)(s+1), &buf);
		if (si < 0)
			si = bam_smpl_rg2smid(t->sm, t->bamfile.c_str(), 0, &buf);

		if (si < 0)
		{
			free(buf.s);
			std::string rogue_rg(bam_aux2Z(s));
			msg = "Problem assigning read group " + rogue_rg + " to a sample.\nPlease check BAM header for correct SM and PO tags";
			fatalError(msg);
		}

		if (depth[si] < t->max_depth)
		{
			p[si][depth[si]] = const_cast<bam_pileup1_t*>(pl+i);
			depth[si]++;
		}
		else
			continue;
	}

	// fill in the base array
	for (j = 0; j < n_smpl; ++j)
	{
		rmsq = 0;
		if (depth[j] > 0)
		{
			try
			{
				bases = new unsigned short [depth[j]]();
			}
			catch (std::bad_alloc& ba)
			{
				std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
			}

			for (i=k=0; i < depth[j]; ++i)
			{
				tmp_baseQ = bam1_qual(p[j][i]->b)[p[j][i]->qpos];

				if (t->flag & BAM_ILLUMINA)
					baseQ = tmp_baseQ > 31 ? tmp_baseQ - 31 : 0;
				else
					baseQ = tmp_baseQ;

				assert(baseQ >= 0);

				if ((baseQ < t->min_baseQ) || (p[j][i]->b->core.qual < t->min_mapQ))
					continue;

				b = bam_nt16_nt4_table[bam1_seqi(bam1_seq(p[j][i]->b), p[j][i]->qpos)];

				if (b > 3)
					continue;

				qq = baseQ < p[j][i]->b->core.qual ? baseQ : p[j][i]->b->core.qual;

				if (qq < 4)
					qq = 4;

				if (qq > 63)
					qq = 63;

				bases[k++] = qq << 5 | (unsigned short)bam1_strand(p[j][i]->b) << 4 | b;
				rmsq += SQ(p[j][i]->b->core.qual);
			}

			// calculate genotype likelihoods
			errmod_cal(t->em, k, NBASES, bases, q);

			// finalize root mean quality score
			rms = (unsigned long long)(sqrt((float)(rmsq) / k) + 0.499);

			// get consensus base call
			cb[j] = gl2cns(q, k);

			// add root-mean map quality score to cb array
			cb[j] |= rms << (CHAR_BIT * 6);

			// take out some garbage
			delete [] bases;
			bases = NULL;
		}
		else
			continue;
	}

	// take out garbage
	for (i = 0; i < n_smpl; i++)
		delete [] p[i];
	delete [] p;
	free(buf.s);

	return cb;
}