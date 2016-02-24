/** \file pop_tree.cpp
 *  \brief Functions for constructing phylogenetic trees file from BAM files
 *  \author Daniel Garrigan
 *  \version 0.4
*/

#include <cstdlib>
#include <cstdint>
#include <cfloat>
#include <cassert>
#include <string>
#include <new>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "pop_global.h"
#include "pop_options.h"
#include "pop_sample.h"
#include "pop_utils.h"
#include "popbam.h"
#include "pop_tree.h"

int makeTree(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data);
void usageTree(const std::string);
template uint64_t* callBase<treeData>(treeData *t, int n, const bam_pileup1_t *pl);
void calcDiffMatrix(treeData*);

int
mainTree(int argc, char *argv[])
{
    int chr = 0;                  //! chromosome identifier
    int beg = 0;                  //! beginning coordinate for analysis
    int end = 0;                  //! end coordinate for analysis
    int ref = 0;                  //! ref
    long nWindows = 1;            //! number of windows
    std::string msg;              //! string for error message
    bam_sample_t *sm =
        nullptr;   //! Pointer to the sample information for the input BAM file
    bam_plbuf_t *buf = nullptr;   //! pileup buffer

    // initialize user command line options
    popbamOptions p(argc, argv);

    if (p.errorCount > 0)
        {
            usageTree(p.errorMsg);
        }

    // check input BAM file for errors
    p.checkBAM();

    // initialize the sample data structure
    sm = bam_smpl_init();

    // add samples
    bam_smpl_add(sm, &p);

    // initialize the tree data structure
    treeData t(p);
    t.sm = sm;
    t.npops = sm->npops;

    // initialize error model
    t.em = errmod_init(0.17);

    // extract name of reference sequence
    t.refid = get_refid(p.h->text);

    // parse genomic region
    int k = bam_parse_region(p.h, p.region, &chr, &beg, &end);
    if (k < 0)
        {
            msg = "Bad genome coordinates: " + p.region;
            fatalError(msg);
        }

    // fetch reference sequence
    t.ref_base = faidx_fetch_seq(p.fai_file, p.h->target_name[chr], 0, 0x7fffffff,
                                 &(t.len));

    // calculate the number of windows
    if (p.flag & BAM_WINDOW)
        {
            nWindows = ((end - beg) - 1) / p.winSize;
        }
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
            winc << ':' << beg + (j * p.winSize) + 1 << '-' << ((j + 1) * p.winSize) +
                 (beg - 1);
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

            // initialize tree-specific variables
            t.allocTree();

            // create population assignments
            t.assignPops(&p);

            // initialize pileup
            buf = bam_plbuf_init(makeTree, &t);

            // fetch region from bam file
            if ((bam_fetch(p.bam_in->x.bam, p.idx, ref, t.beg, t.end, buf, fetch_func)) < 0)
                {
                    msg = "Failed to retrieve region " + p.region +
                          " due to corrupted BAM index file";
                    fatalError(msg);
                }

            // finalize pileup
            bam_plbuf_push(0, buf);

            // count pairwise differences
            calcDiffMatrix(&t);

            // construct distance matrix
            t.calcDistMatrix();

            // construct nj tree
            t.makeNJ(std::string(p.h->target_name[chr]));

            // take out the garbage
            bam_plbuf_destroy(buf);
        }  // end of window interation

    errmod_destroy(t.em);
    samclose(p.bam_in);
    bam_index_destroy(p.idx);
    bam_smpl_destroy(sm);
    delete [] t.refid;
    free(t.ref_base);
    return 0;
}

int
makeTree(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data)
{
    int i = 0;
    int fq = 0;
    uint64_t sample_cov = 0;
    uint64_t *cb = nullptr;
    treeData *t = nullptr;

    // get control data structure
    t = (treeData*)data;

    // only consider sites located in designated region
    if ((t->beg <= (int)pos) && (t->end > (int)pos))
        {
            // call bases
            cb = callBase(t, n, pl);

            // resolve heterozygous sites
            if (!(t->flag & BAM_HETEROZYGOTE))
                {
                    cleanHeterozygotes(t->sm->n, cb, (int)t->ref_base[pos], t->minSNPQ);
                }

            // determine if site is segregating
            fq = segBase(t->sm->n, cb, t->ref_base[pos], t->minSNPQ);

            // determine how many samples pass the quality filters
            sample_cov = qualFilter(t->sm->n, cb, t->minRMSQ, t->minDepth, t->maxDepth);

            for (i = 0; i < t->npops; i++)
                {
                    t->pop_sample_mask[i] = sample_cov & t->pop_mask[i];
                }

            if (bitcount64(sample_cov) == t->sm->n)
                {
                    // calculate the site type
                    t->types[t->num_sites] = calculateSiteType(t->sm->n, cb);

                    if (fq > 0)
                        {
                            t->hap.pos[t->segsites] = pos;
                            t->hap.ref[t->segsites] = bam_nt16_table[(int)t->ref_base[pos]];
                            for (i = 0; i < t->sm->n; i++)
                                {
                                    t->hap.rms[i][t->segsites] = (cb[i] >> (CHAR_BIT * 6)) & 0xffff;
                                    t->hap.snpq[i][t->segsites] = (cb[i] >> (CHAR_BIT * 4)) & 0xffff;
                                    t->hap.num_reads[i][t->segsites] = (cb[i] >> (CHAR_BIT * 2)) & 0xffff;
                                    t->hap.base[i][t->segsites] = bam_nt16_table[(int)iupac[(cb[i] >> CHAR_BIT) &
                                                                  0xff]];
                                    if (cb[i] & 0x2ULL)
                                        {
                                            t->hap.seq[i][t->segsites / 64] |= 0x1ULL << t->segsites % 64;
                                        }
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

int
treeData::makeNJ(const std::string scaffold)
{
    int i = 0;
    tree curtree;
    node **cluster = nullptr;

    if ((num_sites < minSites) || (segsites < 1))
        {
            std::cout << scaffold << '\t' << beg + 1 << '\t' << end + 1 << '\t' <<
                      num_sites;
            std::cout << "\tNA" << std::endl;
            return 0;
        }
    try
        {
            cluster = new node* [ntaxa];
            enterorder = new int [ntaxa];
        }
    catch (std::bad_alloc& ba)
        {
            std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
        }
    for (i = 0; i < ntaxa; i++)
        {
            enterorder[i] = i + 1;
        }
    initTree(&curtree.nodep);
    node *p = nullptr;
    p = curtree.nodep[2 * ntaxa - 2]->next;
    curtree.nodep[2 * ntaxa - 2]->next = curtree.nodep[2 * ntaxa - 2];
    free(p->next);
    free(p);
    setupTree(&curtree);
    for (i = 0; i < ntaxa; i++)
        {
            cluster[i] = curtree.nodep[i];
        }
    joinTree(curtree, cluster);
    curtree.start = curtree.nodep[0]->back;
    std::cout << scaffold << '\t' << beg + 1 << '\t' << end + 1 << '\t' << num_sites
              << '\t';
    printTree(curtree.start, curtree.start);
    freeTree(&curtree.nodep);
    delete [] cluster;
    delete [] enterorder;
}

void
treeData::joinTree(tree curtree, node **cluster)
{
    int nc = 0;
    int nextnode = 0;
    int mini = 0;
    int minj = 0;
    int i = 0;
    int j = 0;
    int jj = 0;
    int ii = 0;
    int ia = 0;
    int ja = 0;
    int nude = 0;
    int iter = 0;
    int el[3];
    int *oc = nullptr;
    double fotu2 = 0;
    double total = 0;
    double tmin = 0;
    double dio = 0;
    double djo = 0;
    double bi, bj, bk;
    double dmin = 0;
    double da = 0;
    double *av = nullptr;
    double **x = nullptr;
    double *R = nullptr;

    try
        {
            av = new double [ntaxa];
            oc = new int [ntaxa];
            R = new double [ntaxa];
            x = new double* [ntaxa];
            for (i = 0; i < ntaxa; i++)
                {
                    x[i] = new double [ntaxa];
                    memcpy(x[i], dist_matrix[i], ntaxa * sizeof(double));
                }
        }
    catch (std::bad_alloc& ba)
        {
            std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
        }
    for (i = 0; i < ntaxa - 1; i++)
        {
            for (j = i + 1; j < ntaxa; j++)
                {
                    da = (x[i][j] + x[j][i]) / 2.0;
                    x[i][j] = da;
                    x[j][i] = da;
                }
        }

    // first initialization
    fotu2 = ntaxa - 2.0;
    nextnode = ntaxa + 1;
    for (i = 0; i < ntaxa; i++)
        {
            av[i] = 0.0;
            oc[i] = 1;
        }

    // enter main cycle
    iter = ntaxa - 3;
    for (nc = 1; nc <= iter; nc++)
        {
            for (j = 2; j <= ntaxa; j++)
                {
                    for (i = 0; i <= j - 2; i++)
                        {
                            x[j - 1][i] = x[i][j - 1];
                        }
                }
            tmin = DBL_MAX;

            // compute sij and minimize
            for (i = 0; i < ntaxa; i++)
                {
                    R[i] = 0.0;
                }
            for (ja = 2; ja <= ntaxa; ja++)
                {
                    jj = enterorder[ja - 1];
                    if (cluster[jj - 1] != 0)
                        {
                            for (ia = 0; ia <= ja - 2; ia++)
                                {
                                    ii = enterorder[ia];
                                    if (cluster[ii - 1] != 0)
                                        {
                                            R[ii - 1] += x[ii - 1][jj - 1];
                                            R[jj - 1] += x[ii - 1][jj - 1];
                                        }
                                }
                        }
                }

            for (ja = 2; ja <= ntaxa; ja++)
                {
                    jj = enterorder[ja - 1];
                    if (cluster[jj - 1] != 0)
                        {
                            for (ia = 0; ia <= ja - 2; ia++)
                                {
                                    ii = enterorder[ia];
                                    if (cluster[ii - 1] != 0)
                                        {
                                            total = fotu2 * x[ii - 1][jj - 1] - R[ii - 1] - R[jj - 1];
                                        }
                                    if (total < tmin)
                                        {
                                            tmin = total;
                                            mini = ii;
                                            minj = jj;
                                        }
                                }
                        }
                }
            dio = 0.0;
            djo = 0.0;
            for (i = 0; i < ntaxa; i++)
                {
                    dio += x[i][mini - 1];
                    djo += x[i][minj - 1];
                }
            dmin = x[mini - 1][minj - 1];
            dio = (dio - dmin) / fotu2;
            djo = (djo - dmin) / fotu2;
            bi = (dmin + dio - djo) * 0.5;
            bj = dmin - bi;
            bi -= av[mini - 1];
            bj -= av[minj - 1];
            hookup(curtree.nodep[nextnode-1]->next, cluster[mini - 1]);
            hookup(curtree.nodep[nextnode-1]->next->next, cluster[minj - 1]);
            cluster[mini - 1]->v = bi;
            cluster[minj - 1]->v = bj;
            cluster[mini - 1]->back->v = bi;
            cluster[minj - 1]->back->v = bj;
            cluster[mini - 1] = curtree.nodep[nextnode - 1];
            cluster[minj - 1] = 0;
            nextnode++;
            av[mini - 1] = dmin * 0.5;

            // re-initialization
            fotu2 -= 1.0;
            for (j = 0; j < ntaxa; j++)
                {
                    if (cluster[j] != 0)
                        {
                            da = (x[mini - 1][j] + x[minj - 1][j]) * 0.5;
                            if ((mini - j - 1) < 0)
                                {
                                    x[mini - 1][j] = da;
                                }
                            if ((mini - j - 1) > 0)
                                {
                                    x[j][mini - 1] = da;
                                }
                        }
                }
            for (j = 0; j < ntaxa; j++)
                {
                    x[minj - 1][j] = 0.0;
                    x[j][minj - 1] = 0.0;
                }
            oc[mini - 1] += oc[minj - 1];
        }

    // the last cycle
    nude = 1;
    for (i = 1; i <= ntaxa; i++)
        {
            if (cluster[i - 1] != 0)
                {
                    el[nude - 1] = i;
                    nude++;
                }
        }
    bi = (x[el[0] - 1][el[1] - 1] + x[el[0] - 1][el[2] - 1] - x[el[1] - 1][el[2] - 1]) * 0.5;
    bj = x[el[0] - 1][el[1] - 1] - bi;
    bk = x[el[0] - 1][el[2] - 1] - bi;
    bi -= av[el[0] - 1];
    bj -= av[el[1] - 1];
    bk -= av[el[2] - 1];
    hookup(curtree.nodep[nextnode-1], cluster[el[0] - 1]);
    hookup(curtree.nodep[nextnode-1]->next, cluster[el[1] - 1]);
    hookup(curtree.nodep[nextnode-1]->next->next, cluster[el[2] - 1]);
    cluster[el[0]-1]->v = bi;
    cluster[el[1]-1]->v = bj;
    cluster[el[2]-1]->v = bk;
    cluster[el[0]-1]->back->v = bi;
    cluster[el[1]-1]->back->v = bj;
    cluster[el[2]-1]->back->v = bk;
    curtree.start = cluster[el[0] - 1]->back;

    // take out the garbage
    delete [] av;
    delete [] oc;
    delete [] R;
    for (i = 0; i < ntaxa; i++)
        {
            delete [] x[i];
        }
    delete [] x;
}

void
treeData::hookup(node *p, node *q)
{
    assert(p != 0);
    assert(q != 0);
    p->back = q;
    q->back = p;
}

void
treeData::printTree(node *p, node *start)
{
    if (p->tip)
        {
            if (p->index == 1)
                {
                    std::cout << refid;
                }
            else
                {
                    std::cout << sm->smpl[p->index-2];
                }
        }
    else
        {
            std::cout << '(';
            printTree(p->next->back, start);
            std::cout << ',';
            printTree(p->next->next->back, start);
            if (p == start)
                {
                    std::cout << ',';
                    printTree(p->back, start);
                }
            std::cout << ')';
        }
    if (p == start)
        {
            std::cout << ';' << std::endl;
        }
    else
        {
            if (p->v < 0)
                {
                    std::cout << ":0.00000";
                }
            else
                {
                    std::cout << ':' << std::fixed << std::setprecision(5) << p->v;
                }
        }
}

void
calcDiffMatrix(treeData *t)
{
    int i = 0;
    int j = 0;
    int k = 0;
    int n = t->sm->n;
    int segs = t->segsites;

    // calculate number of differences with reference sequence
    for (i = 0; i < n; i++)
        {
            for (k = 0; k <= SEG_IDX(segs); k++)
                {
                    t->diff_matrix[i+1][0] += bitcount64(t->hap.seq[i][k]);
                }
            t->diff_matrix[0][i+1] = t->diff_matrix[i+1][0];
        }

    // calculate number of pairwise differences
    for (i = 0; i < n - 1; i++)
        {
            for (j = i + 1; j < n; j++)
                {
                    for (k = 0; k <= SEG_IDX(segs); k++)
                        {
                            t->diff_matrix[j + 1][i + 1] += hamming_distance(t->hap.seq[i][k],
                                                        t->hap.seq[j][k]);
                        }
                    t->diff_matrix[i + 1][j + 1] = t->diff_matrix[j + 1][i + 1];
                }
        }
}

int
treeData::calcDistMatrix(void)
{
    int i = 0;
    int j = 0;

    for (i = 0; i < ntaxa - 1; i++)
        {
            for (j = i + 1; j < ntaxa; j++)
                {
                    // p-distance
                    dist_matrix[i][j] = (double)(diff_matrix[i][j]) / num_sites;
                    dist_matrix[j][i] = dist_matrix[i][j];
                    if (dist == "jc")
                        {
                            // Jukes-Cantor distance
                            dist_matrix[i][j] = -0.75 * log(1.0 - (4.0 * dist_matrix[i][j] / 3.0));
                            dist_matrix[j][i] = dist_matrix[i][j];
                        }
                }
        }
}

void
treeData::initTree(ptarray *treenode)
{
    int i = 0;
    int j = 0;
    int nnodes = 2 * ntaxa - 1;
    node *p = nullptr;
    node *q = nullptr;

    *treenode = (ptarray)malloc(nnodes * sizeof(node*));

    for (i = 0; i < ntaxa; i++)
        {
            (*treenode)[i] = (node*)malloc(sizeof(node));
        }
    for (i = ntaxa; i < nnodes; i++)
        {
            q = 0;
            for (j = 1; j <= 3; j++)
                {
                    p = (node*)malloc(sizeof(node));
                    p->next = q;
                    q = p;
                }
            p->next->next->next = p;
            (*treenode)[i] = p;
        }
}

void treeData::setupTree(tree *curtree)
{
    int i = 0;
    int nnodes = 2 * ntaxa - 1;
    node *p = nullptr;

    for (i = 1; i <= nnodes; i++)
        {
            curtree->nodep[i - 1]->back = 0;
            curtree->nodep[i - 1]->tip = (i <= ntaxa);
            curtree->nodep[i - 1]->index = i;
            curtree->nodep[i - 1]->v = 0.0;
            if (i > ntaxa)
                {
                    p = curtree->nodep[i - 1]->next;
                    while (p != curtree->nodep[i - 1])
                        {
                            p->back = 0;
                            p->tip = 0;
                            p->index = i;
                            p = p->next;
                        }
                }
        }
    curtree->start = curtree->nodep[0];
}

void
treeData::freeTree(ptarray *treenode)
{
    int i = 0;
    node *p = nullptr;
    node *q = nullptr;

    for (i = 0; i < ntaxa; i++)
        {
            free((*treenode)[i]);
        }
    for (i = ntaxa; i < 2 * ntaxa - 1; i++)
        {
            p = (*treenode)[i];
            q = p->next;
            while (q != p)
                {
                    node *r = q;
                    q = q->next;
                    free(r);
                }
            free(p);
        }
    free(*treenode);
}

treeData::treeData(const popbamOptions &p)
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
    dist = p.dist;

    // initialize native variables
    derived_type = TREE;
    refid = NULL;
}

int
treeData::allocTree(void)
{
    int i = 0;
    int length = end - beg;
    int n = sm->n;

    ntaxa = n + 1;
    segsites = 0;
    try
        {
            types = new uint64_t [length]();
            pop_mask = new uint64_t [npops]();
            pop_nsmpl = new uint8_t [npops]();
            pop_sample_mask = new uint64_t [npops]();
            hap.pos = new uint32_t [length]();
            hap.idx = new uint32_t [length]();
            hap.ref = new uint8_t [length]();
            hap.seq = new uint64_t* [n];
            hap.base = new uint8_t* [n];
            hap.rms = new uint16_t* [n];
            hap.snpq = new uint16_t* [n];
            hap.num_reads = new uint16_t* [n];
            diff_matrix = new uint16_t* [ntaxa];
            dist_matrix = new double* [ntaxa];
            for (i = 0; i < n; i++)
                {
                    hap.seq[i] = new uint64_t [length]();
                    hap.base[i] = new uint8_t [length]();
                    hap.rms[i] = new uint16_t [length]();
                    hap.snpq[i] = new uint16_t [length]();
                    hap.num_reads[i] = new uint16_t [length]();
                }
            for (i = 0; i < ntaxa; i++)
                {
                    diff_matrix[i] = new uint16_t [ntaxa]();
                    dist_matrix[i] = new double [ntaxa]();
                }
        }
    catch (std::bad_alloc& ba)
        {
            std::cerr << "bad_alloc caught: " << ba.what() << std::endl;
        }
    return 0;
}

treeData::~treeData(void)
{
    int i = 0;

    delete [] pop_mask;
    delete [] types;
    delete [] pop_nsmpl;
    delete [] pop_sample_mask;
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
    for (i = 0; i < ntaxa; i++)
        {
            delete [] diff_matrix[i];
            delete [] dist_matrix[i];
        }
    delete [] hap.seq;
    delete [] hap.base;
    delete [] hap.snpq;
    delete [] hap.rms;
    delete [] hap.num_reads;
    delete [] diff_matrix;
    delete [] dist_matrix;
}

void
usageTree(const std::string msg)
{
    std::cerr << msg << std::endl << std::endl;
    std::cerr << "Usage:   popbam tree [options] <in.bam> [region]" << std::endl;
    std::cerr << std::endl;
    std::cerr <<
              "Options: -i          base qualities are Illumina 1.3+     [ default: Sanger ]"
              << std::endl;
    std::cerr <<
              "         -h  FILE    Input header file                    [ default: none ]" <<
              std::endl;
    std::cerr <<
              "         -d  STR     distance (pdist or jc)               [ default: pdist ]"
              << std::endl;
    std::cerr << "         -w  INT     use sliding window of size (kb)" <<
              std::endl;
    std::cerr <<
              "         -k  INT     minimum number of sites in window    [ default: 10 ]" <<
              std::endl;
    std::cerr << "         -f  FILE    Reference fastA file" << std::endl;
    std::cerr <<
              "         -m  INT     minimum read coverage                [ default: 3 ]" <<
              std::endl;
    std::cerr <<
              "         -x  INT     maximum read coverage                [ default: 255 ]" <<
              std::endl;
    std::cerr <<
              "         -q  INT     minimum rms mapping quality          [ default: 25 ]" <<
              std::endl;
    std::cerr <<
              "         -s  INT     minimum snp quality                  [ default: 25 ]" <<
              std::endl;
    std::cerr <<
              "         -a  INT     minimum map quality                  [ default: 13 ]" <<
              std::endl;
    std::cerr <<
              "         -b  INT     minimum base quality                 [ default: 13 ]" <<
              std::endl;
    std::cerr << std::endl;
    exit(EXIT_FAILURE);
}
