/** \file pop_options.cpp
 *  \brief Functions for parsing runtime options for popbam
 *  \author Daniel Garrigan
 *  \version 0.5
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "pop_global.h"
#include "pop_options.h"
#include "pop_sample.h"
#include "pop_utils.h"
#include "popbam.h"
#include "getopt_pp.h"

popbamOptions::popbamOptions(int argc, char *argv[])
{
    std::vector<std::string> glob_opts;

    // set default parameter values
    flag = 0;
    minSites = 10;
    winSize = 1;
    minRMSQ = 25;
    minSNPQ = 25;
    minDepth = 3;
    maxDepth = 255;
    minMapQ = 13;
    minBaseQ = 13;
    hetPrior = 0.0001;
    dist = "pdist";
    errorCount = 0;

    // get the popbam function and iterate argv
    popFunc = argv[1];
    argv++;
    argc--;

    GetOpt::GetOpt_pp args(argc, argv);

    // get optioned arguments
    args >> GetOpt::Option('f', reffile);
    args >> GetOpt::Option('h', headfile);
    args >> GetOpt::Option('o', output);
    args >> GetOpt::Option('z', hetPrior);
    args >> GetOpt::Option('m', minDepth);
    args >> GetOpt::Option('x', maxDepth);
    args >> GetOpt::Option('q', minRMSQ);
    args >> GetOpt::Option('s', minSNPQ);
    args >> GetOpt::Option('a', minMapQ);
    args >> GetOpt::Option('b', minBaseQ);
    args >> GetOpt::Option('k', minSites);
    args >> GetOpt::Option('n', minPop);
    args >> GetOpt::Option('w', winSize);
    args >> GetOpt::Option('d', dist);

    // get switches
    if (args >> GetOpt::OptionPresent('w'))
        {
            winSize *= KB;
            flag |= BAM_WINDOW;
        }
    if (args >> GetOpt::OptionPresent('h'))
        {
            flag |= BAM_HEADERIN;
        }
    if (args >> GetOpt::OptionPresent('p'))
        {
            flag |= BAM_OUTGROUP;
        }
    if (args >> GetOpt::OptionPresent('n'))
        {
            flag |= BAM_MINPOPSAMPLE;
        }
    if (args >> GetOpt::OptionPresent('t'))
        {
            flag |= BAM_SUBSTITUTE;
        }
    if (args >> GetOpt::OptionPresent('e'))
        {
            flag |= BAM_NOSINGLETONS;
        }
    if (args >> GetOpt::OptionPresent('i'))
        {
            flag |= BAM_ILLUMINA;
        }
    if (args >> GetOpt::OptionPresent('z'))
        {
            flag |= BAM_HETEROZYGOTE;
        }
    if (args >> GetOpt::OptionPresent('v'))
        {
            flag |= BAM_VARIANT;
        }

    // get non-optioned arguments
    args >> GetOpt::GlobalOption(glob_opts);

    // run some checks on the command line

    if ((dist != "pdist") && (dist != "jc"))
        {
            errorMsg = dist + " is not a valid distance option";
            errorCount++;
        }

    // check if output option is valid
    if ((output < 0) || (output > 2))
        {
            errorMsg = "Not a valid output option";
            errorCount++;
        }

    // if no input BAM file is specified -- print usage and exit
    if (glob_opts.size() < 2)
        {
            errorMsg = "Need to specify BAM file name and region";
            errorCount++;
        }
    else
        {
            bamfile = glob_opts[0];
            region = glob_opts[1];
        }

    // check if specified BAM file exists on disk
    if (!(is_file_exist(bamfile.c_str())))
        {
            errorMsg = "Specified input file: " + bamfile + " does not exist";
            errorCount++;
        }

    // check if fastA reference file is specified
    if (reffile.empty())
        {
            errorMsg = "Need to specify fastA reference file";
            errorCount++;
        }
    else if (!(is_file_exist(reffile.c_str())))
        {
            errorMsg = "Specified reference file: " + reffile + " does not exist";
            errorCount++;
        }

    // check if BAM header input file exists on disk
    if (flag & BAM_HEADERIN)
        {
            if (!(is_file_exist(headfile.c_str())))
                {
                    errorMsg = "Specified header file: " + headfile + " does not exist";
                    errorCount++;
                }
        }
}

int
popbamOptions::checkBAM(void)
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
        {
            h = bam_in->header;
        }

    // read in new header text
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
    return 0;
}
