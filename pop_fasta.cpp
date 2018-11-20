/** \file pop_fasta.cpp
 *  \brief Functions for creating a multi-entry fasta from a BAM file
 *  \author Daniel Garrigan
 *  \version 0.5
*/

#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include "pop_global.h"
#include "pop_options.h"
#include "pop_sample.h"
#include "pop_utils.h"
#include "popbam.h"
#include "pop_fasta.h"

void usageFasta(const std::string);

int
mainFasta(int argc, char *argv[])
{
    usageFasta("Function not yet implemented");
    return 0;
}

void
usageFasta(const std::string msg)
{
    std::cerr << msg << std::endl;
    std::cerr << std::endl;
    std::cerr << "Usage:   popbam fasta [options] <in.bam> [region]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Options: -i          base qualities are Illumina 1.3+               [ default: Sanger ]" << std::endl;
    std::cerr << "         -h  FILE    Input header file                              [ default: none ]" << std::endl;
    std::cerr << "         -z  FLT     output heterozygous base calls                 [ default: consensus ]" << std::endl;
    std::cerr << "         -f  FILE    Reference fastA file" << std::endl;
    std::cerr << "         -m  INT     minimum read coverage                          [ default: 3 ]" << std::endl;
    std::cerr << "         -x  INT     maximum read coverage                          [ default: 255 ]" << std::endl;
    std::cerr << "         -q  INT     minimum rms mapping quality                    [ default: 25 ]" << std::endl;
    std::cerr << "         -s  INT     minimum snp quality                            [ default: 25 ]" << std::endl;
    std::cerr << "         -a  INT     minimum map quality                            [ default: 13 ]" << std::endl;
    std::cerr << "         -b  INT     minimum base quality                           [ default: 13 ]" << std::endl << std::endl;
    exit(EXIT_FAILURE);
}
