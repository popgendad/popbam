POPBAM
======

POPBAM is a tool to perform evolutionary or population-based analyses of next-generation sequencing data. 
POPBAM takes a BAM file as its input and can compute many widely used evolutionary genetics measures in 
sliding windows across a genome.

INTRODUCTION
------------

This is the fifth beta release (0.5b) of the POPBAM program. The source code has primarily been alpha tested 
on Debian and Red Hat based operating systems using the GNU C++ compiler. A Makefile is provided for the 
GNU C++ or clang compilers with three compilation modes: release (default), debug, and profile.

POPBAM is now being once again actively developed for Lummei Analytics LLC. The code base is
being refactored with the aim of incorporating support for VCF and GFF files.

OBTAINING THE SOURCE CODE
-------------------------

The POPBAM source code is available to download freely from GitLab at

[https://gitlab.com/evolgen/popbam]

The only dependencies of the POPBAM program are the zlib compression library and headers
and htslib which is available in the repositories of some linux distributions or
can be obtained from

[https://github.com/samtools/htslib]

Please insure these are installed on your system before attempting to compile POPBAM.

COMPILING THE SOURCE CODE
-------------------------

Once, the repository is cloned, move into the POPBAM source code directory by 

	cd popbam

Then you can simply type

	make

to build the source code.  If you have administrator privileges, you can automatically install 
the POPBAM executable into the /usr/local/bin directory by typing

	sudo make install

Lastly, you can clean the POPBAM directory by using

	make clean

GETTING HELP USING THE PROGRAM
------------------------------

Currently, the primary resource for helping users run POPBAM is a manpage.
If one has system administrator privileges, the popbam.1 file can be installed
into the system manpath or, alternatively, the POBPAM manpage can be viewed in
the current directory by typing

    man ./popbam.1

EXAMPLE DATA SET
----------------

The example data set mentioned in the Evolutionary Bioinformatics paper is 
no longer being made available.
