
#ifndef __POP_DIVERGE_PARSER_HPP__
#define __POP_DIVERGE_PARSER_HPP__

#include <stdint.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>
#include <string.h>
#include <stdexcept>
#include <string>
#include <limits>
#include <vector>
#include <iostream>
#include <sstream>
#include <memory>

class pop_diverge_parser
{
    static bool adjust_double_si_suffix(double &res, const char *suffix)
    {
        if(*suffix == '\0')
            {
                return true;
            }
        if(*(suffix + 1) != '\0')
            {
                return false;
            }

        switch(*suffix)
            {
            case 'a':
                res *= 1e-18;
                break;
            case 'f':
                res *= 1e-15;
                break;
            case 'p':
                res *= 1e-12;
                break;
            case 'n':
                res *= 1e-9;
                break;
            case 'u':
                res *= 1e-6;
                break;
            case 'm':
                res *= 1e-3;
                break;
            case 'k':
                res *= 1e3;
                break;
            case 'M':
                res *= 1e6;
                break;
            case 'G':
                res *= 1e9;
                break;
            case 'T':
                res *= 1e12;
                break;
            case 'P':
                res *= 1e15;
                break;
            case 'E':
                res *= 1e18;
                break;
            default:
                return false;
            }
        return true;
    }

    static double conv_double(const char *str, ::std::string &err, bool si_suffix)
    {
        char *endptr = 0;
        errno = 0;
        double res = strtod(str, &endptr);
        if(errno)
            {
                err.assign(strerror(errno));
                return (double)0.0;
            }
        bool invalid =
            si_suffix ? !adjust_double_si_suffix(res, endptr) : *endptr != '\0';
        if(invalid)
            {
                err.assign("Invalid character");
                return (double)0.0;
            }
        return res;
    }

    static int conv_enum(const char* str, ::std::string& err,
                         const char* const strs[])
    {
        int res = 0;
        for(const char* const* cstr = strs; *cstr; ++cstr, ++res)
            if(!strcmp(*cstr, str))
                {
                    return res;
                }
        err += "Invalid constant '";
        err += str;
        err += "'. Expected one of { ";
        for(const char* const* cstr = strs; *cstr; ++cstr)
            {
                if(cstr != strs)
                    {
                        err += ", ";
                    }
                err += *cstr;
            }
        err += " }";
        return -1;
    }

    template<typename T>
    static bool adjust_int_si_suffix(T &res, const char *suffix)
    {
        if(*suffix == '\0')
            {
                return true;
            }
        if(*(suffix + 1) != '\0')
            {
                return false;
            }

        switch(*suffix)
            {
            case 'k':
                res *= (T)1000;
                break;
            case 'M':
                res *= (T)1000000;
                break;
            case 'G':
                res *= (T)1000000000;
                break;
            case 'T':
                res *= (T)1000000000000;
                break;
            case 'P':
                res *= (T)1000000000000000;
                break;
            case 'E':
                res *= (T)1000000000000000000;
                break;
            default:
                return false;
            }
        return true;
    }

    template<typename T>
    static T conv_int(const char *str, ::std::string &err, bool si_suffix)
    {
        char *endptr = 0;
        errno = 0;
        long long int res = strtoll(str, &endptr, 0);
        if(errno)
            {
                err.assign(strerror(errno));
                return (T)0;
            }
        bool invalid =
            si_suffix ? !adjust_int_si_suffix(res, endptr) : *endptr != '\0';
        if(invalid)
            {
                err.assign("Invalid character");
                return (T)0;
            }
        if(res > ::std::numeric_limits<T>::max() ||
                res < ::std::numeric_limits<T>::min())
            {
                err.assign("Value out of range");
                return (T)0;
            }
        return (T)res;
    }

    template<typename T>
    static T conv_uint(const char *str, ::std::string &err, bool si_suffix)
    {
        char *endptr = 0;
        errno = 0;
        while(isspace(*str))
            {
                ++str;
            }
        if(*str == '-')
            {
                err.assign("Negative value");
                return (T)0;
            }
        unsigned long long int res = strtoull(str, &endptr, 0);
        if(errno)
            {
                err.assign(strerror(errno));
                return (T)0;
            }
        bool invalid =
            si_suffix ? !adjust_int_si_suffix(res, endptr) : *endptr != '\0';
        if(invalid)
            {
                err.assign("Invalid character");
                return (T)0;
            }
        if(res > ::std::numeric_limits<T>::max())
            {
                err.assign("Value out of range");
                return (T)0;
            }
        return (T)res;
    }

    template<typename T>
    static ::std::string vec_str(const std::vector<T> &vec)
    {
        ::std::ostringstream os;
        for(typename ::std::vector<T>::const_iterator it = vec.begin();
                it != vec.end(); ++it)
            {
                if(it != vec.begin())
                    {
                        os << ",";
                    }
                os << *it;
            }
        return os.str();
    }

    class string : public ::std::string
    {
    public:
        string() : ::std::string() {}
        explicit string(const ::std::string &s) : std::string(s) {}
        explicit string(const char *s) : ::std::string(s) {}
        int as_enum(const char* const strs[])
        {
            ::std::string err;
            int res = conv_enum((const char*)this->c_str(), err, strs);
            if(!err.empty())
                {
                    throw ::std::runtime_error(err);
                }
            return res;
        }


        uint32_t as_uint32_suffix() const
        {
            return as_uint32(true);
        }
        uint32_t as_uint32(bool si_suffix = false) const
        {
            ::std::string err;
            uint32_t res = conv_uint<uint32_t>((const char*)this->c_str(), err, si_suffix);
            if(!err.empty())
                {
                    ::std::string msg("Invalid conversion of '");
                    msg += *this;
                    msg += "' to uint32_t: ";
                    msg += err;
                    throw ::std::runtime_error(msg);
                }
            return res;
        }
        uint64_t as_uint64_suffix() const
        {
            return as_uint64(true);
        }
        uint64_t as_uint64(bool si_suffix = false) const
        {
            ::std::string err;
            uint64_t res = conv_uint<uint64_t>((const char*)this->c_str(), err, si_suffix);
            if(!err.empty())
                {
                    ::std::string msg("Invalid conversion of '");
                    msg += *this;
                    msg += "' to uint64_t: ";
                    msg += err;
                    throw ::std::runtime_error(msg);
                }
            return res;
        }
        int32_t as_int32_suffix() const
        {
            return as_int32(true);
        }
        int32_t as_int32(bool si_suffix = false) const
        {
            ::std::string err;
            int32_t res = conv_int<int32_t>((const char*)this->c_str(), err, si_suffix);
            if(!err.empty())
                {
                    ::std::string msg("Invalid conversion of '");
                    msg += *this;
                    msg += "' to int32_t: ";
                    msg += err;
                    throw ::std::runtime_error(msg);
                }
            return res;
        }
        int64_t as_int64_suffix() const
        {
            return as_int64(true);
        }
        int64_t as_int64(bool si_suffix = false) const
        {
            ::std::string err;
            int64_t res = conv_int<int64_t>((const char*)this->c_str(), err, si_suffix);
            if(!err.empty())
                {
                    ::std::string msg("Invalid conversion of '");
                    msg += *this;
                    msg += "' to int64_t: ";
                    msg += err;
                    throw ::std::runtime_error(msg);
                }
            return res;
        }
        int as_int_suffix() const
        {
            return as_int(true);
        }
        int as_int(bool si_suffix = false) const
        {
            ::std::string err;
            int res = conv_int<int>((const char*)this->c_str(), err, si_suffix);
            if(!err.empty())
                {
                    ::std::string msg("Invalid conversion of '");
                    msg += *this;
                    msg += "' to int_t: ";
                    msg += err;
                    throw ::std::runtime_error(msg);
                }
            return res;
        }
        long as_long_suffix() const
        {
            return as_long(true);
        }
        long as_long(bool si_suffix = false) const
        {
            ::std::string err;
            long res = conv_int<long>((const char*)this->c_str(), err, si_suffix);
            if(!err.empty())
                {
                    ::std::string msg("Invalid conversion of '");
                    msg += *this;
                    msg += "' to long_t: ";
                    msg += err;
                    throw ::std::runtime_error(msg);
                }
            return res;
        }
        double as_double_suffix() const
        {
            return as_double(true);
        }
        double as_double(bool si_suffix = false) const
        {
            ::std::string err;
            double res = conv_double((const char*)this->c_str(), err, si_suffix);
            if(!err.empty())
                {
                    ::std::string msg("Invalid conversion of '");
                    msg += *this;
                    msg += "' to double_t: ";
                    msg += err;
                    throw ::std::runtime_error(msg);
                }
            return res;
        }
    };

public:
    const char *                   ref_arg;
    bool                           ref_given;
    const char *                   head_arg;
    bool                           head_given;
    const char *                   outgroup_arg;
    bool                           outgroup_given;
    const char *                   out_arg;
    bool                           out_given;
    int                            format_arg;
    bool                           format_given;
    int                            min_depth_arg;
    bool                           min_depth_given;
    int                            max_depth_arg;
    bool                           max_depth_given;
    int                            min_rms_arg;
    bool                           min_rms_given;
    int                            min_snp_arg;
    bool                           min_snp_given;
    int                            min_map_arg;
    bool                           min_map_given;
    int                            min_base_arg;
    bool                           min_base_given;
    double                         min_sites_arg;
    bool                           min_sites_given;
    double                         min_pops_arg;
    bool                           min_pops_given;
    bool                           illumina_flag;
    bool                           subst_flag;
    bool                           single_flag;
    bool                           pop_flag;
    bool                           clean_hets;
    double                         win_size_arg;
    bool                           win_size_given;
    const char *                   dist_arg;
    bool                           dist_given;
    const char *                   input_arg;
    const char *                   region_arg;

    enum
    {
        START_OPT = 1000,
        HELP_OPT
    };

    pop_diverge_parser() :
        ref_arg(""), ref_given(false),
        head_arg(""), head_given(false),
        outgroup_arg(""), outgroup_given(false),
        out_arg(""), out_given(false),
        format_arg(0), format_given(false),
        min_depth_arg(3), min_depth_given(false),
        max_depth_arg(255), max_depth_given(false),
        min_rms_arg(25), min_rms_given(false),
        min_snp_arg(25), min_snp_given(false),
        min_map_arg(13), min_map_given(false),
        min_base_arg(13), min_base_given(false),
        min_sites_arg(1.0), min_sites_given(false),
        min_pops_arg(1.0), min_pops_given(false),
        illumina_flag(true),
        subst_flag(true),
        single_flag(true),
        pop_flag(true),
        clean_hets(false),
        win_size_arg(1.0), win_size_given(false),
        dist_arg("pdist"), dist_given(false),
        input_arg(""),
        region_arg("")
    { }

    pop_diverge_parser(int argc, char* argv[]) :
        ref_arg(""), ref_given(false),
        head_arg(""), head_given(false),
        outgroup_arg(""), outgroup_given(false),
        out_arg(""), out_given(false),
        format_arg(0), format_given(false),
        min_depth_arg(3), min_depth_given(false),
        max_depth_arg(255), max_depth_given(false),
        min_rms_arg(25), min_rms_given(false),
        min_snp_arg(25), min_snp_given(false),
        min_map_arg(13), min_map_given(false),
        min_base_arg(13), min_base_given(false),
        min_sites_arg(1.0), min_sites_given(false),
        min_pops_arg(1.0), min_pops_given(false),
        illumina_flag(false),
        subst_flag(false),
        single_flag(false),
        pop_flag(false),
        clean_hets(false),
        win_size_arg(1.0), win_size_given(false),
        dist_arg("pdist"), dist_given(false),
        input_arg(""),
        region_arg("")
    {
        parse(argc, argv);
    }

    void parse(int argc, char* argv[])
    {
        static struct option long_options[] =
        {
            {"ref", 1, 0, 'r'},
            {"head", 1, 0, 'h'},
            {"outgroup", 1, 0, 'g'},
            {"out", 1, 0, 'o'},
            {"format", 1, 0, 'F'},
            {"min-depth", 1, 0, 'm'},
            {"max-depth", 1, 0, 'x'},
            {"min-rms", 1, 0, 'q'},
            {"min-snp", 1, 0, 's'},
            {"min-map", 1, 0, 'a'},
            {"min-base", 1, 0, 'b'},
            {"min-sites", 1, 0, 'k'},
            {"min-pops", 1, 0, 'n'},
            {"illumina", 0, 0, 'i'},
            {"subst", 0, 0, 't'},
            {"single", 0, 0, 'e'},
            {"pop", 0, 0, 'p'},
            {"hom", 0, 0, 'c'},
            {"win-size", 1, 0, 'w'},
            {"dist", 1, 0, 'd'},
            {"help", 0, 0, HELP_OPT},
            {"usage", 0, 0, 'U'},
            {"version", 0, 0, 'V'},
            {0, 0, 0, 0}
        };
        static const char *short_options = "VUF:r:h:g:o:m:x:q:s:a:b:k:n:itepcw:d:";

        ::std::string err;
#define CHECK_ERR(type,val,which) if(!err.empty()) { ::std::cerr << "Invalid " #type " '" << val << "' for [" which "]: " << err << "\n"; exit(1); }
        while(true)
            {
                int index = -1;
                int c = getopt_long(argc, argv, short_options, long_options, &index);
                if(c == -1)
                    {
                        break;
                    }
                switch(c)
                    {
                    case ':':
                        ::std::cerr << "Missing required argument for "
                                    << (index == -1 ? ::std::string(1,
                                            (char)optopt) : std::string(long_options[index].name))
                                    << ::std::endl;
                        exit(1);
                    case HELP_OPT:
                        ::std::cout << usage() << "\n\n" << help() << std::endl;
                        exit(0);
                    case 'U':
                        ::std::cout << usage() << "\nUse --help for more information." << std::endl;
                        exit(0);
                    case 'V':
                        print_version();
                        exit(0);
                    case '?':
                        ::std::cerr << "Use --usage or --help for some help\n";
                        exit(1);
                    case 'r':
                        ref_given = true;
                        ref_arg = optarg;
                        break;
                    case 'h':
                        head_given = true;
                        head_arg = optarg;
                        break;
                    case 'g':
                        outgroup_given = true;
                        outgroup_arg = optarg;
                        break;
                    case 'o':
                        out_given = true;
                        out_arg = optarg;
                        break;
                    case 'F':
                        format_given = true;
                        format_arg = conv_int<int>((const char*)optarg, err, false);
                        CHECK_ERR(int_t, optarg, "-F, --format=int")
                        break;
                    case 'm':
                        min_depth_given = true;
                        min_depth_arg = conv_int<int>((const char*)optarg, err, false);
                        CHECK_ERR(int_t, optarg, "-m, --min-depth=int")
                        break;
                    case 'x':
                        max_depth_given = true;
                        max_depth_arg = conv_int<int>((const char*)optarg, err, false);
                        CHECK_ERR(int_t, optarg, "-x, --max-depth=int")
                        break;
                    case 'q':
                        min_rms_given = true;
                        min_rms_arg = conv_int<int>((const char*)optarg, err, false);
                        CHECK_ERR(int_t, optarg, "-q, --min-rms=int")
                        break;
                    case 's':
                        min_snp_given = true;
                        min_snp_arg = conv_int<int>((const char*)optarg, err, false);
                        CHECK_ERR(int_t, optarg, "-s, --min-snp=int")
                        break;
                    case 'a':
                        min_map_given = true;
                        min_map_arg = conv_int<int>((const char*)optarg, err, false);
                        CHECK_ERR(int_t, optarg, "-a, --min-map=int")
                        break;
                    case 'b':
                        min_base_given = true;
                        min_base_arg = conv_int<int>((const char*)optarg, err, false);
                        CHECK_ERR(int_t, optarg, "-b, --min-base=int")
                        break;
                    case 'k':
                        min_sites_given = true;
                        min_sites_arg = conv_double((const char*)optarg, err, false);
                        CHECK_ERR(double_t, optarg, "-k, --min-sites=double")
                        break;
                    case 'n':
                        min_pops_given = true;
                        min_pops_arg = conv_double((const char*)optarg, err, false);
                        CHECK_ERR(double_t, optarg, "-n, --min-pops=double")
                        break;
                    case 'i':
                        illumina_flag = true;
                        break;
                    case 't':
                        subst_flag = true;
                        break;
                    case 'e':
                        single_flag = true;
                        break;
                    case 'p':
                        pop_flag = true;
                        break;
                    case 'c':
                        clean_hets = true;
                        break;
                    case 'w':
                        win_size_given = true;
                        win_size_arg = conv_double((const char*)optarg, err, false);
                        CHECK_ERR(double_t, optarg, "-w, --win-size=double")
                        break;
                    case 'd':
                        dist_given = true;
                        dist_arg = optarg;
                        break;
                    }
            }

        // Check that required switches are present
        if(!ref_given)
            {
                error("[-r, --ref=path] required switch");
            }

        // Parse arguments
        if(argc - optind != 2)
            {
                error("Requires exactly 2 arguments.");
            }
        input_arg = argv[optind];
        ++optind;
        region_arg = argv[optind];
        ++optind;

        // Check access rights
        if(access(input_arg, R_OK))
            {
                ::std::string err("Argument input, access right (read) failed for file '");
                ((err += input_arg) += "': ") += strerror(errno);
                error(err.c_str());
            }
        if(access(ref_arg, R_OK))
            {
                ::std::string
                err("Switch -r, --ref=path, access right (read) failed for file '");
                ((err += ref_arg) += "': ") += strerror(errno);
                error(err.c_str());
            }
    }
    static const char * usage()
    {
        return "Usage: popbam diverge [options] input:path region:string";
    }
    class error
    {
        int code_;
        std::ostringstream msg_;

        // Select the correct version (GNU or XSI) version of
        // strerror_r. strerror_ behaves like the GNU version of strerror_r,
        // regardless of which version is provided by the system.
        static const char* strerror__(char* buf, int res)
        {
            return res != -1 ? buf : "Invalid error";
        }
        static const char* strerror__(char* buf, char* res)
        {
            return res;
        }
        static const char* strerror_(int err, char* buf, size_t buflen)
        {
            return strerror__(buf, strerror_r(err, buf, buflen));
        }
        struct no_t { };

    public:
        static no_t no;
        error(int code = EXIT_FAILURE) : code_(code) { }
        explicit error(const char* msg, int code = EXIT_FAILURE) : code_(code)
        {
            msg_ << msg;
        }
        error(const std::string& msg, int code = EXIT_FAILURE) : code_(code)
        {
            msg_ << msg;
        }
        error& operator<<(no_t)
        {
            char buf[1024];
            msg_ << ": " << strerror_(errno, buf, sizeof(buf));
            return *this;
        }
        template<typename T>
        error& operator<<(const T& x)
        {
            msg_ << x;
            return (*this);
        }
        ~error()
        {
            ::std::cerr << "Error: " << msg_.str() << "\n"
                        << usage() << "\n"
                        << "Use --help for more information"
                        << ::std::endl;
            exit(code_);
        }
    };
    static const char * help()
    {
        return
            "Tools to perform evolutionary analysis from BAM files\n\n\n" \
            "Calculates either the proportion of aligned sites that differ between each individual genome\n"
            \
            "and the reference sequence, or the Jukes-Cantor distance (Jukes and Cantor 1969). Also outputs\n"
            \
            "the mean divergence per population, or the number of substitutions with the outgroup and the\n"
            \
            "number of segregating sites per population, both of which can be used in the HKA test.\n\n"
            "Options (default value in (), *required):\n"
            " -r, --ref=path                          *Reference fasta file name\n"
            " -h, --head=string                        Text file with BAM header\n"
            " -g, --outgroup=string                    Name of outgroup sample\n"
            " -o, --out=string                         Output file name\n"
            " -F, --format                             Output format (0)\n"
            " -m, --min-depth=int                      Minimum read depth to consider a site (3)\n"
            " -x, --max-depth=int                      Maximum read depth to consider a site (255)\n"
            " -q, --min-rms=int                        Minimum root-mean alignment score (25)\n"
            " -s, --min-snp=int                        Minimum SNP quality score to consider site variable (25)\n"
            " -a, --min-map=int                        Minimum mapping score to consider a site (13)\n"
            " -b, --min-base=int                       Minimum base quality (13)\n"
            " -k, --min-sites=double                   Proportion of sites to consider a window (1.0)\n"
            " -n, --min-pops=double                    Minimum proportion of populations with coverage (1.0)\n"
            " -i, --illumina                           Base qualities are Illumina 1.3+ (false)\n"
            " -t, --subst                              Only count substitutions (false)\n"
            " -e, --single                             Exclude singleton sites (false)\n"
            " -p, --pop                                Output population-based divergence statistics (false)\n"
            " -c, --hom                                Make genotypes homozygous derived (false)\n"
            " -w, --win-size=double                    Size of a sliding window in kilobases (1.0)\n"
            " -d, --dist=string                        Sequence distance metric (pdist)\n"
            " -U, --usage                              Usage\n"
            "     --help                               This message\n"
            " -V, --version                            Version";
    }

    static const char* hidden()
    {
        return "";
    }

    void print_version(::std::ostream &os = std::cout) const
    {
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.0"
#endif
        os << "0.5" << "\n";
    }

    void dump(::std::ostream &os = std::cout)
    {
        os << "ref_given:" << ref_given << " ref_arg:" << ref_arg << "\n";
        os << "head_given:" << head_given << " head_arg:" << head_arg << "\n";
        os << "outgroup_given:" << outgroup_given << " outgroup_arg:" << outgroup_arg << "\n";
        os << "out_given:" << out_given << " out_arg:" << out_arg << "\n";
        os << "format_given:" << format_given << " format_arg: " << format_arg << "\n";
        os << "min_depth_given:" << min_depth_given << " min_depth_arg:" << min_depth_arg << "\n";
        os << "max_depth_given:" << max_depth_given << " max_depth_arg:" << max_depth_arg << "\n";
        os << "min_rms_given:" << min_rms_given << " min_rms_arg:" << min_rms_arg << "\n";
        os << "min_snp_given:" << min_snp_given << " min_snp_arg:" << min_snp_arg << "\n";
        os << "min_map_given:" << min_map_given << " min_map_arg:" << min_map_arg << "\n";
        os << "min_base_given:" << min_base_given << " min_base_arg:" << min_base_arg << "\n";
        os << "min_sites_given:" << min_sites_given << " min_sites_arg:" << min_sites_arg << "\n";
        os << "min_pops_given:" << min_pops_given << " min_pops_arg:" << min_pops_arg << "\n";
        os << "illumina_flag:" << illumina_flag << "\n";
        os << "subst_flag:" << subst_flag << "\n";
        os << "single_flag:" << single_flag << "\n";
        os << "pop_flag:" << pop_flag << "\n";
        os << "clean_hets:" << clean_hets << "\n";
        os << "win_size_given:" << win_size_given << " win_size_arg:" << win_size_arg << "\n";
        os << "dist_given:" << dist_given << " dist_arg:" << dist_arg << "\n";
        os << "input_arg:" << input_arg << "\n";
        os << "region_arg:" << region_arg << "\n";
    }
};
#endif // __POP_DIVERGE_PARSER_HPP__"
