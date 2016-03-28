#include "command_line_routines.hpp"

#include <boost/program_options.hpp>
#include <build_info.hpp>

bool command_line_requires_parsing(int argc, char **argv) {
    if(argc == 1)
        return false;
    if(argc > 2)
        return true;
    return std::string(argv[1]) != "--help" or
           std::string(argv[1]) != "--version" or std::string(argv[1]) != "--help-hidden";
}

// cfg contains default values from config file
void parse_command_line_args(vj_finder::vjf_config &cfg, int argc, char** argv) {
    if(!command_line_requires_parsing(argc, argv))
        return;

    namespace po = boost::program_options;
    using std::exit;

    // Declare a group of options that will be
    // allowed only on command line
    po::options_description generic("Generic options");
    generic.add_options()
            ("version,v", "print version string")
            ("help,h", "produce help message")
            ("config-file,c", "name of a file of a configuration")
            ("input-file,i", po::value<std::string>(&cfg.iop.input.input_reads)->required(),
             "name of an input file (FASTA|FASTQ)")
            ("output-dir,o", po::value<std::string>(&cfg.iop.output.of.output_dir)->required(),
             "output directory");

/*        std::string separator;
        auto separator_notifier = [this](const std::string &arg) {
            if (arg == "comma") {
                separator = ",";
            } else if (arg == "semicolon") {
                separator = ";";
            } else if (arg == "tab" || arg == "tabular") {
                separator = "\t";
            } else {
                separator = arg;
            }
        };
        */

    // Declare a group of options that will be
    // allowed both on command line and in
    // config file
    po::options_description config("Configuration");
    config.add_options()
            ("compress,Z", po::value<bool>(&cfg.iop.output.od.compress)->default_value(cfg.iop.output.od.compress),
             "compress output FASTA files using zlib")
            ("verbose,V", po::value<bool>(&cfg.iop.output.od.verbose)->default_value(cfg.iop.output.od.verbose)->implicit_value(true),
             "produce alignment output for each query")
            ("fix-spaces", po::value<bool>(&cfg.iop.output.od.fix_spaces)->default_value(cfg.iop.output.od.fix_spaces),
             "replace spaces in read ids by underline symbol '_'")
            ("separator", po::value<std::string>(&cfg.iop.output.od.separator)->default_value(cfg.iop.output.od.separator),
             "separator for alignment info file: ','")

            ("pseudogenes,P", po::value<bool>(&cfg.algop.gp.pseudogenes)->default_value(cfg.algop.gp.pseudogenes),
             "use pseudogenes along with normal germline genes")
            ("loci,l", po::value<std::string>(&cfg.algop.gp.loci)->default_value(cfg.algop.gp.loci),
             "loci: IGH, IGL, IGK, IG (all BCRs), TRA, TRB, TRG, TRD, TR (all TCRs) or all")
            ("db-directory", po::value<std::string>(&cfg.algop.gp.db_directory)->default_value(cfg.algop.gp.db_directory),
             "directory with germline database")
            ("organism", po::value<std::string>(&cfg.algop.gp.organism)->default_value(cfg.algop.gp.organism),
             "organism ('human', 'mouse', 'pig', 'rabbit', 'rat' and 'rhesus_monkey' are supported)")


            ("threads,t", po::value<size_t>(&cfg.rp.num_threads)->default_value(cfg.rp.num_threads),
             "the number of threads")

            ("word-size,k", po::value<int>(&cfg.algop.ap.K)->default_value(cfg.algop.ap.K),
             "word size for V-genes")
            ("word-size-j", po::value<int>(&cfg.algop.ap.word_size_j)->default_value(cfg.algop.ap.word_size_j),
             "word size for J-genes")
            ("min-k-coverage,n", po::value<int>(&cfg.algop.ap.min_k_coverage)->default_value(cfg.algop.ap.min_k_coverage),
             "minimal k+-coverage")
            ("min-k-coverage-j", po::value<int>(&cfg.algop.ap.min_k_coverage_j)->default_value(cfg.algop.ap.min_k_coverage_j),
             "minimal k+-coverage for J-gene")

            ("max-candidates-v,N", po::value<int>(&cfg.algop.qp.max_candidates)->default_value(cfg.algop.qp.max_candidates),
             "maximal number of candidates for each query")
            ("max-candidates-j", po::value<int>(&cfg.algop.qp.max_candidates_j)->default_value(cfg.algop.qp.max_candidates_j),
             "maximal number of J-gene candidates for each query")
            ("min-len", po::value<size_t>(&cfg.algop.qp.min_len)->default_value(cfg.algop.qp.min_len),
             "minimal length of reported sequence")
            ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden("Hidden options");
    hidden.add_options()
            ("help-hidden", "show all options, including developers' ones")
            ("left-uncoverage-limit", po::value<size_t>(&cfg.algop.ap.left_uncovered_limit)->default_value(cfg.algop.ap.left_uncovered_limit),
             "uncoverage limit of left end")
            ("right-uncoverage-limit", po::value<size_t>(&cfg.algop.ap.right_uncovered_limit)->default_value(cfg.algop.ap.right_uncovered_limit),
             "uncoverage limit of right end")
            ("min-vsegment-length", po::value<size_t>(&cfg.algop.ap.min_v_segment_length)->default_value(cfg.algop.ap.min_v_segment_length),
             "minimal allowed length of V gene segment")
            ("min-jsegment-length", po::value<size_t>(&cfg.algop.ap.min_j_segment_length)->default_value(cfg.algop.ap.min_j_segment_length),
             "minimal allowed length of J gene segment")

            ("max-global-gap", po::value<int>(&cfg.algop.asp.max_global_gap)->default_value(cfg.algop.asp.max_global_gap),
             "maximal allowed size of global gap")
            ("max-local-insertions", po::value<int>(&cfg.algop.asp.max_local_insertions)->default_value(cfg.algop.asp.max_local_insertions),
             "maximal allowed size of local insertion")
            ("max-local-deletions", po::value<int>(&cfg.algop.asp.max_local_deletions)->default_value(cfg.algop.asp.max_local_deletions),
             "maximal allowed size of local deletion")
            ("gap-opening-cost", po::value<int>(&cfg.algop.asp.gap_opening_cost)->default_value(cfg.algop.asp.gap_opening_cost),
             "gap opening cost")
            ("gap-extention-cost", po::value<int>(&cfg.algop.asp.gap_extention_cost)->default_value(cfg.algop.asp.gap_extention_cost),
             "gap extention cost")

            ("fill-left", po::value<bool>(&cfg.algop.fxp.fill_left)->default_value(cfg.algop.fxp.fill_left),
             "fill left cropped positions by germline")
            ("fill-right", po::value<bool>(&cfg.algop.fxp.fill_right)->default_value(cfg.algop.fxp.fill_right),
             "fill right cropped positions by germline")
            ("crop-left", po::value<bool>(&cfg.algop.fxp.crop_left)->default_value(cfg.algop.fxp.crop_left),
             "crop extra left positions")
            ("crop-right", po::value<bool>(&cfg.algop.fxp.crop_right)->default_value(cfg.algop.fxp.crop_right),
             "crop extra right positions")
            ("fix-left", po::value<size_t>(&cfg.algop.fxp.fix_left)->default_value(cfg.algop.fxp.fix_left),
             "the number left read positions which will be fixed by germline")
            ("fix-right", po::value<size_t>(&cfg.algop.fxp.fix_right)->default_value(cfg.algop.fxp.fix_right),
             "the number right read positions which will be fixed by germline")
            ;

    po::options_description cmdline_options("All command line options");
    cmdline_options.add(generic).add(config).add(hidden);

    po::options_description config_file_options;
    config_file_options.add(config).add(hidden);

    po::options_description visible("Allowed options");
    visible.add(generic).add(config);

    po::positional_options_description p;
    p.add("input-file", 1).add("output-dir", 1);

    po::variables_map vm;
    store(po::command_line_parser(argc, argv).
            options(cmdline_options).positional(p).run(), vm);

    if (vm.count("help-hidden")) {
        INFO(cmdline_options);
        exit(0);
    }

    if (vm.count("help")) {
        INFO(visible);
        exit(0);
    }

    if (vm.count("version")) {
        INFO(boost::format("VJ Finder, an algorithm for alignment against VJ immune genes %s; git version %s") %
             build_info::version % build_info::git_hash7);
        exit(0);
    }

    if (vm.count("config-file")) {
        std::string config_file = vm["config-file"].as<std::string>();
        std::ifstream ifs(config_file.c_str());
        if (!ifs) {
            ERROR("Config file " << config_file << " was not found");
            exit(1);
        } else {
            store(parse_config_file(ifs, config_file_options), vm);
            // reparse cmd line again for update config defaults
            store(po::command_line_parser(argc, argv).
                    options(cmdline_options).positional(p).run(), vm);
        }
    }

    try {
        notify(vm);
    } catch (const po::error &e) {
        ERROR("Command-line parser error: " << e.what());
        exit(1);
    } catch (const std::exception &e) {
        ERROR("Unknown exception: " << e.what());
        exit(1);
    }
}