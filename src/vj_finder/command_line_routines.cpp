#include "command_line_routines.hpp"

#include <boost/program_options.hpp>
#include <build_info.hpp>

bool command_line_requires_parsing(int argc, const char* const* argv) {
    if(argc == 1)
        return false;
    if(argc > 2)
        return true;
    return std::string(argv[1]) != "--help" or
           std::string(argv[1]) != "--version" or std::string(argv[1]) != "--help-hidden";
}

// cfg contains default values from config file
void parse_command_line_args(vj_finder::VJFinderConfig &cfg, int argc, const char* const* argv) {
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
            //("input-file,i", po::value<std::string>(&cfg.iop.input.input_reads)->required(),
            ("input-file,i", po::value<std::string>(&cfg.io_params.input_params.input_reads)->default_value(cfg.io_params.input_params.input_reads),
             "name of an input file (FASTA|FASTQ)")
            //("output-dir,o", po::value<std::string>(&cfg.iop.output.of.output_dir)->required(),
            ("output-dir,o", po::value<std::string>(&cfg.io_params.output_params.output_files.output_dir)->default_value(cfg.io_params.output_params.output_files.output_dir),
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
            //("verbose,V", po::value<bool>(&cfg.io_params.output_params.output_details.verbose)->default_value(cfg.io_params.output_params.output_details.verbose)->implicit_value(true),
            // "produce alignment output for each query")
            ("fix-spaces", po::value<bool>(&cfg.io_params.output_params.output_details.fix_spaces)->default_value(cfg.io_params.output_params.output_details.fix_spaces),
             "replace spaces in read headers with underline symbol '_'")

            ("pseudogenes,P", po::value<bool>(&cfg.algorithm_params.germline_params.pseudogenes)->default_value(cfg.algorithm_params.germline_params.pseudogenes),
             "use pseudogenes along with normal germline genes")
            ("loci,l", po::value<std::string>(&cfg.algorithm_params.germline_params.loci)->default_value(cfg.algorithm_params.germline_params.loci),
             "loci: IGH, IGL, IGK, IG (all BCRs), TRA, TRB, TRG, TRD, TR (all TCRs) or all")
            ("organism", po::value<std::string>(&cfg.algorithm_params.germline_params.organism)->default_value(cfg.algorithm_params.germline_params.organism),
             "organism ('human', 'mouse', 'pig', 'rabbit', 'rat' and 'rhesus_monkey' are supported)")

            ("threads,t", po::value<size_t>(&cfg.run_params.num_threads)->default_value(cfg.run_params.num_threads),
             "the number of threads")

            ("word-size,k", po::value<size_t>(&cfg.algorithm_params.aligner_params.word_size_v)->default_value(cfg.algorithm_params.aligner_params.word_size_v),
             "word size for V genes")
            ("word-size-j", po::value<size_t>(&cfg.algorithm_params.aligner_params.word_size_j)->default_value(cfg.algorithm_params.aligner_params.word_size_j),
             "word size for J genes")
            ("min-k-coverage-v,n", po::value<size_t>(&cfg.algorithm_params.aligner_params.min_k_coverage_v)->default_value(cfg.algorithm_params.aligner_params.min_k_coverage_v),
             "minimal block coverage for V gene")
            ("min-k-coverage-j", po::value<size_t>(&cfg.algorithm_params.aligner_params.min_k_coverage_j)->default_value(cfg.algorithm_params.aligner_params.min_k_coverage_j),
             "minimal block coverage for J gene")
            ("max-candidates-v,N", po::value<size_t>(&cfg.algorithm_params.aligner_params.max_candidates_v)->default_value(cfg.algorithm_params.aligner_params.max_candidates_v),
             "maximal number of V gene candidates for each query")
            ("max-candidates-j", po::value<size_t>(&cfg.algorithm_params.aligner_params.max_candidates_j)->default_value(cfg.algorithm_params.aligner_params.max_candidates_j),
             "maximal number of J gene candidates for each query")

            ("min-len", po::value<size_t>(&cfg.algorithm_params.filtering_params.min_aligned_length)->default_value(cfg.algorithm_params.filtering_params.min_aligned_length),
             "minimal length of reported sequence")
            ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden("Hidden options");
    hidden.add_options()
            ("help-hidden", "show all options, including developers options")
            ("left-uncoverage-limit", po::value<int>(&cfg.algorithm_params.filtering_params.left_uncovered_limit)->default_value(cfg.algorithm_params.filtering_params.left_uncovered_limit),
             "uncoverage limit of left end")
            ("right-uncoverage-limit", po::value<int>(&cfg.algorithm_params.filtering_params.right_uncovered_limit)->default_value(cfg.algorithm_params.filtering_params.right_uncovered_limit),
             "uncoverage limit of right end")
            ("min-vsegment-length", po::value<size_t>(&cfg.algorithm_params.filtering_params.min_v_segment_length)->default_value(cfg.algorithm_params.filtering_params.min_v_segment_length),
             "minimal allowed length of V gene segment")
            ("min-jsegment-length", po::value<size_t>(&cfg.algorithm_params.filtering_params.min_j_segment_length)->default_value(cfg.algorithm_params.filtering_params.min_j_segment_length),
             "minimal allowed length of J gene segment")

            ("db-directory", po::value<std::string>(&cfg.algorithm_params.germline_params.germline_dir)->default_value(cfg.algorithm_params.germline_params.germline_dir),
             "directory with germline database")

            ("max-global-gap-v", po::value<int>(&cfg.algorithm_params.scoring_params.v_scoring.max_global_gap)->default_value(cfg.algorithm_params.scoring_params.v_scoring.max_global_gap),
             "maximal allowed size of global gap for V gene")
            ("max-local-insertions-v", po::value<int>(&cfg.algorithm_params.scoring_params.v_scoring.max_local_insertions)->default_value(cfg.algorithm_params.scoring_params.v_scoring.max_local_insertions),
             "maximal allowed size of local insertion for V gene")
            ("max-local-deletions-v", po::value<int>(&cfg.algorithm_params.scoring_params.v_scoring.max_local_deletions)->default_value(cfg.algorithm_params.scoring_params.v_scoring.max_local_deletions),
             "maximal allowed size of local deletion for V gene")
            ("gap-opening-cost-v", po::value<int>(&cfg.algorithm_params.scoring_params.v_scoring.gap_opening_cost)->default_value(cfg.algorithm_params.scoring_params.v_scoring.gap_opening_cost),
             "gap opening cost for V gene alignement")
            ("gap-extention-cost-v", po::value<int>(&cfg.algorithm_params.scoring_params.v_scoring.gap_extention_cost)->default_value(cfg.algorithm_params.scoring_params.v_scoring.gap_extention_cost),
             "gap extention cost for V gene alignment")
            ("mismatch-extention-cost-v", po::value<int>(&cfg.algorithm_params.scoring_params.v_scoring.mismatch_extention_cost)->default_value(cfg.algorithm_params.scoring_params.v_scoring.mismatch_extention_cost),
             "mismatch extention cost for V gene alignment")
            ("mismatch-opening-cost-v", po::value<int>(&cfg.algorithm_params.scoring_params.v_scoring.mismatch_opening_cost)->default_value(cfg.algorithm_params.scoring_params.v_scoring.mismatch_opening_cost),
             "mismatch opening cost for V gene alignment")
            //("min-vsegment-length", po::value<size_t>(&cfg.algorithm_params.asp.v_scoring.min_v_segment_length)->default_value(min_v_segment_length),
            // "minimal allowed length of V gene segment")

            ("max-global-gap-j", po::value<int>(&cfg.algorithm_params.scoring_params.j_scoring.max_global_gap)->default_value(cfg.algorithm_params.scoring_params.j_scoring.max_global_gap),
             "maximal allowed size of global gap for J gene")
            ("max-local-insertions-j", po::value<int>(&cfg.algorithm_params.scoring_params.j_scoring.max_local_insertions)->default_value(cfg.algorithm_params.scoring_params.j_scoring.max_local_insertions),
             "maximal allowed size of local insertion for J gene")
            ("max-local-deletions-j", po::value<int>(&cfg.algorithm_params.scoring_params.j_scoring.max_local_deletions)->default_value(cfg.algorithm_params.scoring_params.j_scoring.max_local_deletions),
             "maximal allowed size of local deletion for J gene")
            ("gap-opening-cost-j", po::value<int>(&cfg.algorithm_params.scoring_params.j_scoring.gap_opening_cost)->default_value(cfg.algorithm_params.scoring_params.j_scoring.gap_opening_cost),
             "gap opening cost for J gene alignment")
            ("gap-extention-cost-j", po::value<int>(&cfg.algorithm_params.scoring_params.j_scoring.gap_extention_cost)->default_value(cfg.algorithm_params.scoring_params.j_scoring.gap_extention_cost),
             "gap extention cost for J gene alignment")
            ("mismatch-extention-cost-j", po::value<int>(&cfg.algorithm_params.scoring_params.j_scoring.mismatch_extention_cost)->default_value(cfg.algorithm_params.scoring_params.j_scoring.mismatch_extention_cost),
             "mismatch extention cost for J gene alignment")
            ("mismatch-opening-cost-j", po::value<int>(&cfg.algorithm_params.scoring_params.j_scoring.mismatch_opening_cost)->default_value(cfg.algorithm_params.scoring_params.j_scoring.mismatch_opening_cost),
             "mismatch opening cost for J gene alignment")
            //("min-jsegment-length", po::value<size_t>(&min_j_segment_length)->default_value(min_j_segment_length),
            // "minimal allowed length of J gene segment")

            ("fill-left", po::value<bool>(&cfg.algorithm_params.fix_crop_fill_params.fill_left)->default_value(cfg.algorithm_params.fix_crop_fill_params.fill_left),
             "fill left cropped positions by germline")
            ("fill-right", po::value<bool>(&cfg.algorithm_params.fix_crop_fill_params.fill_right)->default_value(cfg.algorithm_params.fix_crop_fill_params.fill_right),
             "fill right cropped positions by germline")
            ("crop-left", po::value<bool>(&cfg.algorithm_params.fix_crop_fill_params.crop_left)->default_value(cfg.algorithm_params.fix_crop_fill_params.crop_left),
             "crop extra left positions")
            ("crop-right", po::value<bool>(&cfg.algorithm_params.fix_crop_fill_params.crop_right)->default_value(cfg.algorithm_params.fix_crop_fill_params.crop_right),
             "crop extra right positions")
            ("fix-left", po::value<size_t>(&cfg.algorithm_params.fix_crop_fill_params.fix_left)->default_value(cfg.algorithm_params.fix_crop_fill_params.fix_left),
             "the number left read positions which will be fixed by germline")
            ("fix-right", po::value<size_t>(&cfg.algorithm_params.fix_crop_fill_params.fix_right)->default_value(cfg.algorithm_params.fix_crop_fill_params.fix_right),
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

    std::cout << cfg.io_params.output_params.output_files.output_dir << std::endl;
    update_output_files_config(cfg.io_params.output_params.output_files);
}