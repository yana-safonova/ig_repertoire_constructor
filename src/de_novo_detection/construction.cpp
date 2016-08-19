

namespace dnd {

void construct(std::string const& filename_in, std::string const& output_directory,
               bool use_forward_reads, bool use_rc_reads) {
    auto& params = pog::pog_parameters::instance;

    params.set_kmer_size(1);
    pog::partial_order_graph graph = pog::from_file(filename_in, use_forward_reads, use_rc_reads);
    pog::save_graph(output_directory + "/")
}

} // namespace dnd
