#include "mutation_map_writer.hpp"

namespace antevolo {
    std::string MutationMapWriter::EdgesToString(std::vector<std::pair<size_t, size_t> > edges) const {
        std::stringstream ss;
        for(size_t i = 0; i < edges.size() - 1; i++)
            ss << edges[i].first << "->" << edges[i].second << ";";
        ss << edges[edges.size() - 1].first << "->" << edges[edges.size() - 1].second;
        return ss.str();
    }

    std::string MutationMapWriter::TripletPairToString(TreeSHM shm) const {
        std::stringstream ss;
        ss << shm.src_triplet << "->" << shm.dst_triplet;
        return ss.str();
    }

    void MutationMapWriter::operator()(std::string output_fname) const {
        std::ofstream out(output_fname);
        out << "#V_gene_name:" << map_.VGeneName() << std::endl;
        out << "#vnPos\tvaPos\tAA\tREG\tMult\tHS\tSC\tSyn\tEdges\tTri" << std::endl;
        for(auto it = map_.shm_mult_cbegin(); it != map_.shm_mult_cend(); it++) {
            //if(!map_.SHMIsNonTrivial(it->first))
            //    continue;
            auto shm = it->first;
            auto edges = map_.GetEdgesBySHM(it->first);
            out << shm.gene_pos << "\t" << shm.gene_aa_pos << "\t" <<
                shm.dst_aa << "\t" << shm.region << "\t" << it->second << "\t" <<
                map_.SHMIsHotSpot(shm) << "\t" << shm.ToStopCodon() << "\t" << shm.Synonymous() << "\t" <<
                EdgesToString(edges) << "\t" << TripletPairToString(shm) << std::endl;
        }
        out.close();
    }
}