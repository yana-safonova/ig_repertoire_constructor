//
// Created by Andrew Bzikadze on 4/26/17.
//

#include "custom_gene_chooser.hpp"

namespace ig_simulator {

std::vector<double>
CustomGeneChooser::ReadProbabilities(const std::string& filename,
                                     const germline_utils::CustomGeneDatabase& db)
{
    using boost::tokenizer;
    using boost::escaped_list_separator;
    using Tokenizer = tokenizer<escaped_list_separator<char>>;

    VERIFY(db.cbegin() + 1 == db.cend());
    germline_utils::ImmuneGeneType igtype { *db.cbegin() };
    const germline_utils::ImmuneGeneDatabase& igdb = db.GetConstDbByGeneType(igtype);

    std::ifstream in;
    in.open(filename);
    VERIFY(in.is_open());

    std::vector<double> probs(db.size());
    std::string line;

    std::vector <std::string> parsed_vector;
    while (getline(in, line)) {
        if (line.empty())
            break;
        Tokenizer tokenizer(line);
        parsed_vector.assign(tokenizer.begin(), tokenizer.end());
        assert(parsed_vector.size() == 2);
        size_t index_of_current_gene = igdb.GetIndexByName(parsed_vector.front());
        probs[index_of_current_gene] = std::stod(parsed_vector.back());
        VERIFY(probs[index_of_current_gene] >= 0 and
            probs[index_of_current_gene] <= 1);
    }
    in.close();
    return probs;
}

std::discrete_distribution<size_t>
CustomGeneChooser::GetDistr(const std::string& filename,
                            const germline_utils::CustomGeneDatabase& db)
{
    std::vector<double> probs { ReadProbabilities(filename, db) };
    return { probs.begin(), probs.end() };
}

CustomGeneChooser::CustomGeneChooser(const std::vector<germline_utils::CustomGeneDatabase>& db,
                                     const std::string& v_genes_probs,
                                     const std::string& d_genes_probs,
                                     const std::string& j_genes_probs):
        AbstractVDJGeneChooser(db),
        v_distr(GetDistr(v_genes_probs, db.front())),
        d_distr(),
        j_distr(GetDistr(j_genes_probs, db.back()))
{
    if (db.size() == 3) {
        d_distr = GetDistr(d_genes_probs, db[1]);
    }
}

CustomGeneChooser::CustomGeneChooser(const std::vector<germline_utils::CustomGeneDatabase>& db,
                                     const GeneChooserParams::CustomGeneChooserParams& config):
        CustomGeneChooser(db, config.v_genes_probs, config.d_genes_probs, config.j_genes_probs)
{ }

VDJ_GenesIndexTuple CustomGeneChooser::ChooseGenes() const
{
    VDJ_GenesIndexTuple result(size_t(-1), size_t(-1), size_t(-1));

    VERIFY(v_db_p_ != nullptr);
    std::get<0>(result) = v_distr(MTSingleton::GetInstance());

    if (is_vdj) {
        VERIFY(d_db_p_ != nullptr);
        std::get<1>(result) = d_distr(MTSingleton::GetInstance());
    }

    VERIFY(j_db_p_ != nullptr);
    std::get<2>(result) = j_distr(MTSingleton::GetInstance());

    return result;
}

} // End namespace ig_simulator