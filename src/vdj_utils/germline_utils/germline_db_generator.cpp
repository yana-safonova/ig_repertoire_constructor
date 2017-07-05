#include <logger/logger.hpp>

#include "germline_db_generator.hpp"

namespace germline_utils {
    bool LociParam::LociIncludeIg(std::string loci) {
        if(loci.size() < 2)
            return false;
        return loci == "all" or loci.substr(0, 2) == "IG";
    }

    bool LociParam::LociIncludeTr(std::string loci) {
        if(loci.size() < 2)
            return false;
        return loci == "all" or loci.substr(0, 2) == "TR";
    }

    bool LociParam::LociIncludeIgh(std::string loci) {
        return loci == "all" or loci == "IGH" or loci == "IG";
    }

    bool LociParam::LociIncludeIgk(std::string loci) {
        return loci == "all" or loci == "IGK" or loci == "IG";
    }

    bool LociParam::LociIncludeIgl(std::string loci) {
        return loci == "all" or loci == "IGL" or loci == "IG";
    }

    bool LociParam::LociIncludeTra(std::string loci) {
        return loci == "all" or loci == "TRA" or loci == "TR";
    }

    bool LociParam::LociIncludeTrb(std::string loci) {
        return loci == "all" or loci == "TRB" or loci == "TR";
    }

    bool LociParam::LociIncludeTrg(std::string loci) {
        return loci == "all" or loci == "TRG" or loci == "TR";
    }

    bool LociParam::LociIncludeTrd(std::string loci) {
        return loci == "all" or loci == "TRD" or loci == "TR";
    }

    std::vector<germline_utils::ChainType> LociParam::ConvertIntoChainTypes(std::string loci) {
        std::vector<germline_utils::ChainType> chain_types;
        if (loci.size() < 2)
            return chain_types;
        if (LociIncludeIgh(loci))
            chain_types.push_back(germline_utils::ChainType("IGH"));
        if (LociIncludeIgk(loci))
            chain_types.push_back(germline_utils::ChainType("IGK"));
        if (LociIncludeIgl(loci))
            chain_types.push_back(germline_utils::ChainType("IGL"));
        if (LociIncludeTra(loci))
            chain_types.push_back(germline_utils::ChainType("TRA"));
        if (LociIncludeTrb(loci))
            chain_types.push_back(germline_utils::ChainType("TRB"));
        if (LociIncludeTrg(loci))
            chain_types.push_back(germline_utils::ChainType("TRG"));
        if (LociIncludeTrd(loci))
            chain_types.push_back(germline_utils::ChainType("TRD"));
        return chain_types;
    }

    class GermlineFilesConfig {
        struct ExtendedImmuneGeneType {
            germline_utils::ImmuneGeneType gene_type;
            bool use_pseudogenes;

            ExtendedImmuneGeneType(germline_utils::ImmuneGeneType gene_type, bool use_pseudogenes) :
                    gene_type(gene_type),
                    use_pseudogenes(use_pseudogenes) { }

            bool operator==(const ExtendedImmuneGeneType& obj) const {
                return gene_type == obj.gene_type and use_pseudogenes == obj.use_pseudogenes;
            }
        };

        struct ExtendedImmuneGeneTypeHasher {
            size_t operator()(const ExtendedImmuneGeneType &obj) const {
                return germline_utils::ImmuneGeneTypeHasher()(obj.gene_type) * std::hash<bool>()(obj.use_pseudogenes);
            }
        };

        std::unordered_map<ExtendedImmuneGeneType, std::string, ExtendedImmuneGeneTypeHasher> gene_type_fname_map_;

        germline_utils::ImmuneGeneType ExtractImmuneGeneType(const std::vector<std::string> &splits) {
            return germline_utils::ImmuneGeneType(splits[0]);
        }

        bool FileContainsPseudogenes(const std::vector<std::string> &splits) {
            return splits[1] == "+";
        }

        std::string ExtractFilename(const std::vector<std::string> &splits) {
            return splits[2];
        }

        void AddRecordToFileMap(std::string line) {
            std::vector<std::string> splits;
            boost::split(splits, line, boost::is_any_of("\t "), boost::token_compress_on);
            VERIFY_MSG(splits.size() == 3, "Line of germline config file has unexpected format: " << line <<
                    ", split size: " << splits.size());
            ExtendedImmuneGeneType ext_immune_gene_type(ExtractImmuneGeneType(splits),
                                                        FileContainsPseudogenes(splits));
            gene_type_fname_map_[ext_immune_gene_type] = ExtractFilename(splits);
        }

    public:
        GermlineFilesConfig(std::string config_fname) {
            path::CheckFileExistenceFATAL(config_fname);
            std::ifstream fhandler(config_fname);
            std::string tmp;
            std::getline(fhandler, tmp);
            while(!fhandler.eof()) {
                std::getline(fhandler, tmp);
                if(tmp == "")
                    continue;
                AddRecordToFileMap(tmp);
            }
            INFO(gene_type_fname_map_.size() << " germline filenames were extracted from " << config_fname);
            fhandler.close();
        }

        std::string GetFilenameByImmuneGeneType(germline_utils::ImmuneGeneType gene_type,
                                                bool use_pseudogenes) const {
            ExtendedImmuneGeneType egene_type(gene_type, use_pseudogenes);
            VERIFY_MSG(gene_type_fname_map_.find(egene_type) != gene_type_fname_map_.end(),
                       "Gene type " << gene_type << " was not found in germline files config");
            return gene_type_fname_map_.at(egene_type);
        }
    };

    class ChainDirectoryParam {
        const germline_utils::GermlineInput &gi_;

    public:
        ChainDirectoryParam(const germline_utils::GermlineInput &gi) :
                gi_(gi) { }

        std::string GetDirByChainType(germline_utils::ChainType chain_type) {
            bool lymphocyte_type_unknown = false;
            if(chain_type.Lymphocyte() == germline_utils::LymphocyteType::BLymphocyte)
                return gi_.ig_dir;
            if(chain_type.Lymphocyte() == germline_utils::LymphocyteType::TLymphocyte)
                return gi_.tcr_dir;
            VERIFY_MSG(lymphocyte_type_unknown, "Directory is unknown for lymphocyte type " << chain_type.Lymphocyte());
            return "";
        }
    };

    void GermlineDbGenerator::GenerateGeneFnames() {
        using namespace germline_utils;
        std::string germline_db = path::append_path(germ_params_.germline_dir, germ_params_.organism);
        chain_types_ = LociParam::ConvertIntoChainTypes(germ_params_.loci);
        GermlineFilesConfig germline_files_config(germ_input_.germline_filenames_config);
        ChainDirectoryParam chain_dir_param(germ_input_);
        for(auto it = chain_types_.begin(); it != chain_types_.end(); it++) {
            std::string lymph_dir = path::append_path(germline_db, chain_dir_param.GetDirByChainType(*it));
            v_genes_fnames_.push_back(path::append_path(lymph_dir,
                                                        germline_files_config.GetFilenameByImmuneGeneType(
                                                                ImmuneGeneType(*it, SegmentType::VariableSegment),
                                                                germ_params_.pseudogenes)));
            if (it->IsVDJ())
                d_genes_fnames_.push_back(path::append_path(lymph_dir,
                                                            germline_files_config.GetFilenameByImmuneGeneType(
                                                                ImmuneGeneType(*it, SegmentType::DiversitySegment),
                                                                germ_params_.pseudogenes)));
            j_genes_fnames_.push_back(path::append_path(lymph_dir,
                                                        germline_files_config.GetFilenameByImmuneGeneType(
                                                                ImmuneGeneType(*it, SegmentType::JoinSegment),
                                                                germ_params_.pseudogenes)));
        }
        INFO(v_genes_fnames_.size() << " V gene segment files will be used for DB: ");
        for(size_t i = 0; i < v_genes_fnames_.size(); i++)
            INFO(chain_types_[i] << ": " << v_genes_fnames_[i]);

        INFO(d_genes_fnames_.size() << " D gene segment files will be used for DB: ");
        for(size_t i = 0; i < d_genes_fnames_.size(); i++)
            INFO(chain_types_[i] << ": " << d_genes_fnames_[i]);

        INFO(j_genes_fnames_.size() << " J gene segment files will be used for DB: ");
        for(size_t i = 0; i < j_genes_fnames_.size(); i++)
            INFO(chain_types_[i] << ": " << j_genes_fnames_[i]);
    }

    germline_utils::CustomGeneDatabase GermlineDbGenerator::GenerateVariableDb() {
        germline_utils::CustomGeneDatabase v_custom_db(germline_utils::SegmentType::VariableSegment);
        for(size_t i = 0; i < v_genes_fnames_.size(); i++)
            v_custom_db.AddDatabase(germline_utils::ImmuneGeneType(chain_types_[i],
                                                                   germline_utils::SegmentType::VariableSegment),
                                    v_genes_fnames_[i]);
        return v_custom_db;
    }

    germline_utils::CustomGeneDatabase GermlineDbGenerator::GenerateDiversityDb() {
        germline_utils::CustomGeneDatabase d_custom_db(germline_utils::SegmentType::DiversitySegment);
        for(size_t i = 0; i < d_genes_fnames_.size(); i++)
            d_custom_db.AddDatabase(germline_utils::ImmuneGeneType(chain_types_[i],
                                                                   germline_utils::SegmentType::DiversitySegment),
                                    d_genes_fnames_[i]);
        return d_custom_db;
    }

    germline_utils::CustomGeneDatabase GermlineDbGenerator::GenerateJoinDb() {
        germline_utils::CustomGeneDatabase j_custom_db(germline_utils::SegmentType::JoinSegment);
        for(size_t i = 0; i < j_genes_fnames_.size(); i++)
            j_custom_db.AddDatabase(germline_utils::ImmuneGeneType(chain_types_[i],
                                                                   germline_utils::SegmentType::JoinSegment),
                                    j_genes_fnames_[i]);
        return j_custom_db;
    }
}