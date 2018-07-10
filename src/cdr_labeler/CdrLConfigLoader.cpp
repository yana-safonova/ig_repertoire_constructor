#include "CdrLConfigLoader.hpp"

namespace cdr_labeler {
    std::string CdrLConfigLoader::GetDefaultCfgFilename() const {
        return "configs/cdr_labeler/config.info";
    }

    void CdrLConfigLoader::FillConfigFromCommandline(CDRLabelerConfig&, int, const char* const*) const {}

    std::string CdrLConfigLoader::GetOutputDirPath(const CDRLabelerConfig& config) const {
        return config.output_params.output_dir;
    }
}