#pragma once

#include "cdr_config.hpp"
#include "../config/config_utils.hpp"

namespace cdr_labeler {
    class CdrLConfigLoader : public config_utils::ConfigLoader<CDRLabelerConfig> {
    protected:
        std::string GetDefaultCfgFilename() const override;

        void FillConfigFromCommandline(CDRLabelerConfig&, int, const char* const*) const override;

        std::string GetOutputDirPath(const CDRLabelerConfig& config) const override;
    };
}
