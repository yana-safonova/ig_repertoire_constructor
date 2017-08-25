#pragma once

#include <string>
#include <path_helper.hpp>
#include <config_singl.hpp>
#include <copy_file.hpp>

namespace config_utils {
    template<typename Config>
    class ConfigLoader {
    public:
        virtual std::string LoadConfig(int argc, const char* const* argv) const {
            std::string cfg_filename = GetConfigFilename(argc, argv);
            path::CheckFileExistenceFATAL(cfg_filename);
            ConfigHelper::create_instance(cfg_filename);
            FillConfigFromCommandline(ConfigHelper::get_writable(), argc, argv);
            PrepareOutputDir(GetOutputDirPath(ConfigHelper::get()));
            CopyConfigs(cfg_filename, ConfigHelper::get());
            return cfg_filename;
        }

    protected:
        typedef config_common::config<Config> ConfigHelper;

        virtual std::string GetDefaultCfgFilename() const = 0;

        virtual std::string GetConfigFilename(int argc, const char* const* argv) const {
            if(argc == 2 and (std::string(argv[1]) != "--help" and std::string(argv[1]) != "--version" and
                              std::string(argv[1]) != "--help-hidden"))
                return std::string(argv[1]);
            return GetDefaultCfgFilename();
        }

        virtual void FillConfigFromCommandline(Config& config, int argc, const char* const* argv) const = 0;

        virtual std::string GetOutputDirPath(const Config& config) const = 0;

        virtual void PrepareOutputDir(const std::string& output_dir_path) const {
            path::make_dir(output_dir_path);
        }

        virtual void CopyConfigs(const std::string& cfg_filename, const Config& config) const {
            std::string to_dir = path::append_path(GetOutputDirPath(config), "configs");
            path::make_dir(to_dir);
            path::copy_files_by_ext(path::parent_path(cfg_filename), to_dir, ".info", true);
            path::copy_files_by_ext(path::parent_path(cfg_filename), to_dir, ".properties", true);
        }
    };

}