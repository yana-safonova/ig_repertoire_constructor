try:
    from .build_info import *
except:
    version = ""
    git_hash = ""
    git_hash7 = ""
    git_refspec = ""
    git_describe = ""
    git_tag = ""

    cmake = ""
    system = ""
    cpu = ""

    c_compiler = ""
    c_compiler_id = ""
    c_compiler_version = ""
    cxx_compiler = ""
    cxx_compiler_id = ""
    cxx_compiler_version = ""

    build_type = ""
    c_flags = ""
    cxx_flags = ""


def Log(log):
    log.info("version: " + version)
    # log.info("git hash: " + git_hash)
    log.info("git hash7: " + git_hash7)
    # log.info("git refspec: " + git_refspec)
    log.info("git describe: " + git_describe)
    # log.info("git tag: " + git_tag)
    log.info("cmake version: " + cmake)
    log.info("Platform: " + system)
    log.info("Architecture: " + cpu)

    log.info("C compiler: " + c_compiler)
    log.info("C compiler id: " + c_compiler_id)
    log.info("C compiler version: " + c_compiler_version)
    log.info("C flags: " + c_flags)
    log.info("C++ compiler: " + cxx_compiler)
    log.info("C++ compiler id: " + cxx_compiler_id)
    log.info("C++ compiler version: " + cxx_compiler_version)
    log.info("C++ flags: " + cxx_flags)

    log.info("Build type: " + build_type)

    import sys
    log.info("Python version: " + sys.version.replace("\n", "\t"))


    def package_version(package):
        import importlib
        try:
            return importlib.import_module(package).__version__
        except:
            return "<failed to import, try reinstalling>"


    def check_package(package):
        import importlib
        try:
            importlib.import_module(package)
            return True
        except:
            return False


    log.info("Packages versions:")
    for package in ["matplotlib", "seaborn", "numpy", "scipy", "Bio", "pandas"]:
        log.info("\t%s: %s" % (package, package_version(package)))
