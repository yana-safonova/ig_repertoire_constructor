#include <assert.h>
#include <ostream>
#include <iostream>
#include <fstream>
#include <string.h>

int main(int argc, char * argv[]){
    if(argc != 5 || strncmp(argv[1], "-i", 3) != 0 || strncmp(argv[3], "-o", 3) != 0) {
        std::cout << "Extracts UMIs from fastq file into a separate one." << std::endl;
        std::cout << "Usage: -i <input file> -o <output file>";
        return 1;
    }
    std::cout << "Extracting barcodes." << std::endl;
    char *src_name = argv[2];
    std::ifstream src(src_name);
    assert(!src.fail());
    char *dst_name = argv[4];
    std::ofstream dst(dst_name);
    std::size_t num_reads = 0;
    while(!src.eof()){
        std::string s;
        std::getline(src, s);
        if(s == "")
            break;
        std::size_t space = s.find(' ');
        assert(space != std::string::npos);
        space = s.find(' ', space + 1);
        assert(space != std::string::npos);
        std::string id = s.substr(0, space);
        std::string umi_info = s.substr(space + 1);
        assert(!umi_info.empty());
        dst << id << std::endl;
        std::size_t colon = umi_info.find(':');
        assert(colon != std::string::npos);
        umi_info = umi_info.substr(colon + 1);
        colon = umi_info.find(':');
        assert(colon != std::string::npos);
        assert(colon * 2 + 1 == umi_info.length());
        dst << umi_info.substr(0, colon) << std::endl;
        dst << "+" << std::endl;
        dst << umi_info.substr(colon + 1) << std::endl;
        num_reads++;

        std::getline(src, s);
        std::getline(src, s);
        std::getline(src, s);
    }
    src.close();
    dst.close();
    std::cout << num_reads << " barcodes were extracted from " << src_name << " to " << dst_name << std::endl;
}
