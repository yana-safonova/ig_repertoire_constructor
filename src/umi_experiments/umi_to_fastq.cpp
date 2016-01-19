#include <assert.h>
#include <ostream>
#include <iostream>
#include <fstream>

int main(int argc, char * argv[]){
    if(argc != 3) {
        std::cout << "Invalid input parameters" << std::endl <<
        "\targv[1] - input FASTQ file with barcoded reads" << std::endl <<
        "\targv[2] - output FASTQ file made from barcodes" << std::endl;
        return 1;
    }
    char *src_name = argv[1];
    std::ifstream src(src_name);
    assert(!src.fail());
    char *dst_name = argv[2];
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
