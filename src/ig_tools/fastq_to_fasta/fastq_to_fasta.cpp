#include <fstream>
#include <verify.hpp>

int main(int argc, char * argv[]){

    if(argc != 3) {
        std::cout << "Invalid input parameters" << std::endl <<
            "\targv[1] - input FASTQ file" << std::endl <<
            "\targv[2] - output FASTA file" << std::endl;
        return 1;
    }

	char * src_name = argv[1];
    std::ifstream src(src_name);
    VERIFY(!src.fail());

    std::ofstream dst(std::string(argv[2]).c_str());
    size_t num_reads = 0;
    while(!src.eof()){
        std::string tmp;
        getline(src, tmp);
        if(tmp == "")
            break;
		dst << ">" << tmp.substr(1, tmp.size() - 1) << std::endl;
		getline(src, tmp);
		dst << tmp << std::endl;
		getline(src, tmp);
		getline(src, tmp);
        num_reads++;
    }
	
	src.close();
	dst.close();

    std::cout << num_reads << " reads were rewritten from " << argv[1] << " to " << argv[2] << std::endl;
}
