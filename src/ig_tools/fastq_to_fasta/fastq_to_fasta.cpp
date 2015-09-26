#include <utils/include_me.hpp>
using namespace std;

int main(int argc, char * argv[]){

    if(argc != 3) {
        cout << "Invalid input parameters" << endl <<
            "\targv[1] - input FASTQ file" << endl <<
            "\targv[2] - output FASTA file" << endl;
        return 1;
    }

	char * src_name = argv[1];
	ifstream src(src_name);
    assert(!src.fail());

	ofstream dst(string(argv[2]).c_str());
    size_t num_reads = 0;
    while(!src.eof()){
        string tmp;
        getline(src, tmp);
        if(tmp == "")
            break;
		dst << ">" << tmp.substr(1, tmp.size() - 1) << endl;
		getline(src, tmp);
		dst << tmp << endl;
		getline(src, tmp);
		getline(src, tmp);
        num_reads++;
    }
	
	src.close();
	dst.close();

    cout << num_reads << " reads were rewritten from " << argv[1] << " to " << argv[2] << endl;
}
