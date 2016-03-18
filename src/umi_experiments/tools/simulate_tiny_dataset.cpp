#include <seqan/seq_io.h>

const size_t BASE_READS = 100;
const size_t AMPLIFIED_READS = 40;
const size_t SEED = 34956;
const size_t MAX_MUTATIONS = 10;

void add_read(std::vector<seqan::Dna5String>& reads,
              std::vector<seqan::CharString>& read_ids,
              size_t& cnt,
              const seqan::Dna5String& base_read,
              const std::string& umi) {
    seqan::Dna5String read(base_read);
    size_t mutations = rand() % MAX_MUTATIONS;
    for (size_t i = 0; i < mutations; i ++) {
        size_t position = rand() % length(read);
        size_t nt = rand() % 4;
        read[position] = nt;
    }
    reads.push_back(read);
    read_ids.emplace_back("read_id__" + std::to_string(cnt ++) + "__UMI:" + umi);
}

int main() {
    std::srand(SEED);
    const std::string umi("GGGGGGGGGGGG");
    std::vector<seqan::Dna5String> reads;
    std::vector<seqan::CharString> read_ids;
    size_t read_cnt = 0;
    seqan::Dna5String base_read("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
    for (size_t i = 0; i < BASE_READS; i ++) {
        for (size_t j = 0; j < length(base_read); j ++) {
            base_read[j] = rand() % 4;
        }
        for (size_t j = 0; j < AMPLIFIED_READS; j++) {
            add_read(reads, read_ids, read_cnt, base_read, umi);
        }
    }
    for (size_t i = 1; i < read_cnt; i ++) {
        size_t p = rand() % ( i + 1 );
        std::swap(reads[i], reads[p]);
        std::swap(read_ids[i], read_ids[p]);
    }
    seqan::SeqFileOut out_file("tiny_dataset.fasta");
    writeRecords(out_file, read_ids, reads);
    return 0;
}
