// Independent streaming reconstruction of VANILLA Columba build-tag 21 text.
// Index layout: little-endian uint32 length, concatenated uppercase/transformed
// FASTA sequences, '$'. Source: Columba edcc7ce buildindex.cpp, seed length 0.
#include <zlib.h>
#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cstring>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
int main(int argc, char** argv) {
    try {
        if(argc!=3) throw std::runtime_error("usage: compare_index_text FASTA.gz PREFIX.txt.bin");
        gzFile input=gzopen(argv[1],"rb"); if(!input) throw std::runtime_error("cannot open FASTA");
        std::ifstream index(argv[2],std::ios::binary); if(!index) throw std::runtime_error("cannot open indexed text");
        unsigned char hdr[4]; index.read((char*)hdr,4); if(!index) throw std::runtime_error("truncated length");
        uint64_t length=uint32_t(hdr[0])|(uint32_t(hdr[1])<<8)|(uint32_t(hdr[2])<<16)|(uint32_t(hdr[3])<<24);
        std::minstd_rand gen(42); std::uniform_int_distribution<size_t> distribution(0,3);
        const std::string bases="ACGT"; const size_t capacity=1<<20;
        std::vector<char> buffer(capacity), expected; expected.reserve(capacity);
        uint64_t compared=0, ambiguous=0, records=0; bool line_start=true, header=false;
        auto compare=[&]() { if(expected.empty())return;std::vector<char> actual(expected.size());errno=0;index.read(actual.data(),actual.size());const int read_errno=errno;if(!index)throw std::runtime_error("indexed text read failed at "+std::to_string(compared)+" requested "+std::to_string(expected.size())+" received "+std::to_string(index.gcount())+" length "+std::to_string(length)+" contigs "+std::to_string(records)+" errno "+std::to_string(read_errno)+" "+std::strerror(read_errno));
            auto mismatch=std::mismatch(expected.begin(),expected.end(),actual.begin());
            if(mismatch.first!=expected.end())throw std::runtime_error("sequence mismatch at concatenated position "+std::to_string(compared+(mismatch.first-expected.begin()))+" expected "+*mismatch.first+" indexed "+*mismatch.second);
            compared+=expected.size();expected.clear(); };
        int n;
        while((n=gzread(input,buffer.data(),buffer.size()))>0) {
            for(int i=0;i<n;++i) { unsigned char c=buffer[i];
                if(c=='\n'||c=='\r'){line_start=true;header=false;continue;}
                if(line_start){header=(c=='>');line_start=false;if(header)++records;}
                if(header)continue;
                c=std::toupper(c);if(c!='A'&&c!='C'&&c!='G'&&c!='T'){c=bases[distribution(gen)];++ambiguous;}
                expected.push_back(c);if(expected.size()==capacity)compare();
            }
        }
        if(n<0)throw std::runtime_error("gzip decoding failed"); if(gzclose(input)!=Z_OK)throw std::runtime_error("gzip close failed");
        compare();char sentinel=0;index.get(sentinel);if(sentinel!='$'||length!=compared+1||index.peek()!=EOF)throw std::runtime_error("invalid indexed length/sentinel/trailing bytes");
        std::cout<<"{\"status\":\"PASS\",\"bases\":"<<compared<<",\"contigs\":"<<records<<",\"ambiguous_replacements\":"<<ambiguous<<",\"normalization\":\"uppercase; minstd_rand(42), uniform ACGT, seed length 0\"}\n";
        return 0;
    }catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}
}
