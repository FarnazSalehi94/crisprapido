// Reconstruct only a VANILLA tag-21 .txt.bin in a NEW O_EXCL output.
// Original BWT/SA/rank components are not computed or modified here.
#include <zlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
void write_all(int fd, const char* data, size_t size) {
    while(size) { ssize_t n=::write(fd,data,size);if(n<=0)throw std::runtime_error("output write failed");data+=n;size-=n; }
}
int main(int argc, char** argv) {
    try {
        if(argc!=4)throw std::runtime_error("usage: regenerate_text FASTA.gz NEW.txt.bin EXPECTED_BASES");
        uint64_t bases=std::stoull(argv[3]);if(bases+1>UINT32_MAX)throw std::runtime_error("text exceeds 32-bit index format");
        gzFile input=gzopen(argv[1],"rb");if(!input)throw std::runtime_error("FASTA open failed");
        int fd=::open(argv[2],O_WRONLY|O_CREAT|O_EXCL,0644);if(fd<0)throw std::runtime_error("output exists or cannot be created");
        uint32_t length=bases+1;char hdr[4];for(int i=0;i<4;++i)hdr[i]=(length>>(8*i))&255;write_all(fd,hdr,4);
        std::minstd_rand gen(42);std::uniform_int_distribution<size_t> distribution(0,3);const std::string alphabet="ACGT";
        const size_t capacity=1<<20;std::vector<char> buffer(capacity),output;output.reserve(capacity);
        uint64_t count=0,ambiguous=0,records=0;bool line_start=true,header=false;int n;
        while((n=gzread(input,buffer.data(),buffer.size()))>0){for(int i=0;i<n;++i){unsigned char c=buffer[i];
            if(c=='\n'||c=='\r'){line_start=true;header=false;continue;}
            if(line_start){header=c=='>';line_start=false;if(header)++records;}if(header)continue;
            c=std::toupper(c);if(c!='A'&&c!='C'&&c!='G'&&c!='T'){c=alphabet[distribution(gen)];++ambiguous;}
            output.push_back(c);++count;if(output.size()==capacity){write_all(fd,output.data(),output.size());output.clear();}
        }}
        if(n<0||gzclose(input)!=Z_OK)throw std::runtime_error("gzip integrity failure");
        if(count!=bases)throw std::runtime_error("FASTA base count differs from frozen dictionary");
        write_all(fd,output.data(),output.size());write_all(fd,"$",1);
        if(::fsync(fd)!=0||::close(fd)!=0)throw std::runtime_error("output flush failed");
        std::cout<<"{\"status\":\"REGENERATED_NOT_YET_INDEPENDENTLY_VERIFIED\",\"bases\":"<<count<<",\"contigs\":"<<records<<",\"ambiguous_replacements\":"<<ambiguous<<"}\n";
    }catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}
}
