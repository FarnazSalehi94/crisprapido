// Research prototype; links to Columba 2.0.3 (AGPL-3.0-or-later).
#include "fmindex/fmindex.h"
#include <algorithm>
#include <cctype>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
using namespace std;
struct Seq { string name, seq; size_t offset; };
vector<Seq> fasta(const string& path) {
    ifstream f(path); if (!f) throw runtime_error("cannot open FASTA " + path);
    vector<Seq> out; string line; size_t offset=0;
    while (getline(f,line)) {
        if (!line.empty() && line[0]=='>') {
            if (!out.empty()) offset += out.back().seq.size();
            out.push_back({line.substr(1,line.find_first_of(" \t")-1),"",offset});
        } else if (!line.empty()) {
            if(out.empty()) throw runtime_error("invalid FASTA");
            for(char c:line) if(!isspace((unsigned char)c)) out.back().seq+=toupper((unsigned char)c);
        }
    } return out;
}
char comp(char c) { switch(c){case 'A':return 'T';case 'C':return 'G';case 'G':return 'C';case 'T':return 'A';case 'N':return 'N';default:throw runtime_error("guide must be ACGT");} }
struct Hit { string guide,contig,cigar,pam; char strand; size_t start,end,anchor; int mm,ins,del; };
struct Search {
 FMIndex& idx; vector<Seq>& ref; Counters counts; string name, guide; char strand; long long states=0; int maxmm, maxgap; double fraction;
 map<tuple<string,char,size_t>,Hit> hits;
 void emit(const SARangePair& range, int span, int /* indexed_mm */, int ins, int del, string ops) {
   reverse(ops.begin(),ops.end());
   auto r=range.getRangeSA(); for(length_t j=r.getBegin();j<r.getEnd();j++) {
    size_t p=idx.findSA(j); for(const auto& s:ref) {
      if(p<s.offset || p-s.offset+span+3>s.seq.size()) continue;
      size_t q=p-s.offset, start=strand=='+'?q:q+3, end=strand=='+'?q+span:q+3+span, anchor=strand=='+'?end:start;
      string pam=strand=='+'?s.seq.substr(end,3):s.seq.substr(q,3);
      if(pam.find_first_not_of("ACGT")!=string::npos) continue;
      if(strand=='-') { reverse(pam.begin(),pam.end()); for(char& c:pam)c=comp(c); }
      if(pam.size()!=3 || pam[1]!='G'||pam[2]!='G'||string("ACGT").find(pam[0])==string::npos) continue;
      // The indexed letter at a reference N is a replacement, so recompute
      // every aligned column and both eligibility counters from the FASTA.
      string target=s.seq.substr(start,end-start);
      if(target.find_first_not_of("ACGTN")!=string::npos) continue;
      if(strand=='-') {reverse(target.begin(),target.end());for(char& c:target)c=comp(c);}
      size_t gi=0,ti=0; int true_mm=0,true_matches=0; bool valid=true; string true_ops;
      for(char op:ops) {
        if(op=='I') {if(gi>=guide.size()){valid=false;break;} ++gi;true_ops+='I';}
        else if(op=='D') {if(ti>=target.size()){valid=false;break;} ++ti;true_ops+='D';}
        else {if(gi>=guide.size()||ti>=target.size()){valid=false;break;}
          bool same=guide[gi]==target[ti];true_ops+=same?'=':'X';true_matches+=same;true_mm+=!same;++gi;++ti;}
      }
      if(!valid || gi!=guide.size() || ti!=target.size() || true_mm>maxmm || true_matches+1e-9<fraction*guide.size()) continue;
      string cigar; for(size_t i=0;i<true_ops.size();) {size_t j=i+1;while(j<true_ops.size()&&true_ops[j]==true_ops[i])++j;cigar+=to_string(j-i)+true_ops[i];i=j;}
      Hit h{name,s.name,cigar,pam,strand,start,end,anchor,true_mm,ins,del}; auto key=make_tuple(s.name,strand,anchor);
      auto it=hits.find(key); if(it==hits.end() || make_tuple(true_mm+ins+del,true_mm,cigar,start)<make_tuple(it->second.mm+it->second.ins+it->second.del,it->second.mm,it->second.cigar,it->second.start)) hits[key]=h;
    }
   }
 }
 void walk(SARangePair r,int i,int mm,int gaptype,int gaplen,int span,int matches,string ops) {
   // On resolved guide positions a reference N is always a true mismatch.
   // Its A/C/G/T index replacement can only lower that mismatch count and
   // increase matches, so these pruning tests remain conservative for N.
   states++; if(mm>maxmm || matches+i+1 < fraction*guide.size()-1e-9) return;
   if(i<0) emit(r,span,mm,(gaptype==1 || gaptype==3)?gaplen:0,(gaptype==2 || gaptype==4)?gaplen:0,ops);
   // Gap in guide (D): consumes reference only. A single contiguous event.
   if((gaptype==0 || gaptype==2) && gaplen<maxgap) for(char c:string("ACGT")) { auto n=r; if(idx.addChar(c,n,counts)) walk(n,i,mm,2,gaplen+1,span+1,matches,ops+'D'); }
   if(i<0) return;
   // Gap in reference (I): consumes guide only.
   if((gaptype==0 || gaptype==1) && gaplen<maxgap) walk(r,i-1,mm,1,gaplen+1,span,matches,ops+'I');
   char expected=strand=='+'?guide[i]:comp(guide[i]);
   for(char c:string("ACGT")) { auto n=r; if(idx.addChar(c,n,counts)) {int m=(c!=expected); walk(n,i-1,mm+m,gaptype==1?3:gaptype==2?4:gaptype,gaplen,span+1,matches+!m,ops+(m?'X':'='));} }
 }
 void run() {
   for(char n:string("ACGT")) { idx.setDirection(strand=='+'?BACKWARD:FORWARD,false); SARangePair r=idx.getCompleteRange(); string seed=strand=='+'?string("GG")+n:string("CC")+comp(n); bool ok=true; for(char c:seed) if(!idx.addChar(c,r,counts)){ok=false;break;} if(ok) walk(r,guide.size()-1,0,0,0,0,0,""); }
 }
};
int main(int argc,char** argv) { try {
 if(argc!=7) {cerr<<"usage: poc INDEX REF.fa GUIDES.fa MAX_MISMATCH MAX_GAP MIN_MATCH_FRACTION\n";return 2;}
 auto ref=fasta(argv[2]), guides=fasta(argv[3]); FMIndex idx(argv[1],0,true,4,false,0); long long states=0; double search_seconds=0;
 cout<<"guide\tcontig\tstrand\tstart\tend\tanchor\tpam\tcigar\tmismatches\tinserted_guide_bases\tdeleted_guide_bases\n";
 for(auto& g:guides) for(char strand:string("+-")) { for(char c:g.seq) if(string("ACGT").find(c)==string::npos) throw runtime_error("guide must be ACGT"); Search s{idx,ref,Counters(),g.name,g.seq,strand,0,stoi(argv[4]),stoi(argv[5]),stod(argv[6]),{}}; auto started=chrono::steady_clock::now(); s.run(); search_seconds+=chrono::duration<double>(chrono::steady_clock::now()-started).count(); states+=s.states; for(auto& kv:s.hits){auto& h=kv.second;cout<<h.guide<<'\t'<<h.contig<<'\t'<<h.strand<<'\t'<<h.start<<'\t'<<h.end<<'\t'<<h.anchor<<'\t'<<h.pam<<'\t'<<h.cigar<<'\t'<<h.mm<<'\t'<<h.ins<<'\t'<<h.del<<'\n';} }
 cerr<<"states="<<states<<" search_wall_seconds="<<search_seconds<<'\n';
 } catch(const exception& e){cerr<<e.what()<<'\n';return 1;} }
