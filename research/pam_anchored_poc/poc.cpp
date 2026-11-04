// Research prototype; links to Columba 2.0.3 (AGPL-3.0-or-later).
#include "fmindex/fmindex.h"
#include <algorithm>
#include <cctype>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
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
struct Hit { string guide,contig,cigar,pam; char strand; size_t start,end,anchor; int mm,ins,del,matches; };
struct Search {
 FMIndex& idx; vector<Seq>& ref; Counters counts; string name, guide, requested_pam; char strand; long long states=0; int maxmm, maxgap; double fraction;
 map<tuple<string,char,size_t>,Hit> hits;
 bool passes_fraction(int matches) const { return (float(matches)/float(guide.size())*100.0f) >= (float(fraction)*100.0f); }
 void emit(const SARangePair& range, int span, int /* indexed_mm */, int ins, int del, string ops) {
   reverse(ops.begin(),ops.end());
   auto r=range.getRangeSA(); for(length_t j=r.getBegin();j<r.getEnd();j++) {
    size_t p=idx.findSA(j); for(const auto& s:ref) {
      if(p<s.offset || p-s.offset+span+3>s.seq.size()) continue;
      size_t q=p-s.offset, start=strand=='+'?q:q+3, end=strand=='+'?q+span:q+3+span, anchor=strand=='+'?end:start;
      string pam=strand=='+'?s.seq.substr(end,3):s.seq.substr(q,3);
      if(pam.find_first_not_of("ACGT")!=string::npos) continue;
      if(strand=='-') { reverse(pam.begin(),pam.end()); for(char& c:pam)c=comp(c); }
      if(pam.size()!=3 || pam[1]!=requested_pam[0]||pam[2]!=requested_pam[1]||string("ACGT").find(pam[0])==string::npos) continue;
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
      if(!valid || gi!=guide.size() || ti!=target.size() || true_mm>maxmm || true_matches==0 || !passes_fraction(true_matches)) continue;
      string cigar; for(size_t i=0;i<true_ops.size();) {size_t j=i+1;while(j<true_ops.size()&&true_ops[j]==true_ops[i])++j;cigar+=to_string(j-i)+true_ops[i];i=j;}
      Hit h{name,s.name,cigar,pam,strand,start,end,anchor,true_mm,ins,del,true_matches}; auto key=make_tuple(s.name,strand,anchor);
      auto it=hits.find(key); if(it==hits.end() || make_tuple(true_mm+ins+del,true_mm,cigar,start)<make_tuple(it->second.mm+it->second.ins+it->second.del,it->second.mm,it->second.cigar,it->second.start)) hits[key]=h;
    }
   }
 }
 // gap_state: 0=no event, 1=currently in one contiguous I/D event,
 // 2=event closed by an aligned column. Production counts adjacent I/D
 // columns as one event and applies -z to their combined length.
 void walk(SARangePair r,int i,int mm,int gap_state,int gaplen,int ins,int del,int span,int matches,string ops) {
   // Substituted reference N can only lower indexed mismatch cost.
   states++; if(mm>maxmm || !passes_fraction(matches+i+1)) return;
   if(i<0) emit(r,span,mm,ins,del,ops);
   if(gap_state!=2 && gaplen<maxgap) for(char c:string("ACGT")) {auto n=r;if(idx.addChar(c,n,counts))walk(n,i,mm,1,gaplen+1,ins,del+1,span+1,matches,ops+'D');}
   if(i<0) return;
   if(gap_state!=2 && gaplen<maxgap)walk(r,i-1,mm,1,gaplen+1,ins+1,del,span,matches,ops+'I');
   char expected=strand=='+'?guide[i]:comp(guide[i]);
   for(char c:string("ACGT")) {auto n=r;if(idx.addChar(c,n,counts)){int m=(c!=expected);walk(n,i-1,mm+m,gap_state?2:0,gaplen,ins,del,span+1,matches+!m,ops+(m?'X':'='));}}
 }
 void run() {
   for(char n:string("ACGT")) { idx.setDirection(strand=='+'?BACKWARD:FORWARD,false); SARangePair r=idx.getCompleteRange(); string seed=strand=='+'?string(1,requested_pam[1])+requested_pam[0]+n:string(1,comp(requested_pam[1]))+comp(requested_pam[0])+comp(n); bool ok=true; for(char c:seed) if(!idx.addChar(c,r,counts)){ok=false;break;} if(ok) walk(r,guide.size()-1,0,0,0,0,0,0,0,""); }
 }
};
int main(int argc,char** argv) { try {
 if(argc!=7 && argc!=9) {cerr<<"usage: poc INDEX REF.fa GUIDES.fa MAX_MISMATCH MAX_GAP MIN_MATCH_FRACTION [PAM MAX_BULGE_EVENTS]\n";return 2;}
 auto number=[](const char* x){size_t p=0;int n=stoi(x,&p);if(x[p] || n<0)throw runtime_error("limits must be non-negative integers");return n;};
 int maxmm=number(argv[4]), maxgap=number(argv[5]), maxbulges=argc==9?number(argv[8]):1;
 size_t parsed=0; double fraction=stod(argv[6],&parsed);if(argv[6][parsed] || !(fraction>=0.0 && fraction<=1.0))throw runtime_error("match fraction must be between 0 and 1");
 string pam=argc==9?argv[7]:"GG";for(char& c:pam)c=toupper((unsigned char)c);
 if(pam.size()!=2 || pam.find_first_not_of("ACGT")!=string::npos)throw runtime_error("PAM must be exactly two A/C/G/T bases");
 if(maxbulges>1)throw runtime_error("unsupported maximum bulge events: prototype supports only 0 or 1");
 if(maxgap>2)throw runtime_error("unsupported maximum bulge size: prototype supports 0, 1 or 2");
 if(maxbulges==0)maxgap=0;
 auto ref=fasta(argv[2]), guides=fasta(argv[3]);
 if(ref.empty() || guides.empty())throw runtime_error("reference and guide FASTA must be nonempty");
 set<string> ref_names, guide_names;
 for(const auto& r:ref) {if(r.name.empty() || r.seq.empty() || r.seq.find_first_not_of("ACGTN")!=string::npos)throw runtime_error("reference must be nonempty A/C/G/T/N");if(!ref_names.insert(r.name).second)throw runtime_error("duplicate reference name");}
 for(const auto& g:guides) {if(g.name.empty() || g.seq.empty() || g.seq.find_first_not_of("ACGT")!=string::npos)throw runtime_error("guide must be nonempty ACGT");if(!guide_names.insert(g.name).second)throw runtime_error("duplicate guide name");}
 FMIndex idx(argv[1],0,true,4,false,0); long long states=0; double search_seconds=0;
 cout<<"guide\tcontig\tstrand\tstart\tend\tanchor\tpam\tcigar\tmismatches\tinserted_guide_bases\tdeleted_guide_bases\tbulge_events\tmatch_count\tmatch_fraction\n";
 for(auto& g:guides) for(char strand:string("+-")) { Search s{idx,ref,Counters(),g.name,g.seq,pam,strand,0,maxmm,maxgap,fraction,{}}; auto started=chrono::steady_clock::now(); s.run(); search_seconds+=chrono::duration<double>(chrono::steady_clock::now()-started).count(); states+=s.states; for(auto& kv:s.hits){auto& h=kv.second;cout<<h.guide<<'\t'<<h.contig<<'\t'<<h.strand<<'\t'<<h.start<<'\t'<<h.end<<'\t'<<h.anchor<<'\t'<<h.pam<<'\t'<<h.cigar<<'\t'<<h.mm<<'\t'<<h.ins<<'\t'<<h.del<<'\t'<<((h.ins||h.del)?1:0)<<'\t'<<h.matches<<'\t'<<(double(h.matches)/g.seq.size())<<'\n';} }
 cerr<<"states="<<states<<" search_wall_seconds="<<search_seconds<<'\n';
 } catch(const exception& e){cerr<<e.what()<<'\n';return 1;} }
