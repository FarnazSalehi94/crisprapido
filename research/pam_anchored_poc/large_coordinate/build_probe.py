#!/usr/bin/env python3
"""Compile a bounded exact-index probe using the prototype's FM and emit code."""
import argparse,subprocess,sys
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--scratch',required=True);p.add_argument('--columba-dir',required=True);p.add_argument('--columba-build',required=True);p.add_argument('--cxx',default='g++');p.add_argument('--sa-cap',type=int,default=64);a=p.parse_args();assert 1<=a.sa_cap<=256;root=Path(__file__).resolve().parents[3];out=Path(a.scratch);out.mkdir(exist_ok=True,parents=True);proto=(root/'research/pam_anchored_poc/poc.cpp').read_text();source=('#define PROBE_SA_CAP %d\n#include <thread>\n'%a.sa_cap)+proto[:proto.index('int main(')]+r'''
static vector<string> split(const string& s) {vector<string> r;size_t b=0;while(true){size_t p=s.find('\t',b);if(p==string::npos){r.push_back(s.substr(b));return r;}r.push_back(s.substr(b,p-b));b=p+1;}}
static string revcomp(string s) {reverse(s.begin(),s.end());for(char& c:s)c=comp(c);return s;}
int main(int argc,char** argv) {try {
 if(argc!=4)throw runtime_error("usage: probe INDEX REF CASES.tsv");
 cerr<<"PHASE\tverification"<<endl;require_index_binding(argv[2],argv[1]);
 cerr<<"PHASE\tindex"<<endl;auto t=chrono::steady_clock::now();FMIndex idx(argv[1],0,true,4,false,0);cerr<<"EVENT\tindex_seconds\t"<<chrono::duration<double>(chrono::steady_clock::now()-t).count()<<endl;
 cerr<<"PHASE\treference"<<endl;t=chrono::steady_clock::now();auto ref=fasta(argv[2]);cerr<<"EVENT\treference_seconds\t"<<chrono::duration<double>(chrono::steady_clock::now()-t).count()<<endl;
 cerr<<"PHASE\tlookup"<<endl;t=chrono::steady_clock::now();ifstream input(argv[3]);if(!input)throw runtime_error("cannot read cases");string line;int total=0;
 while(getline(input,line)){if(line.empty()||line[0]=='#')continue;auto f=split(line);if(f.size()!=11)throw runtime_error("invalid case row");
  string label=f[0],kind=f[1],name=f[3],word=f[7],trueword=f[8],guide=f[9];char strand=f[2][0];size_t global=stoull(f[4]),start=stoull(f[5]),anchor=stoull(f[6]);int expected_mm=stoi(f[10]);
  size_t ci=0;while(ci<ref.size()&&ref[ci].name!=name)++ci;if(ci==ref.size()||ref[ci].offset+start!=global)throw runtime_error("contig/global offset mismatch "+label);
  string observed;if(start+23<=ref[ci].seq.size())observed=ref[ci].seq.substr(start,23);else{if(ci+1>=ref.size())throw runtime_error("missing next contig");observed=ref[ci].seq.substr(start)+ref[ci+1].seq.substr(0,23-(ref[ci].seq.size()-start));}
  if(observed!=trueword)throw runtime_error("loaded FASTA differs from independent sample "+label);
  if(word.size()!=23||guide.size()!=20||(strand!='+'&&strand!='-'))throw runtime_error("invalid probe sequence");
  string oriented=strand=='+'?word:revcomp(word);if(oriented.substr(0,20)!=guide)throw runtime_error("guide orientation mismatch");
  idx.setDirection(strand=='+'?BACKWARD:FORWARD,false);SARangePair r=idx.getCompleteRange();Counters counts;bool ok=true;
  if(strand=='+'){for(auto it=word.rbegin();it!=word.rend();++it)if(!idx.addChar(*it,r,counts)){ok=false;break;}}
  else{for(char c:word)if(!idx.addChar(c,r,counts)){ok=false;break;}}
  if(!ok)throw runtime_error("exact indexed pattern absent "+label);
  auto sr=r.getRangeSA();if(sr.width()>PROBE_SA_CAP){cout<<"SKIP\t"<<label<<"\t"<<kind<<"\t"<<strand<<"\t"<<global<<"\t"<<sr.width()<<"\n";continue;}bool located=false;
  for(length_t j=sr.getBegin();j<sr.getEnd();++j)if(size_t(idx.findSA(j))==global)located=true;
  if(!located)throw runtime_error("findSA did not locate expected position "+label);
  auto evaluate=[&](string control)->pair<bool,int>{Search s{idx,ref,Counters(),"probe",guide,"GG",strand,0,5,2,0.75,{}};s.emit(r,20,0,0,0,string(20,'='));auto it=s.hits.find(make_tuple(name,strand,anchor));return {it!=s.hits.end(),it==s.hits.end()?-1:it->second.mm};};
  auto actual=evaluate("");bool crossing=kind=="cross_boundary";if(crossing?actual.first:(!actual.first||actual.second!=expected_mm))throw runtime_error("prototype emit failed "+label);
  cout<<"CASE\t"<<label<<"\t"<<kind<<"\t"<<strand<<"\t"<<global<<"\t"<<name<<"\t"<<start<<"\t"<<anchor<<"\t"<<sr.width()<<"\t"<<(actual.first?1:0)<<"\t"<<actual.second<<"\n";
  if(label=="before_2_31"){
   for(string mode:{string("target_n"),string("pam_wildcard_n"),string("pam_required_n")}){
    size_t pos=mode=="target_n"?(strand=='+'?start+10:start+3+10):mode=="pam_wildcard_n"?(strand=='+'?start+20:start+2):(strand=='+'?start+21:start+1);
    char saved=ref[ci].seq.at(pos);if(saved=='N')throw runtime_error("N control requires resolved baseline");ref[ci].seq[pos]='N';auto x=evaluate(mode);ref[ci].seq[pos]=saved;
    bool expect=mode=="target_n";if(x.first!=expect||(expect&&x.second!=1))throw runtime_error("in-memory N control failed "+mode);
    cout<<"CONTROL\t"<<mode<<"\t"<<strand<<"\t"<<global<<"\t"<<(x.first?1:0)<<"\t"<<x.second<<"\n";
   }
  }
  ++total;
 }
 if(total==0)throw runtime_error("empty probe input");cerr<<"EVENT\tlookup_seconds\t"<<chrono::duration<double>(chrono::steady_clock::now()-t).count()<<endl;std::this_thread::sleep_for(std::chrono::milliseconds(300));cerr<<"PHASE\tpostverification"<<endl;require_index_binding(argv[2],argv[1]);cerr<<"PHASE\tdone"<<endl;
 }catch(const exception& e){cerr<<"ERROR\t"<<e.what()<<endl;return 1;} }
''';src=out/'exact_probe.cpp';src.write_text(source);c=Path(a.columba_dir);b=Path(a.columba_build);helper=root/'research/pam_anchored_poc/compatibility/index_manifest.py'
cmd=[a.cxx,'-std=gnu++14','-O3','-mpopcnt','-DALPHABET=5','-DTHIRTY_TWO','-DHAVE_ZLIB',f'-DINDEX_BINDING_HELPER="{helper}"',f'-DINDEX_BINDING_PYTHON="{Path(sys.executable).resolve()}"','-include','limits','-include','cstdint','-I'+str(c/'src'),'-I'+str(b/'_deps/parallel-hashmap-src'),'-I'+str(b/'_deps/fmtlib-src/include'),str(src),str(b/'libAuxiliaryClasses.a'),str(b/'_deps/fmtlib-build/libfmt.a'),'-lz','-pthread','-o',str(out/'exact_probe')];subprocess.run(cmd,check=True);print(out/'exact_probe')
