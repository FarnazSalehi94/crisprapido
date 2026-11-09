// Isolated research executable. NOT accepted by the integrated backend protocol.
// Reuse unchanged emission/eligibility and unpruned search for the control.
#define main original_unpruned_main
#include "../poc.cpp"
#undef main
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <sys/resource.h>

struct Score { bool available; double value; };
struct Tables {
 std::map<std::string,double> mm, pam;
 Score factors[21][5][5]; bool contract=true;
 static int base(char c) { auto p=std::string("ACGT-").find(c); return p==std::string::npos?-1:int(p); }
 static std::map<std::string,double> read(const std::string& path) {
  std::ifstream in(path); if(!in) throw std::runtime_error("cannot read score table");
  std::map<std::string,double> result; std::string line;
  while(std::getline(in,line)) { std::istringstream row(line); std::string key,value;
   if(!(row>>key>>value)) continue;
   size_t used=0; double v=std::stod(value,&used); if(used!=value.size()) throw std::runtime_error("invalid score value");
   result[key]=v;
  } return result;
 }
 Tables(const std::string& mismatch,const std::string& pams):mm(read(mismatch)),pam(read(pams)) {
  for(const auto& table:{mm,pam}) for(const auto& kv:table)
   if(!std::isfinite(kv.second)||kv.second<0||kv.second>1) contract=false;
  for(int p=1;p<=20;++p) for(int g=0;g<5;++g) for(int t=0;t<5;++t) {
   char gb=std::string("ACGT-")[g],tb=std::string("ACGT-")[t]; Score f{true,1};
   if(gb!=tb && !(p==1 && (gb=='-'||tb=='-'))) {
    char rna=gb=='T'?'U':gb, dna=tb=='-'?'-':comp(tb);
    auto it=mm.find(std::string("r")+rna+":d"+dna+","+std::to_string(p));
    f=it==mm.end()?Score{false,0}:Score{true,it->second};
   }
   factors[p][g][t]=f;
   if(g!=t && !(g==4&&t==4) && !f.available) contract=false;
  }
  for(char a:std::string("ACGT")) for(char b:std::string("ACGT"))
   if(!pam.count(std::string(1,a)+b)) contract=false;
 }
 Score column(char g,char t,int pos) const {
  if(g==t || (pos==1 && (g=='-'||t=='-'))) return {true,1};
  int gi=base(g),ti=base(t); if(gi<0||ti<0||pos<1||pos>20) return {false,0};
  return factors[pos][gi][ti];
 }
 Score score(const std::string& guide,const std::string& target,const std::string& ops,const std::string& suffix) const {
  size_t q=0,t=0; int dna=0; bool has_i=false; std::string ag,at;
  for(char op:ops) {
   if(op=='I') { if(q>=guide.size()) return {false,0}; ag+=guide[q++];at+='-';has_i=true; }
   else if(op=='D') {if(t>=target.size()) return {false,0};ag+='-';at+=target[t++];++dna;}
   else {if(q>=guide.size()||t>=target.size()) return {false,0};ag+=guide[q++];at+=target[t++];}
  }
  if(q!=guide.size()||t!=target.size()) return {false,0};
  if(dna&&!has_i) {ag=ag.substr(dna);at=at.substr(dna);}
  if(ag.size()!=20||at.size()!=20||!pam.count(suffix)) return {false,0};
  double value=1;
  for(size_t i=0;i<20;++i) { auto f=column(ag[i],at[i],i+1);if(!f.available)return {false,0}; value*=f.value; }
  value*=pam.at(suffix);return {true,value};
 }
};

std::string expand(const std::string& cigar) {
 std::string ops;size_t n=0;for(char c:cigar) {
  if(std::isdigit(static_cast<unsigned char>(c)))n=n*10+(c-'0');
  else {if(!n||std::string("=XID").find(c)==std::string::npos)throw std::runtime_error("invalid CIGAR");ops.append(n,c);n=0;}
 } if(n)throw std::runtime_error("unfinished CIGAR"); return ops;
}
std::string compress(const std::string& ops) {
 std::string cigar;for(size_t i=0;i<ops.size();) {size_t j=i+1;while(j<ops.size()&&ops[j]==ops[i])++j;cigar+=std::to_string(j-i)+ops[i];i=j;}return cigar;
}
bool preferred(const Hit& a,const Hit& b) {
 return std::make_tuple(a.mm+a.ins+a.del,a.mm,a.cigar,a.start)<std::make_tuple(b.mm+b.ins+b.del,b.mm,b.cigar,b.start);
}
std::string oriented_target(const Seq& seq,const Hit& h) {
 auto t=seq.seq.substr(h.start,h.end-h.start);if(h.strand=='-') {std::reverse(t.begin(),t.end());for(char& c:t)c=comp(c);}return t;
}
const Seq& sequence(const std::vector<Seq>& refs,const std::string& name) {
 for(const auto& s:refs)if(s.name==name)return s;throw std::runtime_error("unknown contig");
}

// No CFD pruning or optimization in reconstruction: enumerate every one-event
// gap word at every forward guide boundary, including both ID and DI.
Hit reconstruct(const Hit& anchor,const Seq& seq,const std::string& guide,int m,int z,double fraction,long long& trials) {
 std::vector<std::string> words{""};for(int n=1;n<=z;++n)for(int bits=0;bits<(1<<n);++bits) {
  std::string word;for(int j=0;j<n;++j)word+=(bits&(1<<j))?'I':'D';words.push_back(word);
 }
 bool found=false;Hit best;
 for(const auto& word:words) {
  int ins=std::count(word.begin(),word.end(),'I'),del=std::count(word.begin(),word.end(),'D');
  int span=int(guide.size())-ins+del;if(span<=0||ins>int(guide.size()))continue;
  size_t start,end;
  if(anchor.strand=='+') {if(anchor.anchor<size_t(span))continue;end=anchor.anchor;start=end-span;}
  else {start=anchor.anchor;if(start>seq.seq.size()||size_t(span)>seq.seq.size()-start)continue;end=start+span;}
  Hit h{anchor.guide,anchor.contig,"",anchor.pam,anchor.strand,start,end,anchor.anchor,0,ins,del,0};
  std::string target=oriented_target(seq,h);
  const size_t last_prefix=word.empty()?0:guide.size()-ins;
  for(size_t prefix=0;prefix<=last_prefix;++prefix) {
   ++trials;std::string ops=std::string(prefix,'=')+word+std::string(guide.size()-ins-prefix,'=');
   size_t gi=0,ti=0;int mm=0,matches=0;std::string actual;
   for(char op:ops) {
    if(op=='I') {++gi;actual+='I';} else if(op=='D') {++ti;actual+='D';}
    else {bool same=guide.at(gi++)==target.at(ti++);mm+=!same;matches+=same;actual+=same?'=':'X';}
   }
   if(mm>m||matches==0||!(float(matches)/float(guide.size())*100.0f>=float(fraction)*100.0f))continue;
   h.cigar=compress(actual);h.mm=mm;h.matches=matches;
   if(!found||preferred(h,best)) {found=true;best=h;}
  }
 }
 if(!found)throw std::runtime_error("discovered anchor has no reconstructible eligible alignment");
 return best;
}

struct Discovery:Search {
 const Tables& tables;bool active,return_on_bound;double tau;
 long long pruned=0,would_prune=0,bound_checks=0,bound_updates=0;
 Discovery(FMIndex& idx,std::vector<Seq>& refs,const Seq& g,const std::string& pam,char strand,int m,int z,double fraction,
           const Tables& tables,bool active,bool prune,double tau):
  Search{idx,refs,Counters(),g.name,g.seq,pam,strand,0,m,z,fraction,{}},tables(tables),active(active),return_on_bound(prune),tau(tau) {}
 bool unknown(int i,int state,int gaplen,int ins,int del) const {
  if(ins&&del)return true;
  if(state==2||maxgap<2)return false;
  if(state==0)return i>=0;
  if(gaplen>=maxgap)return false;
  return ins||(del&&i>=0);
 }
 double update(double upper,size_t columns,char g,char t) {
  if(!active||columns>=20)return upper;
  ++bound_updates;auto f=tables.column(g,t,20-columns);
  if(!f.available)throw std::runtime_error("validated bound factor missing");
  return std::nextafter(upper*f.value,std::numeric_limits<double>::infinity());
 }
 void walk_bound(SARangePair r,int i,int mm,int state,int gaplen,int ins,int del,int span,int matches,std::string ops,double upper) {
  ++states;if(mm>maxmm||!passes_fraction(matches+i+1))return;
  if(active&&!unknown(i,state,gaplen,ins,del)) {
   ++bound_checks;
   double safe=std::nextafter(upper*(1+64*std::numeric_limits<double>::epsilon())+64*std::numeric_limits<double>::denorm_min(),std::numeric_limits<double>::infinity());
   if(safe<tau) {++would_prune;if(return_on_bound){++pruned;return;}}
  }
  if(i<0)emit(r,span,mm,ins,del,ops);
  if(state!=2&&gaplen<maxgap)for(char c:std::string("ACGT")) {auto n=r;if(idx.addChar(c,n,counts))
   walk_bound(n,i,mm,1,gaplen+1,ins,del+1,span+1,matches,ops+'D',update(upper,ops.size(),'-',strand=='+'?c:comp(c)));}
  if(i<0)return;
  if(state!=2&&gaplen<maxgap)
   walk_bound(r,i-1,mm,1,gaplen+1,ins+1,del,span,matches,ops+'I',update(upper,ops.size(),guide[i],'-'));
  char expected=strand=='+'?guide[i]:comp(guide[i]);
  for(char c:std::string("ACGT")) {auto n=r;if(idx.addChar(c,n,counts)) {int mismatch=c!=expected;
   walk_bound(n,i-1,mm+mismatch,state?2:0,gaplen,ins,del,span+1,matches+!mismatch,ops+(mismatch?'X':'='),update(upper,ops.size(),guide[i],strand=='+'?c:comp(c)));}}
 }
 void discover() {
  for(char n:std::string("ACGT")) {
   idx.setDirection(strand=='+'?BACKWARD:FORWARD,false);SARangePair r=idx.getCompleteRange();
   std::string seed=strand=='+'?std::string(1,requested_pam[1])+requested_pam[0]+n:std::string(1,comp(requested_pam[1]))+comp(requested_pam[0])+comp(n);
   bool ok=true;for(char c:seed)if(!idx.addChar(c,r,counts)){ok=false;break;}
   if(ok)walk_bound(r,guide.size()-1,0,0,0,0,0,0,0,"",active?tables.pam.at(requested_pam):1);
  }
 }
};

using Clock=std::chrono::steady_clock;
double seconds(Clock::time_point t) {return std::chrono::duration<double>(Clock::now()-t).count();}
int main(int argc,char** argv) {try {
 if(argc!=13)throw std::runtime_error("RESEARCH ONLY: native INDEX REF GUIDES M Z FRACTION PAM B MODE(unpruned|overhead|pruned) MIN_CFD MM_TABLE PAM_TABLE");
 auto integer=[](const char* s){size_t n=0;int v=std::stoi(s,&n);if(s[n]||v<0)throw std::runtime_error("invalid nonnegative integer");return v;};
 auto real=[](const char* s){size_t n=0;double v=std::stod(s,&n);if(s[n]||!std::isfinite(v)||v<0||v>1)throw std::runtime_error("invalid finite fraction/threshold");return v;};
 int m=integer(argv[4]),z=integer(argv[5]),b=integer(argv[8]);double frac=real(argv[6]),tau=real(argv[10]);
 std::string pam=argv[7],mode=argv[9];
 if(z>2||b>1||pam.size()!=2||pam.find_first_not_of("ACGT")!=std::string::npos)throw std::runtime_error("unsupported research parameters");
 if(mode!="unpruned"&&mode!="pruned"&&mode!="overhead")throw std::runtime_error("invalid mode");if(!b)z=0;
 auto pre=Clock::now();require_index_binding(argv[2],argv[1]);double pre_seconds=seconds(pre);
 auto refs=fasta(argv[2]),guides=fasta(argv[3]);if(refs.empty()||guides.empty())throw std::runtime_error("empty input");
 bool has_n=false;std::set<std::string> rn,gn;
 for(const auto& r:refs){if(r.name.empty()||r.seq.empty()||r.seq.find_first_not_of("ACGTN")!=std::string::npos||!rn.insert(r.name).second)throw std::runtime_error("invalid reference");has_n|=r.seq.find('N')!=std::string::npos;}
 for(const auto& g:guides)if(g.name.empty()||g.seq.empty()||g.seq.find_first_not_of("ACGT")!=std::string::npos||!gn.insert(g.name).second)throw std::runtime_error("invalid guide");
 Tables tables(argv[11],argv[12]);auto load=Clock::now();FMIndex idx(argv[1],0,true,4,false,0);double index_seconds=seconds(load);
 long long states=0,pruned=0,checks=0,updates=0,would=0,trials=0,candidates=0,refined=0,reported=0,unavailable=0;
 double search_seconds=0,refine_seconds=0,filter_seconds=0;bool any_active=false;
 std::ostringstream output;
 output<<"guide\tcontig\tstrand\tstart\tend\tanchor\tpam\tcigar\tmismatches\tinserted_guide_bases\tdeleted_guide_bases\tbulge_events\tmatch_count\tmatch_fraction\tcfd\n";
 auto algorithm_start=Clock::now();
 for(const auto& g:guides)for(char strand:std::string("+-")) {
  std::map<std::tuple<std::string,char,size_t>,Hit> hits;bool lost_paths=false;
  bool active=!has_n&&g.seq.size()==20&&tables.contract&&tau>0;any_active|=active;
  auto started=Clock::now();
  if(mode=="unpruned") {Search s{idx,refs,Counters(),g.name,g.seq,pam,strand,0,m,z,frac,{}};s.run();states+=s.states;hits=std::move(s.hits);}
  else {Discovery s(idx,refs,g,pam,strand,m,z,frac,tables,active,mode=="pruned",tau);s.discover();states+=s.states;pruned+=s.pruned;would+=s.would_prune;checks+=s.bound_checks;updates+=s.bound_updates;lost_paths=s.pruned>0;hits=std::move(s.hits);}
  search_seconds+=seconds(started);candidates+=hits.size();
  for(auto& kv:hits) {
   auto& h=kv.second;const auto& seq=sequence(refs,h.contig);
   if(lost_paths) {auto t=Clock::now();h=reconstruct(h,seq,g.seq,m,z,frac,trials);refine_seconds+=seconds(t);++refined;}
   auto t=Clock::now();auto score=tables.score(g.seq,oriented_target(seq,h),expand(h.cigar),h.pam.substr(1));
   bool unknown=!score.available||!std::isfinite(score.value);bool keep=unknown||score.value>=tau;
   filter_seconds+=seconds(t);if(!keep)continue;
   ++reported;unavailable+=unknown;
   output<<h.guide<<'\t'<<h.contig<<'\t'<<h.strand<<'\t'<<h.start<<'\t'<<h.end<<'\t'<<h.anchor<<'\t'<<h.pam<<'\t'<<h.cigar<<'\t'<<h.mm<<'\t'<<h.ins<<'\t'<<h.del<<'\t'<<((h.ins||h.del)?1:0)<<'\t'<<h.matches<<'\t'<<(double(h.matches)/g.seq.size())<<'\t';
   if(unknown)output<<"NA";else output<<std::setprecision(17)<<score.value;output<<'\n';
  }
 }
 double algorithm_seconds=seconds(algorithm_start);
 auto post=Clock::now();require_index_binding(argv[2],argv[1]);double post_seconds=seconds(post);
 std::cout<<output.str();std::cout.flush();
 struct rusage usage;getrusage(RUSAGE_SELF,&usage);
 std::cerr<<std::setprecision(12)<<"CFD_METRICS {\"states\":"<<states<<",\"pruned\":"<<pruned<<",\"would_prune\":"<<would<<",\"bound_checks\":"<<checks<<",\"bound_updates\":"<<updates<<",\"reconstruction_trials\":"<<trials<<",\"discovery_anchors\":"<<candidates<<",\"refined_anchors\":"<<refined<<",\"reported\":"<<reported<<",\"unavailable\":"<<unavailable<<",\"search_seconds\":"<<search_seconds<<",\"refinement_seconds\":"<<refine_seconds<<",\"filter_seconds\":"<<filter_seconds<<",\"algorithm_seconds\":"<<algorithm_seconds<<",\"native_peak_rss_kib\":"<<usage.ru_maxrss<<",\"index_load_seconds\":"<<index_seconds<<",\"precheck_seconds\":"<<pre_seconds<<",\"postcheck_seconds\":"<<post_seconds<<",\"bound_eligible\":"<<(any_active?"true":"false")<<",\"has_reference_n\":"<<(has_n?"true":"false")<<",\"matrix_contract\":"<<(tables.contract?"true":"false")<<"}\n";
 }catch(const std::exception& e){std::cerr<<"research CFD prototype: "<<e.what()<<'\n';return 1;}}
