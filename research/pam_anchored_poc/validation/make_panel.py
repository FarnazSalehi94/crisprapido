#!/usr/bin/env python3
"""Build nonoverlapping real-reference extracts and direct PAM-first exhaustive oracle."""
import argparse,collections,hashlib,json
from pathlib import Path
DNA='ACGT'; rc=lambda s:s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def read_fa(p):
    d={}; key=None
    for line in Path(p).read_text().splitlines():
        if line.startswith('>'): key=line[1:].split()[0];d[key]=''
        else:d[key]+=line.strip().upper()
    return d

def main():
    ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('out');args=ap.parse_args();src=Path(args.source);out=Path(args.out);out.mkdir(parents=True,exist_ok=False)
    ws=json.loads((src/'windows.json').read_text()); seqs={}; info=[]; guides={}; design=[]
    for wi,w in enumerate(ws):
        original=next(iter(read_fa(src/w['name']/'ref.fa').values())); offset=({1:13000,3:1865000}.get(wi,400000+wi*10000)); seq=original[offset:offset+1000]; assert len(seq)==1000 and set(seq)<=set(DNA)
        name=w['name'];seqs[name]=seq;info.append({'name':name,'parent':w['region'],'parent_offset':offset,'size':len(seq),'sha256':hashlib.sha256(seq.encode()).hexdigest(),'artificial_N':False})
        for strand in '+-':
            ori=seq if strand=='+' else rc(seq); pams=[i for i in range(60,940) if ori[i] in DNA and ori[i+1:i+3]=='GG']; assert pams
            anchor=pams[len(pams)//2]
            for case in ['exact','mm5','I1','I2','D1','D2','I2mm3','D2mm5']:
                d=0 if case in ('exact','mm5') else (-int(case[1]) if case[0]=='I' else int(case[1]))
                target=ori[anchor-20-d:anchor]
                if d<0: guide=target[:8]+'AC'[:(-d)]+target[8:]
                elif d>0: guide=target[:8]+target[8+d:]
                else: guide=target
                if 'mm' in case:
                    chars=list(guide)
                    for k in [1,4,7,11,15][:3 if case.endswith('mm3') else 5]:chars[k]=DNA[(DNA.index(chars[k])+1)%4]
                    guide=''.join(chars)
                assert len(guide)==20
                gid=name+'_'+strand.replace('+','plus').replace('-','minus')+'_'+case;guides[gid]=guide
                design.append({'guide':gid,'contig':name,'strand':strand,'anchor':anchor if strand=='+' else len(seq)-anchor-3,'case':case})
        # An explicitly modified extract is a stress fixture, not an unmodified genome.
        edited=list(seq); chosen=[d for d in design if d['contig']==name and d['case']=='exact']; positions=[40,41,42,43,44,950,951,952]
        for d in chosen: positions.append(d['anchor']-5 if d['strand']=='+' else d['anchor']+3+4)
        for i in positions: edited[i]='N'
        seqs[name+'_N']=''.join(edited)
        info.append({'name':name+'_N','parent':w['region'],'parent_offset':offset,'size':len(seq),'artificial_N':True,'N_positions':positions})
    (out/'ref.fa').write_text(''.join('>'+n+'\n'+s+'\n' for n,s in seqs.items()));(out/'guides.fa').write_text(''.join('>'+n+'\n'+s+'\n' for n,s in guides.items()));(out/'panel.json').write_text(json.dumps({'windows':info,'designs':design},indent=2)+'\n')
    truth={}
    for contig,seq in seqs.items():
        for strand in '+-':
            ori=seq if strand=='+' else rc(seq)
            for end in range(18,len(ori)-2):
                pam=ori[end:end+3]
                if len(pam)!=3 or pam[0] not in DNA or pam[1:]!='GG': continue
                anchor=end if strand=='+' else len(seq)-end-3
                for delta in (-2,-1,0,1,2):
                    length=20+delta
                    if end<length:continue
                    t=ori[end-length:end]
                    for gid,g in guides.items():
                        if delta==0: variants=[(g,t)]
                        elif delta<0: variants=((g[:i]+g[i-delta:],t) for i in range(length+1))
                        else: variants=((g,t[:i]+t[i+delta:]) for i in range(21))
                        for q,r in variants:
                            mm=sum(a!=b for a,b in zip(q,r)); matches=len(q)-mm
                            if mm<=5 and matches>=15:
                                cat=('exact' if mm==0 and delta==0 else 'mismatch_'+str(mm) if delta==0 else ('I' if delta<0 else 'D')+str(abs(delta))+('_mm'+str(mm) if mm else ''))
                                truth.setdefault((gid,contig,strand,anchor),set()).add(cat)
    missing=[d for d in design if (d['guide'],d['contig'],d['strand'],d['anchor']) not in truth]; assert not missing,missing
    names={'mm5':'mismatch_5','I2mm3':'I2_mm3','D2mm5':'D2_mm5'}
    missing_classes=[d for d in design if names.get(d['case'],d['case']) not in truth[(d['guide'],d['contig'],d['strand'],d['anchor'])]]
    assert not missing_classes, missing_classes
    (out/'oracle.jsonl').write_text(''.join(json.dumps({'guide':g,'contig':c,'strand':s,'anchor':p,'categories':sorted(v)},sort_keys=True)+'\n' for (g,c,s,p),v in sorted(truth.items())))
    print(json.dumps({'guides':len(guides),'contigs':len(seqs),'bases':sum(map(len,seqs.values())),'anchors':len(truth),'categories':dict(collections.Counter(c for cats in truth.values() for c in cats))},sort_keys=True))
if __name__=='__main__':main()
