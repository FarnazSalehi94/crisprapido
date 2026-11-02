"""Check real-panel provenance, all planted structures, and PAM-key comparison."""
import hashlib
import json
from pathlib import Path
import unittest
from audit_emitted import load_guides, parse_row, rc

ROOT=Path(__file__).resolve().parents[2]/'docs/real_defaults_validation_49f20fc_2026-11-02_v1'
def read_fasta(path):
    seqs={};name=None
    for line in path.read_text().splitlines():
        if line.startswith('>'):name=line[1:];seqs[name]=''
        else:seqs[name]+=line
    return seqs

class RealDefaultsPanelTests(unittest.TestCase):
    def test_reference_provenance_and_designed_alignments(self):
        panel=json.loads((ROOT/'panel.json').read_text());refs=read_fasta(ROOT/'ref.fa');guides=load_guides(ROOT/'guides.fa')
        for window in panel['windows']:
            seq=refs[window['name']];self.assertEqual(len(seq),2000)
            self.assertEqual(hashlib.sha256(seq.encode()).hexdigest(),window['sha256'])
        self.assertEqual(len(panel['designs']),96)
        for d in panel['designs']:
            ref=refs[d['contig']];strand=d['strand'];ori=ref if strand=='+' else rc(ref)
            end=d['anchor'] if strand=='+' else len(ref)-d['anchor']-3
            self.assertEqual(ori[end:end+3],d['pam']);self.assertEqual(d['pam'][1:],'GG')
            case=d['case'];guide=guides[d['guide']]
            if case.startswith('mm'):
                target=ori[end-20:end];self.assertEqual(sum(a!=b for a,b in zip(guide,target)),int(case[2:]))
            elif case=='exact':self.assertEqual(guide,ori[end-20:end])
            else:
                direction=case[0];gap=int(case[1]);target=ori[end-(20-gap if direction=='I' else 20+gap):end]
                comparisons=([guide[:i]+guide[i+gap:] for i in range(21-gap)] if direction=='I' else [target[:i]+target[i+gap:] for i in range(21)])
                other=target if direction=='I' else guide
                want=2 if 'mm2' in case else 0
                self.assertTrue(any(sum(a!=b for a,b in zip(candidate,other))==want for candidate in comparisons),d)
    def test_reported_keys_and_cigars(self):
        refs=read_fasta(ROOT/'ref.fa');guides=load_guides(ROOT/'guides.fa')
        truth={(x['guide'],x['contig'],x['strand'],x['anchor']) for x in map(json.loads,(ROOT/'oracle.jsonl').read_text().splitlines())}
        found=[]
        for line in (ROOT/'actual.paf').read_text().splitlines():
            f,_,_,_,_=parse_row(line,guides,{k:len(v) for k,v in refs.items()})
            found.append((f[0],f[5],f[4],int(f[8]) if f[4]=='+' else int(f[7])-3))
        self.assertEqual(len(truth),204);self.assertEqual(len(found),len(set(found)))
        self.assertEqual(set(found),truth)
        audit=[json.loads(x) for x in (ROOT/'audit.jsonl').read_text().splitlines()]
        self.assertEqual(len(audit),len(found));self.assertTrue(all(x['pam'][1:]=='GG' for x in audit))

if __name__=='__main__':unittest.main()
