import json,random,subprocess,tempfile,unittest
from pathlib import Path

HERE=Path(__file__).resolve().parent
SAMTOOLS='/export/local/home/salehi/verkko_env/bin/samtools'
class CandidateGeometryTests(unittest.TestCase):
 def test_strict_pam_geometry_against_brute_force(self):
  rng=random.Random(281);seq=''.join(rng.choice('ACGTN') for _ in range(450))
  with tempfile.TemporaryDirectory() as temp:
   p=Path(temp);ref=p/'ref.fa';ref.write_text('>chr\n'+seq+'\n')
   subprocess.run([SAMTOOLS,'faidx',str(ref)],check=True)
   rows=[];unique=direct=feasible=0
   for pos in list(range(0,28))+list(range(50,300,7))+list(range(420,430)):
    if pos+20>len(seq):continue
    for rev in (False,True):
     flag=16 if rev else 0;row=f'g\t{flag}\tchr\t{pos+1}\t60\t20M\t*\t0\t0\t*\t*\tAS:i:0\tNM:i:0\n';rows.append(row)
     if pos==0 and not rev:rows.append(row)
     unique+=1;end=pos+20
     if rev:good=pos>=3 and seq[pos-3:pos-1]=='CC' and seq[pos-1] in 'ACGT'
     else:good=end+3<=len(seq) and seq[end+1:end+3]=='GG' and seq[end] in 'ACGT'
     direct+=good;possible=False
     for start in range(len(seq)):
      for length in range(18,23):
       stop=start+length
       if stop>len(seq) or not (start<end and stop>pos):continue
       if rev:ok=start>=3 and seq[start-3:start-1]=='CC' and seq[start-1] in 'ACGT'
       else:ok=stop+3<=len(seq) and seq[stop+1:stop+3]=='GG' and seq[stop] in 'ACGT'
       possible|=ok
     feasible+=possible
   sam=p/'in.sam';sam.write_text(''.join(rows));exact=p/'exact.sam';exact.write_text(rows[0])
   subprocess.run(['python3',str(HERE/'analyze_pam_candidate_geometry.py'),'--sam',str(sam),'--exact-sam',str(exact),'--fasta',str(ref),'--samtools',SAMTOOLS,'--output',str(p/'out')],check=True,stdout=subprocess.DEVNULL)
   stats=json.loads((p/'out/summary.json').read_text());self.assertEqual(stats['mapped_raw'],unique+1);self.assertEqual(stats['mapped_unique'],unique)
   self.assertEqual(stats['duplicate_keys'],1);self.assertEqual(stats['exact_pass_keys_already_in_e7'],1)
   self.assertEqual(stats['unique']['direct_pam_valid'],direct);self.assertEqual(stats['unique']['no_possible_overlapping_pam'],unique-feasible)
if __name__=='__main__':unittest.main()
