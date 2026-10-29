import sys,subprocess,tempfile,pathlib,json
sys.path.insert(0,'scripts/release_validation')
from coverage_oracle import enumerate_sites
binary='/export/local/home/salehi/crisprapido-pam-target/release/crisprapido'
with tempfile.TemporaryDirectory() as d:
 p=pathlib.Path(d)
 g='ACGTCAGTACGATCGTACGA'
 # Deliberately isolated 10-base example: four substitutions + one missing guide base.
 short='ACGTCAGTAC'; target='TGGTAGTAC'
 records={'twenty':(g,'TCGTAAGTACGATAGTACGTAGG'),'short':(short,target+'AGG')}
 fasta='/export/local/home/salehi/crisprapido-index-recovery-v1-2026-10-28/recovered_text_CHM13/input/chm13v2.0_maskedY_rCRS.fa.PanSN.fa.gz'
 samtools='/export/local/home/salehi/verkko_env/bin/samtools'
 raw=subprocess.check_output([samtools,'faidx',fasta,'CHM13#0#chr1:100001-100500'],text=True)
 seq=''.join(raw.splitlines()[1:]).upper()
 for i in range(len(seq)-23):
  if set(seq[i:i+23])<=set('ACGT') and seq[i+21:i+23]=='GG':
   records['real']=(seq[i:i+20],seq);break
 for label,(guide,seq) in records.items():
  truth=list(enumerate_sites(seq,guide,4,True,0))
  keys=sorted({(x['start'],x['end'],x['strand']) for x in truth})
  (p/'ref.fa').write_text('>chr\n'+seq+'\n')
  for f in (.75,0):
   rows=[]
   for start,end,strand in keys:
    # Supply each oracle locus independently, avoiding overlap selection between loci.
    (p/'in.sam').write_text(f'guide\t{16 if strand=="-" else 0}\tchr\t{start+1}\t60\t{end-start}M\t*\t0\t0\t*\t*\tAS:i:0\tNM:i:0\n')
    proc=subprocess.run([binary,'-r',str(p/'ref.fa'),'-g',guide,'--columba-sam',str(p/'in.sam'),'-m','4','-b','1','-z','1','-f',str(f)],capture_output=True,text=True,check=True)
    if proc.stdout: rows.append({'input':[start,end,strand],'paf':proc.stdout.strip().splitlines()})
   print(json.dumps({'case':label,'guide':guide,'f':f,'truth_loci':len(keys),'rows':rows}))
