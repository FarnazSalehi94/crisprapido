#!/usr/bin/env python3
"""Check Columba vanilla 32-bit index text, sequence names, and offsets against FASTA.

The index replaces ambiguous input bases, so identity at those positions cannot
be established from the index. The returned ambiguity count makes that explicit.
"""
import hashlib,json,struct,sys
from pathlib import Path

def fasta(p):
    seqs=[]
    for line in Path(p).read_text().splitlines():
        if line.startswith('>'): seqs.append([line[1:].split()[0],''])
        else:
            if not seqs: raise ValueError('invalid FASTA')
            seqs[-1][1]+=line.strip().upper()
    if not seqs or len({n for n,_ in seqs})!=len(seqs) or any(not s for _,s in seqs): raise ValueError('empty or duplicate FASTA records')
    return seqs

def check(prefix,ref):
    prefix=str(prefix); seqs=fasta(ref); concat=''.join(s for _,s in seqs).encode('ascii'); blob=Path(prefix+'.txt.bin').read_bytes()
    if len(blob)<4: raise ValueError('short text file')
    length=struct.unpack_from('<I',blob)[0]; text=blob[4:]
    if length!=len(text) or text!=text[:len(concat)]+b'$' or length!=len(concat)+1: raise ValueError('index text length/sentinel mismatch')
    bad=sum(a!=b for a,b in zip(concat,text) if a in b'ACGT')
    if bad: raise ValueError('index/FASTA resolved base mismatch: '+str(bad))
    pos=Path(prefix+'.pos').read_bytes(); actual=list(struct.unpack('<'+'I'*(len(pos)//4),pos))
    expected=[0];
    for _,s in seqs: expected.append(expected[-1]+len(s))
    if actual!=expected: raise ValueError('index/FASTA contig offsets differ')
    raw=Path(prefix+'.sna').read_bytes(); off=0; names=[]
    while off<len(raw):
        if off+8>len(raw): raise ValueError('truncated name length')
        n=struct.unpack_from('<Q',raw,off)[0]; off+=8
        if off+n>len(raw): raise ValueError('truncated name')
        names.append(raw[off:off+n].decode()); off+=n
    if names!=[n for n,_ in seqs]: raise ValueError('index/FASTA names or order differ')
    return {'fasta_sha256':hashlib.sha256(Path(ref).read_bytes()).hexdigest(),'indexed_text_sha256':hashlib.sha256(text).hexdigest(),'bases':len(concat),'contigs':len(seqs),'ambiguous_positions_unverifiable':sum(c not in b'ACGT' for c in concat),'resolved_base_mismatches':bad,'scope':'32-bit little-endian vanilla index; does not prove BWT/SA integrity or original ambiguous bases'}
if __name__=='__main__': print(json.dumps(check(sys.argv[1],sys.argv[2]),sort_keys=True))
