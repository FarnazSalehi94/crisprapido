#!/usr/bin/env python3
"""Check Columba vanilla 32-bit index text, sequence names, and offsets against FASTA.

The index replaces ambiguous input bases, so identity at those positions cannot
be established from the index. The returned ambiguity count makes that explicit.
"""
import hashlib,json,struct,sys
from pathlib import Path

def fasta(p):
    # Build each record once: repeated string concatenation is quadratic for
    # wrapped chromosome-scale FASTA. Preserve the existing parser's rules.
    records=[]
    with Path(p).open() as source:
        for line in source:
            if line.startswith('>'):
                records.append([line[1:].split()[0], []])
            else:
                if not records: raise ValueError('invalid FASTA')
                records[-1][1].append(line.strip().upper())
    seqs=[[name, ''.join(chunks)] for name,chunks in records]
    if not seqs or len({n for n,_ in seqs})!=len(seqs) or any(not s for _,s in seqs): raise ValueError('empty or duplicate FASTA records')
    return seqs

def check(prefix,ref):
    """Stream the 32-bit indexed text and FASTA; never materialize a genome."""
    prefix=str(prefix); names=[]; lengths=[]; ambiguous=0; total=0
    ref_hash=hashlib.sha256(); text_hash=hashlib.sha256()
    with Path(ref).open('rb') as source, Path(prefix+'.txt.bin').open('rb') as index:
        raw=index.read(4)
        if len(raw)!=4: raise ValueError('short text file')
        length=struct.unpack('<I',raw)[0]
        if Path(prefix+'.txt.bin').stat().st_size!=length+4: raise ValueError('index text length/sentinel mismatch')
        for line in source:
            ref_hash.update(line)
            if line.startswith(b'>'):
                fields=line[1:].split()
                if not fields: raise ValueError('empty FASTA record ID')
                names.append(fields[0].decode('ascii')); lengths.append(0)
                continue
            if not names: raise ValueError('invalid FASTA')
            seq=line.strip().upper(); indexed=index.read(len(seq))
            if len(indexed)!=len(seq): raise ValueError('index text length/sentinel mismatch')
            text_hash.update(indexed); lengths[-1]+=len(seq); total+=len(seq)
            ambiguous+=len(seq)-sum(seq.count(bytes([b])) for b in b'ACGT')
            if seq!=indexed and any(a!=b for a,b in zip(seq,indexed) if a in b'ACGT'):
                raise ValueError('index/FASTA resolved base mismatch')
        if not names or len(set(names))!=len(names) or any(n==0 for n in lengths): raise ValueError('empty or duplicate FASTA records')
        sentinel=index.read(1);text_hash.update(sentinel)
        if sentinel!=b'$' or index.read(1) or length!=total+1: raise ValueError('index text length/sentinel mismatch')
    pos=Path(prefix+'.pos').read_bytes()
    if len(pos)%4: raise ValueError('invalid position array')
    actual=list(struct.unpack('<'+'I'*(len(pos)//4),pos)); expected=[0]
    for n in lengths: expected.append(expected[-1]+n)
    if actual!=expected: raise ValueError('index/FASTA contig offsets differ')
    raw=Path(prefix+'.sna').read_bytes(); off=0; index_names=[]
    while off<len(raw):
        if off+8>len(raw): raise ValueError('truncated name length')
        n=struct.unpack_from('<Q',raw,off)[0]; off+=8
        if off+n>len(raw): raise ValueError('truncated name')
        index_names.append(raw[off:off+n].decode()); off+=n
    if index_names!=names: raise ValueError('index/FASTA names or order differ')
    return {'fasta_sha256':ref_hash.hexdigest(),'indexed_text_sha256':text_hash.hexdigest(),'bases':total,'contigs':len(names),'ambiguous_positions_unverifiable':ambiguous,'resolved_base_mismatches':0,'scope':'32-bit little-endian vanilla index; does not prove BWT/SA integrity or original ambiguous bases'}
if __name__=='__main__': print(json.dumps(check(sys.argv[1],sys.argv[2]),sort_keys=True))
