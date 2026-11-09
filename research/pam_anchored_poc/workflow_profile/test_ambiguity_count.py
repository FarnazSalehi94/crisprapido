#!/usr/bin/env python3
"""Differential test for the verifier's one-pass ambiguous-byte count."""
import random
r=random.Random(567839);cases=[b'',bytes(range(256)),b'ACGT',b'NNNN',b'ACGTNacgtn',b'\x00\xffACGT']
for n in range(512):cases.append(bytes(r.randrange(256) for _ in range(n)))
for s in cases:
 old=len(s)-sum(s.count(bytes([b])) for b in b'ACGT');new=len(s.translate(None,b'ACGT'));assert old==new,(s,old,new)
print('ambiguous-byte count matches prior algorithm on 518 cases including all byte values')
