#!/usr/bin/env python3
"""Checked research-only entry point; no production backend integration."""
import argparse,math,subprocess,sys
from index_manifest import verify
p=argparse.ArgumentParser();p.add_argument('--poc',required=True);p.add_argument('-r','--reference',required=True);p.add_argument('-G','--guides',required=True);p.add_argument('--index',required=True);p.add_argument('-p','--pam',default='GG');p.add_argument('-m',type=int,default=5);p.add_argument('-b',type=int,default=1);p.add_argument('-z',type=int,default=2);p.add_argument('-f',type=float,default=.75)
a=p.parse_args()
try:
 if a.m<0 or a.b not in (0,1) or a.z not in (0,1,2):raise ValueError('unsupported limits: m >= 0, b in {0,1}, z in {0,1,2} required')
 if not math.isfinite(a.f) or not 0<=a.f<=1:raise ValueError('fraction must be finite and between 0 and 1')
 pam=a.pam.upper()
 if len(pam)!=2 or set(pam)-set('ACGT'):raise ValueError('PAM must be exactly two A/C/G/T bases')
 verify(a.reference,a.index)
 sys.exit(subprocess.call([a.poc,a.index,a.reference,a.guides,str(a.m),str(a.z),str(a.f),pam,str(a.b)]))
except (OSError,ValueError,KeyError,TypeError) as e:print('checked search error: '+str(e),file=sys.stderr);sys.exit(2)
