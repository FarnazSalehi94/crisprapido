#!/usr/bin/env python3
"""Run the unmodified official ensemble_predict on identical aligned inputs.

This reference-only path loads upstream trusted pickle files. The production-
independent adapter worker deliberately does not unpickle models/import training.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from adapter import digest


def main():
    p=argparse.ArgumentParser()
    for key in ('upstream','input','output','metrics'):
        p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();start=time.perf_counter()
    a.upstream=a.upstream.resolve();a.input=a.input.resolve()
    a.output=a.output.resolve();a.metrics=a.metrics.resolve()
    os.environ.update(CUDA_VISIBLE_DEVICES='-1', TF_CPP_MIN_LOG_LEVEL='2', TF_NUM_INTRAOP_THREADS='1',
                      TF_NUM_INTEROP_THREADS='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
    spec=json.loads(Path(__file__).with_name('model_manifest.json').read_text())
    commit=subprocess.check_output(['git','-C',a.upstream,'rev-parse','HEAD'],text=True).strip()
    if commit!=spec['commit']:
        raise ValueError('wrong upstream commit')
    for name,value in spec['source_sha256'].items():
        if digest(a.upstream/name)!=value:raise ValueError('changed upstream source: '+name)
    if digest(a.upstream/'LICENSE')!=spec['license_sha256']:raise ValueError('wrong upstream license')
    for model in spec['models']:
        for path,key in [('path','sha256'),('pickle_path','pickle_sha256')]:
            if digest(a.upstream/model[path])!=model[key]:raise ValueError('changed upstream model')
    # Upstream imports third-party auxiliary model modules; pin those too before
    # executing their imports (Nuclea-seq reads/evaluates its parameter files).
    for path in (a.upstream/'models').rglob('*'):
        if path.suffix not in ('.py','.txt'):continue
        rel=str(path.relative_to(a.upstream))
        expected=subprocess.check_output(['git','-C',a.upstream,'rev-parse','HEAD:'+rel],text=True).strip()
        data=path.read_bytes()
        actual=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
        if expected!=actual:raise ValueError('changed auxiliary upstream source')
    sys.path.insert(0,str(a.upstream.resolve()))
    import pandas as pd
    import numpy as np
    import tensorflow as tf
    from train_and_predict_scripts.utilities import ensemble_predict
    from OT_deep_score_src.data_processing_utilities import build_sequence_features
    from OT_deep_score_src.general_utilities import Padding_type
    rows=[json.loads(line) for line in a.input.open()]
    if not rows:raise ValueError('reference requires at least one eligible input')
    frame=pd.DataFrame({'Align.sgRNA':[r['guide'] for r in rows],
                        'Align.off-target':[r['target'] for r in rows]})
    features=build_sequence_features(frame.copy(),bulges=True,padding_type=Padding_type.GAP,
                                     aligned=True,flat_encoding=False)
    categories=np.argmax(features,axis=2).astype('<i8')
    feature_hash=hashlib.sha256(categories.tobytes()).hexdigest()
    paths=[str(a.upstream/m['path'])[:-3] for m in spec['models']]
    predict_start=time.perf_counter()
    os.chdir(a.output.parent)  # upstream writes predictions.csv; never write it into CRISPRapido
    result=ensemble_predict(paths,frame)
    predict_seconds=time.perf_counter()-predict_start
    with a.output.open('x') as f:
        for row,score in zip(rows,result['pred_averege_ensemble'].values):
            f.write(json.dumps(dict(id=row['id'],score=float(score)),allow_nan=False)+'\n')
    # Characterize encoder limits, not unsupported model predictions.
    g='GAAACAGTCGATTTTATCAC'
    dna2_error=None
    try:
        build_sequence_features(pd.DataFrame({'Align.sgRNA':[g[:5]+'--'+g[5:]+'NGG'],
            'Align.off-target':[g[:5]+'TT'+g[5:]+'AGG']}),bulges=True,padding_type=Padding_type.GAP,flat_encoding=False)
    except ValueError as e:dna2_error=str(e)
    rna2=build_sequence_features(pd.DataFrame({'Align.sgRNA':[g+'NGG'],
        'Align.off-target':[g[:5]+'--'+g[7:]+'AGG']}),bulges=True,padding_type=Padding_type.GAP,flat_encoding=False)
    metrics=dict(official_commit=commit,rows=len(rows),tensorflow=tf.__version__,numpy=np.__version__,
        categorical_sha256=feature_hash,official_ensemble_call_seconds=predict_seconds,
        elapsed_seconds=time.perf_counter()-start,peak_reference_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        dna_two_base_encoder_error=dna2_error,rna_two_base_encoder_shape=list(rna2.shape),
        scope='unmodified official function plus imports, hash validation and feature/limit probes; CPU')
    with a.metrics.open('x') as f:json.dump(metrics,f,indent=2);f.write('\n')


if __name__=='__main__':main()
