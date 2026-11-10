#!/usr/bin/env python3
"""CPU-only batched inference for the pinned official ensemble. No pickle loading."""
import argparse
import json
import os
from pathlib import Path
import resource
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from adapter import digest, encode_tokens

MANIFEST = Path(__file__).with_name('model_manifest.json')


def verify_models(root, manifest):
    paths = []
    for spec in manifest['models']:
        path = root/spec['path']
        if digest(path) != spec['sha256']:
            raise ValueError('model checksum mismatch: '+str(path))
        paths.append(path)
    return paths


class Ensemble:
    def __init__(self, root):
        started = time.perf_counter()
        # CPU/thread settings are explicit, not inherited GPU/autothreading defaults.
        os.environ.update(CUDA_VISIBLE_DEVICES='-1', TF_NUM_INTRAOP_THREADS='1',
                          TF_NUM_INTEROP_THREADS='1', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
        os.environ.setdefault('TF_CPP_MIN_LOG_LEVEL', '2')
        import tensorflow as tf
        import numpy as np
        if tf.__version__ != '2.12.0':
            raise ValueError('this model integration is validated with TensorFlow 2.12.0 only')
        self.tf, self.np = tf, np
        self.import_seconds = time.perf_counter()-started
        self.manifest = json.loads(MANIFEST.read_text())
        self.root = Path(root)
        started = time.perf_counter()
        self.paths = verify_models(self.root, self.manifest)
        self.verify_pre_seconds = time.perf_counter()-started
        started = time.perf_counter()
        self.models = [tf.keras.models.load_model(path, compile=False) for path in self.paths]
        for model in self.models:
            if model.input_shape != (None, 24) or model.output_shape != (None, 1):
                raise ValueError('unexpected model shape')
        self.load_seconds = time.perf_counter()-started

    def predict(self, rows):
        started = time.perf_counter()
        features = self.np.asarray(encode_tokens(rows), dtype=self.np.int64)
        encoding_seconds = time.perf_counter()-started
        started = time.perf_counter()
        components = self.np.stack([m.predict(features, batch_size=4096, verbose=0).reshape(-1)
                                   for m in self.models], axis=1)
        predictions = components.mean(axis=1)
        if not self.np.isfinite(predictions).all() or not ((predictions >= 0) & (predictions <= 1)).all():
            raise ValueError('invalid model prediction')
        return predictions, dict(rows=len(rows), encoding_seconds=encoding_seconds,
                                 predict_seconds=time.perf_counter()-started)


def main():
    p = argparse.ArgumentParser()
    for name in ('model-root', 'input', 'output', 'metrics'):
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--batch-size', type=int, default=512)
    a = p.parse_args()
    if not 1 <= a.batch_size <= 4096:
        raise ValueError('batch size must be in 1..4096')
    started = time.perf_counter()
    ensemble = Ensemble(a.model_root)
    timings = []
    with a.input.open() as source, a.output.open('x') as target:
        batch = []
        def flush():
            if not batch:
                return
            scores, elapsed = ensemble.predict(batch)
            for row, score in zip(batch, scores):
                target.write(json.dumps(dict(id=row['id'], score=float(score)), allow_nan=False)+'\n')
            timings.append(elapsed)
            batch.clear()
        for line in source:
            batch.append(json.loads(line))
            if len(batch) == a.batch_size:
                flush()
        flush()
    post = time.perf_counter()
    verify_models(a.model_root, ensemble.manifest)
    result = dict(model_id=ensemble.manifest['model_id'], model_load_count=len(ensemble.models),
        tensorflow=ensemble.tf.__version__, numpy=ensemble.np.__version__, python=sys.version,
        import_seconds=ensemble.import_seconds, verify_pre_seconds=ensemble.verify_pre_seconds,
        model_load_seconds=ensemble.load_seconds, verify_post_seconds=time.perf_counter()-post,
        batches=timings, elapsed_seconds=time.perf_counter()-started,
        peak_worker_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        scope='CPU worker peak, not aggregate memory; includes TF and models; not search runtime')
    with a.metrics.open('x') as f:
        json.dump(result, f, indent=2)
        f.write('\n')


if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print('CRISPR-Bulge worker failed: '+str(e), file=sys.stderr)
        sys.exit(1)
