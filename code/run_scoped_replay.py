"""Run and compare the bounded frozen-input replay; never retrain any model."""
from pathlib import Path
import argparse, csv, hashlib, importlib.metadata, json, os, platform, subprocess, sys, time


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest().upper()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument('--out', type=Path, required=True, help='New directory outside the package')
    p.add_argument('--rscript', required=True, help='Rscript executable, preferably R 4.5.2')
    p.add_argument('--r-library', type=Path, help='Explicit installed R package library')
    a = p.parse_args()
    root, out = a.root.resolve(), a.out.resolve()
    if out.is_relative_to(root) or root.is_relative_to(out) or out.exists():
        p.error('--out must be a new directory outside and not containing the package')
    out.mkdir(parents=True)
    env = os.environ.copy()
    env.update(PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1',
               MKL_NUM_THREADS='1', MPLBACKEND='Agg', MPLCONFIGDIR=str(out/'mplconfig'),
               TEMP=str(out/'tmp'), TMP=str(out/'tmp'))
    (out/'tmp').mkdir()
    if a.r_library:
        env['R_LIBS_USER'] = str(a.r_library.resolve())
    report = dict(status='RUNNING', scope='Frozen scoring, clinical tables, archived Fig6/7 panel rendering',
                  fresh_full_environment_install=False, scientific_end_to_end=False,
                  model_training=False, bootstrap_refits=0, random_pool_refits=0,
                  external_random_null_evaluated=False, uploaded=False, checks=[])

    def run(name, cmd):
        start = time.time()
        v = subprocess.run([str(x) for x in cmd], cwd=out, env=env, capture_output=True,
                           text=True, encoding='utf-8', errors='replace')
        (out/(name+'.log')).write_text(v.stdout+'\nSTDERR\n'+v.stderr, encoding='utf-8')
        report['checks'].append(dict(name=name, returncode=v.returncode,
                                    passed=v.returncode == 0, seconds=round(time.time()-start, 3)))
        if v.returncode:
            raise RuntimeError(name+' failed; inspect '+name+'.log')
        return v.stdout

    try:
        runtime = dict(python=sys.version, executable=sys.executable, platform=platform.platform(),
                       rscript=a.rscript, r_library=env.get('R_LIBS_USER'), packages={})
        for line in (root/'requirements-lock.txt').read_text().splitlines():
            if '==' not in line:
                continue
            name, required = line.strip().split('==')
            actual = importlib.metadata.version(name)
            runtime['packages'][name] = dict(expected=required, actual=actual, match=actual == required)
        (out/'python_environment.json').write_text(json.dumps(runtime, indent=2), encoding='utf-8')
        if not all(x['match'] for x in runtime['packages'].values()):
            raise RuntimeError('Python versions differ from the scoped replay lock')
        run('R_environment', [a.rscript, '--vanilla', root/'code/check_R_environment.R',
                              root/'environments/R_packages_expected.tsv'])
        run('input_hashes_before', [sys.executable, root/'code/check_inputs.py', '--root', root])
        run('frozen_scoring', [sys.executable, root/'code/reproduce_frozen_scores.py',
                               '--root', root, '--out', out/'scoring'])
        run('clinical', [sys.executable, root/'code/reproduce_clinical.py', '--root', root,
                         '--out', out/'clinical', '--rscript', a.rscript])
        for fig in ['Fig6', 'Fig7']:
            run('render_'+fig, [sys.executable, root/f'figures/{fig}/renderers/render_{fig}_panels.py',
                                 '--output', out/fig])
        for ref in sorted((root/'source_data/clinical').glob('[1][0-5]_*.tsv')):
            actual = out/'clinical/outputs'/ref.name
            match = sha(actual) == sha(ref)
            report['checks'].append(dict(name='clinical_table_identity', file=ref.name,
                                        passed=match, sha256=sha(actual)))
        from PIL import Image
        import numpy as np
        for fig in ['Fig6', 'Fig7']:
            for ref in sorted((root/f'figures/{fig}/panels').glob(fig+'?.png')):
                with Image.open(ref) as x, Image.open(out/fig/ref.name) as y:
                    match = np.array_equal(np.asarray(x.convert('RGBA')), np.asarray(y.convert('RGBA')))
                report['checks'].append(dict(name='archived_panel_pixel_identity', file=ref.name,
                                            passed=bool(match),
                                            operation='historical image copy' if ref.stem == 'Fig7c' else 'render frozen aggregates'))
        run('input_hashes_after', [sys.executable, root/'code/check_inputs.py', '--root', root])
        if sum(x['name']=='clinical_table_identity' for x in report['checks']) != 6:
            raise RuntimeError('Expected six clinical reference tables')
        if sum(x['name']=='archived_panel_pixel_identity' for x in report['checks']) != 9:
            raise RuntimeError('Expected nine archived panels')
        if not all(x['passed'] for x in report['checks']):
            raise RuntimeError('Output comparison failed')
        report['scoring'] = json.loads((out/'scoring/validation.json').read_text())
        report['status'] = 'PASS_SCOPED_REPLAY'
    except Exception as e:
        report['status'] = 'FAILED'
        report['error'] = str(e)
    finally:
        (out/'replay_report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'status':report['status'], 'report':str(out/'replay_report.json')}))
    return 0 if report['status'] == 'PASS_SCOPED_REPLAY' else 1


if __name__ == '__main__':
    sys.exit(main())
