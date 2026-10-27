#!/usr/bin/env python3
"""Build the self-contained article with restricted, deterministic pdfLaTeX."""
from pathlib import Path
import hashlib
import json
import os
import shutil
import subprocess

root = Path(__file__).resolve().parent
paper = root / 'paper'
build = paper / 'build'
build.mkdir(exist_ok=True)
env = os.environ.copy()
env.update(SOURCE_DATE_EPOCH='1792108800', FORCE_SOURCE_DATE='1', TZ='UTC')
command = ['pdflatex', '-no-shell-escape', '-interaction=nonstopmode',
           '-halt-on-error', '-file-line-error', '-recorder',
           '-output-directory=build', 'main.tex']
for number in range(1, 4):
    result = subprocess.run(command, cwd=paper, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (build / f'pass-{number}.txt').write_bytes(result.stdout)
    if result.returncode:
        print(result.stdout.decode(errors='replace'))
        raise SystemExit(result.returncode)
log = (build / 'main.log').read_text(errors='replace')
errors = [line for line in log.splitlines() if
          'undefined references' in line or 'undefined citations' in line or
          'multiply defined' in line or 'Rerun to get cross-references right' in line]
if errors:
    raise SystemExit('\n'.join(errors))
shutil.copyfile(build / 'main.pdf', paper / 'main.pdf')
pdf = (paper / 'main.pdf').read_bytes()
info = subprocess.check_output(['pdfinfo', str(paper / 'main.pdf')], text=True)
pages = int(next(line.split(':', 1)[1] for line in info.splitlines()
                 if line.startswith('Pages:')))
receipt = {'command': command, 'passes': 3,
           'source_date_epoch': env['SOURCE_DATE_EPOCH'],
           'engine': subprocess.check_output(['pdflatex', '--version'], text=True).splitlines()[0],
           'pages': pages, 'pdf_bytes': len(pdf),
           'pdf_sha256': hashlib.sha256(pdf).hexdigest()}
(root / 'BUILD_RECEIPT.json').write_text(json.dumps(receipt, indent=2) + '\n')
print(json.dumps(receipt, indent=2))
