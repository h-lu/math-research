#!/usr/bin/env python3
"""Build this independent candidate with installed local TeX tools."""
from pathlib import Path
import os,subprocess,shutil,sys
root=Path(__file__).resolve().parent
build=root/'build'
build.mkdir(exist_ok=True)
env=dict(os.environ)
env['BIBINPUTS']=str(root)+os.pathsep+env.get('BIBINPUTS','')
env['TEXINPUTS']=str(root)+os.pathsep+env.get('TEXINPUTS','')
subprocess.run(['latexmk','-pdf','-interaction=nonstopmode','-halt-on-error',
                '-outdir='+str(build),'main.tex'],cwd=root,env=env,check=True)
shutil.copyfile(build/'main.pdf',root/'DF1b-exact2.pdf')
subprocess.run([sys.executable,str(root/'tools/sanitize-pdf-metadata.py'),str(root/'DF1b-exact2.pdf')],check=True)
print('PDF:',root/'DF1b-exact2.pdf')
