from pathlib import Path
import subprocess,sys,json,platform
ROOT=Path(__file__).resolve().parent
for name in ['analyze.py','plot_figure2.py','plot_figure3.py','plot_figure4.py','plot_figureS1.py','plot_figureS2.py']:
 print('Running',name,flush=True)
 subprocess.run([sys.executable,str(ROOT/'scripts'/name)],check=True,cwd=ROOT)
import numpy,matplotlib
(ROOT/'results/environment.json').write_text(json.dumps({'python':platform.python_version(),'numpy':numpy.__version__,'matplotlib':matplotlib.__version__},indent=2))
print('All analysis and plotting scripts completed. Outputs: results/')
