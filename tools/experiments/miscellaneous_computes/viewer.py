"""Build native Viewer dashboards from time series or phase curves and render selected pages."""
import json
import subprocess
import sys
from pathlib import Path
import numpy as np

experiments=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(experiments/'trajectory_resolution'))
from resolution import DisplayVertices

def dashboard(folder, pages, title, source, selector_library, tolerance=3e-5):
    """pages maps a tab label to panels. Each curve is (name, x, y, color, style, width)."""
    folder=Path(folder).resolve()
    folder.mkdir(parents=True,exist_ok=True)
    selector=DisplayVertices(selector_library)
    figures=[]
    for page,panels in enumerate(pages.values()):
        for panel in panels:
            traces=[]
            for label,x,y,color,style,width in panel['curves']:
                x,y=np.asarray(x,dtype=float),np.asarray(y,dtype=float)
                if x.ndim!=1 or y.shape!=x.shape or not len(x) or not np.all(np.isfinite(x+y)):
                    raise ValueError('Each trace requires paired finite coordinates')
                keep=selector(x,y,tolerance) if len(x)>4096 else np.arange(len(x))
                traces.append(dict(label=label,x=x[keep].tolist(),y=y[keep].tolist(),color=color,
                    line_style=style,line_width_points=width,source_samples=len(x),displayed_vertices=len(keep)))
            figures.append(dict(kind=panel.get('kind','phase'),title=panel['title'],x_label=panel['x_label'],
                         y_label=panel['y_label'],value_scale='linear',dashboard_page=page,traces=traces))
    payload=dict(schema='dynsyslab.viewer.payload/v1',source=str(source),view='phase',
                 tile=dict(name=title,columns=2),page_labels=list(pages),figures=figures)
    catalog=dict(schema='dynsyslab.viewer.catalog/v1',items=[dict(id='study',label=title,kind='trajectory',
                 source='payload.json',views=['phase'])])
    for name,obj in [('payload',payload),('catalog',catalog)]:
        (folder/(name+'.json')).write_text(json.dumps(obj,allow_nan=False),encoding='utf-8')
    return folder/'payload.json'

def render(executable, payload, page, output, width=1250,height=720):
    """Render one dashboard page through the installed native Viewer CLI using a page-specific cache."""
    payload=Path(payload).resolve()
    folder=payload.parent
    (folder/'archive').mkdir(exist_ok=True)
    subprocess.run([str(executable),'render','--source',str(payload),'--catalog',str(folder/'catalog.json'),
        '--project',str(folder),'--archives',str(folder/'archive'),'--cache',str(folder/f'cache-{page}'),
        '--page',str(page),'--width',str(width),'--height',str(height),'--output',str(Path(output).resolve())],check=True)
