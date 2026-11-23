"""Optional local notebook packaging; no scientific computation or server setup."""
from __future__ import annotations
from html import escape
from pathlib import Path


def standalone_fragment(fragment, *, title="Interactive figure"):
    """Give an existing Viewer fragment a self-contained light-theme document."""
    style = """
    :root {--font-size-base:16px;--foreground:#17232e;--background:#fff;
      --border:#d8dee5;--red:#d62728;--green:#00865a;--blue:#0072b2;
      --viz-series-1:#17232e; color:#17232e;background:#fff;
      font:16px system-ui,sans-serif}
    body{margin:16px}.viz-controls{display:flex;gap:20px;flex-wrap:wrap}
    .form-label{display:grid;gap:5px}.form-select{min-width:80px;font:inherit}
    .viz-row{display:flex;gap:10px;align-items:center;flex-wrap:wrap}
    .text-small{font-size:13px}.tabular-nums{font-variant-numeric:tabular-nums}
    """
    return ("<!doctype html><html><head><meta charset='utf-8'>"
            f"<title>{escape(title)}</title><style>{style}</style></head>"
            f"<body>{fragment}</body></html>")


def notebook_html_frame(document, *, height=780, title="Interactive figure"):
    """Embed a document without leaking its CSS into the notebook interface.

    Trust the locally generated notebook to activate HTML/JavaScript outputs.
    The document has no server dependency; restrictive notebook hosts may still
    disable scripts. An exported standalone document is a useful fallback.
    """
    if not isinstance(height, int) or not 100 <= height <= 10000:
        raise ValueError("Provide a frame height between 100 and 10000 pixels")
    return (f'<iframe title="{escape(title, quote=True)}" '
            f'srcdoc="{escape(document, quote=True)}" '
            f'style="width:100%;height:{height}px;border:1px solid #d8dee5" '
            'sandbox="allow-scripts"></iframe>')


def write_executed_notebook(cells, path, *, title="DynSysLab analysis",
                            kernel_name="python3", timeout=120):
    """Execute caller-authored cells locally and persist outputs; errors raise.

    Requires optional nbformat/nbclient in the calling Python environment.
    This does not start a Jupyter server, install packages, or upload content.
    """
    import nbformat
    from nbclient import NotebookClient
    path = Path(path).resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    notebook = nbformat.v4.new_notebook(cells=cells, metadata={
        "title": title,
        "kernelspec": {"display_name": "Python 3", "language": "python", "name": kernel_name},
        "language_info": {"name": "python"},
    })
    NotebookClient(notebook, kernel_name=kernel_name, timeout=timeout,
                   resources={"metadata": {"path": str(path.parent)}}).execute()
    nbformat.validate(notebook)
    nbformat.write(notebook, path)
    return path
