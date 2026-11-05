"""Compact local GUI for configuring and browsing dynsyslab runs."""
from __future__ import annotations
import json
from dataclasses import asdict
from threading import Thread
import webbrowser
from pathlib import Path
import tkinter as tk
from tkinter import ttk, messagebox
from app import DynsysLabService, RunConfiguration
from archiver.gateway import VIEWER_ARCHIVE
from job_manager import JobManager, JobManifest
from postprocessing.koopman import dmd_archive, edmd_archive, hankel_dmd_archive
from rom import sindy_archive
from viewer import ArchiveSource, PlotRequest, TilePlan, Viewer

class DynsysLabDesktop(ttk.Frame):
    def __init__(self,master=None):
        super().__init__(master,padding=12); self.grid(sticky="nsew"); master.columnconfigure(0,weight=1); master.rowconfigure(0,weight=1)
        self.service=DynsysLabService(); self.jobs=JobManager(); self.active=None
        systems=self.service.available_systems(); self.system=tk.StringVar(value=next(iter(systems))); self.preset=tk.StringVar(); self.start=tk.StringVar(value="0"); self.stop=tk.StringVar(value="1"); self.step=tk.StringVar(value="0.01"); self.status=tk.StringVar(value="Ready")
        form=ttk.LabelFrame(self,text="Solve configuration",padding=10); form.grid(row=0,column=0,sticky="ew"); form.columnconfigure(1,weight=1)
        for row,(label,variable) in enumerate((("System",self.system),("Preset",self.preset),("Start",self.start),("Stop",self.stop),("Uniform output step",self.step))):
            ttk.Label(form,text=label).grid(row=row,column=0,sticky="w",padx=(0,8),pady=3)
            if label in {"System","Preset"}: widget=ttk.Combobox(form,textvariable=variable,state="readonly")
            else: widget=ttk.Entry(form,textvariable=variable)
            widget.grid(row=row,column=1,sticky="ew",pady=3)
            if label=="System": self.system_box=widget; widget["values"]=tuple(systems); widget.bind("<<ComboboxSelected>>",lambda _:self._presets())
            if label=="Preset": self.preset_box=widget
        self._presets(); buttons=ttk.Frame(form); buttons.grid(row=5,column=0,columnspan=2,sticky="ew",pady=(8,0)); ttk.Button(buttons,text="Run",command=self.run).pack(side="left"); ttk.Button(buttons,text="Cancel",command=self.cancel).pack(side="left",padx=6)
        ttk.Label(self,textvariable=self.status).grid(row=1,column=0,sticky="w",pady=8)
        browser=ttk.LabelFrame(self,text="Archive, postprocessing, and Viewer",padding=8); browser.grid(row=2,column=0,sticky="nsew"); browser.columnconfigure(0,weight=1); self.rowconfigure(2,weight=1)
        self.archives=tk.Listbox(browser,height=12); self.archives.grid(row=0,column=0,sticky="nsew"); browser.rowconfigure(0,weight=1)
        actions=ttk.Frame(browser); actions.grid(row=1,column=0,sticky="ew",pady=(6,0)); self.analysis=tk.StringVar(value="DMD")
        ttk.Combobox(actions,textvariable=self.analysis,values=("DMD","Hankel DMD","EDMD","SINDy"),state="readonly",width=14).pack(side="left")
        ttk.Button(actions,text="Analyze selected",command=self.analyze).pack(side="left",padx=5)
        self.view_kind=tk.StringVar(value="time")
        ttk.Combobox(actions,textvariable=self.view_kind,values=("time","fourier","phase","three_d"),state="readonly",width=10).pack(side="left",padx=(10,0))
        ttk.Button(actions,text="View selected",command=self.view_selected).pack(side="left")
        self.refresh_archives(); self.after(250,self.poll)
    def _presets(self):
        values=self.service.available_systems()[self.system.get()]; self.preset_box["values"]=values; self.preset.set(values[0] if values else "")
    def run(self):
        try:
            cfg=RunConfiguration(self.system.get(),self.preset.get() or None,time_span=(float(self.start.get()),float(self.stop.get())),interpolation=__import__('solver').InterpolationOptions(enabled=True,time_step=float(self.step.get())),label=f"gui-{self.system.get()}")
            initial="__system_initial_state__" if cfg.initial_state is None else list(cfg.initial_state)
            job=JobManifest(label=cfg.label,systems=({"system_name":cfg.system_name,"preset":cfg.preset,"parameters":cfg.parameters,"driver":cfg.driver,"time_span":list(cfg.time_span),"output":{"delta_t":float(self.step.get())},"solver":{"options":asdict(cfg.solver),"history":asdict(cfg.history),"interpolation_method":cfg.interpolation.method},"hierarchy":{"clusters":[{"label":"gui","initial_conditions":[initial]}]}},))
            self.jobs.submit(job); self.active=job.job_id; Thread(target=lambda:self.jobs.run([job.job_id]),daemon=True).start(); self.status.set(f"Running {job.job_id[:8]}")
        except Exception as error: messagebox.showerror("Configuration error",str(error))
    def cancel(self):
        if self.active: self.status.set("Cancellation is not yet available for durable jobs; the current job will complete or record failure.")
    def poll(self):
        if self.active and self.jobs.status(self.active).get("status") in {"completed","failed","blocked"}:
            outcome=self.jobs.status(self.active); self.status.set(f"{outcome.get('status')}: {outcome.get('reason','completed')}"); self.active=None; self.refresh_archives()
        self.after(250,self.poll)
    def refresh_archives(self):
        self.archives.delete(0,"end")
        for path in sorted(VIEWER_ARCHIVE.solution_paths(),key=lambda item:item.stat().st_mtime,reverse=True): self.archives.insert("end",str(path))
    def selected_archive(self):
        choice=self.archives.curselection()
        if not choice: raise ValueError("Select a solution archive first.")
        return self.archives.get(choice[0])
    def analyze(self):
        try:
            path=self.selected_archive(); methods={"DMD":dmd_archive,"Hankel DMD":hankel_dmd_archive,"EDMD":edmd_archive,"SINDy":sindy_archive}
            output=methods[self.analysis.get()](path)
            self.status.set(f"Saved {output}")
        except Exception as error: messagebox.showerror("Postprocessing error",str(error))
    def view_selected(self):
        try:
            source=self.selected_archive(); kind=self.view_kind.get(); report=Viewer(tile_plan=TilePlan.preset("single")).add(PlotRequest(kind,ArchiveSource(source))).render()
            target=VIEWER_ARCHIVE.archive_root.parent/"viewer_exports"/(Path(source).stem+f"-{kind}.html"); Viewer().write_html(report,target,interactive=True)
            webbrowser.open(target.resolve().as_uri()); self.status.set(f"Viewer opened: {target.name}")
        except Exception as error: messagebox.showerror("Viewer error",str(error))

def launch():
    root=tk.Tk(); root.title("dynsyslab"); DynsysLabDesktop(root); root.minsize(560,420); root.mainloop()

if __name__=="__main__": launch()
