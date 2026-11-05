"""Pure Viewer source classification and automatic request routing."""
from __future__ import annotations
from typing import Any

def classify(source:Any)->str:
    return {"ArchiveSource":"solution","KoopmanSource":"koopman_rom","PODSource":"pod_rom","ControlErrorSource":"control_error","ObservabilitySource":"observability","KLQESource":"klqe","BifurcationSource":"bifurcation_campaign","CampaignSource":"campaign"}.get(type(source).__name__,"unknown")
def required(kind:str)->tuple[str,...]:
    return {"time":("solution","campaign"),"fourier":("solution","campaign"),"phase":("solution","campaign"),"three_d":("solution","campaign"),
            "observable_time":("solution",),"observable_heatmap":("solution",),
            "koopman_spectrum":("koopman_rom",),"eigenfunction":("koopman_rom",),"krom":("koopman_rom",),"pod_mode":("pod_rom",),"control_error":("control_error",),"observability":("observability",),"klqe":("klqe",),"bifurcation":("bifurcation_campaign","campaign")}.get(kind,())
def matching(entries:list[Any],kind:str)->list[Any]: return [entry for entry in entries if getattr(entry,"data_type",None) in required(kind)]
