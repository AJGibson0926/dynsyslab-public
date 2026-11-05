from pathlib import Path
from ..core import load_solution, save_derived
from .algorithms import bagged_pca, havok, pca, robust_pca, snapshot_pod, svd

def _save(source,result):
    return save_derived(source,result.method,mean=result.mean,components=result.components,singular_values=result.singular_values,scores=result.scores,explained_variance_ratio=result.explained_variance_ratio)
def pca_archive(path: str|Path,rank: int|None=None):
    source=load_solution(path); return _save(source,pca(source.values,rank=rank))
def pod_archive(path: str|Path,rank: int|None=None):
    source=load_solution(path); return _save(source,snapshot_pod(source.values,rank=rank))
def bagged_pca_archive(path: str|Path,rank: int|None=None,bags: int=32,seed: int|None=None):
    source=load_solution(path); return _save(source,bagged_pca(source.values,rank=rank,bags=bags,seed=seed))
def robust_pca_archive(path: str|Path,**kwargs):
    source=load_solution(path); result=robust_pca(source.values,**kwargs)
    return save_derived(source,"robust-pca",low_rank=result["low_rank"],sparse=result["sparse"],iterations=result["iterations"],relative_residual=result["relative_residual"],lambda_value=result["lambda"])
def svd_archive(path: str|Path,*,center: bool=False,full_matrices: bool=False):
    source=load_solution(path); result=svd(source.time_series,center=center,full_matrices=full_matrices,source_metadata=source.metadata,source_archive=str(source.source))
    return save_derived(source,"svd",analysis_metadata=result.metadata,left_singular_vectors=result.left_singular_vectors,singular_values=result.singular_values,right_singular_vectors=result.right_singular_vectors,mean=result.mean)
def havok_archive(path: str|Path,**settings):
    source=load_solution(path); result=havok(source.time_series,**settings)
    return save_derived(source,"havok",analysis_metadata=result.metadata,A=result.A,B=result.B,forcing=result.forcing,coordinates=result.coordinates,coordinate_time=result.coordinate_time,left_singular_vectors=result.left_singular_vectors,singular_values=result.singular_values,right_singular_vectors=result.right_singular_vectors,observable=result.observable)
def pod_processor(rank: int|None=None):
    def processor(path: str): return pod_archive(path,rank)
    processor.__name__="snapshot_pod"; return processor
