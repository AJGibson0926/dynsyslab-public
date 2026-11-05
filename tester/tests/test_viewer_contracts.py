from __future__ import annotations
import tempfile
import unittest
import json
from pathlib import Path
from unittest.mock import patch
import numpy as np
from time_series import TimeSeries
from viewer.sources import axis_context, display_series
from viewer.core import _automatic_time_detail_inset, _display_analysis_series, _krom_colors, _krom_display_indices, PlotStyle
from viewer.routing import classify, required
from viewer.routing import matching
from viewer import (ArchiveSource, ControlErrorSource, KLQESource, KoopmanSource,
                    ObservabilitySource, PlotRequest,
                    TilePlan, Viewer, ViewerFigure, ViewerReport)
from rom import (ConstantControlMatrix, ControlApproximationErrorField,
                 ControlApproximationEvaluation,
                 save_control_approximation_error_fields,
                 save_control_approximation_evaluations)
from archives.store import serializable_system
from dynsyslib import load_system

class ViewerContractTests(unittest.TestCase):
    def test_observability_analysis_has_native_strength_direction_and_participation_panels(self):
        from rom import analyze_observability
        analysis=analyze_observability(np.diag((-.5,-1.)),np.eye(2),
            method="continuous-gramian",horizon=2.,coordinate_groups=((0,),(1,)))
        report=Viewer().add(PlotRequest("observability",ObservabilitySource(analysis))).render()
        self.assertEqual(len(report.figures),3)
        self.assertEqual(report.figures[0].metadata["observability"]["numerical_rank"],2)
        self.assertTrue(report.figures[1].metadata["directions_are_columns"])
        self.assertTrue(report.figures[1].traces[0]["heatmap"])
        self.assertEqual(required("observability"),("observability",))
        with tempfile.TemporaryDirectory() as directory:
            html=Viewer(renderer="plotly").write_html(report,Path(directory)/"observability.html").read_text(encoding="utf-8")
        self.assertIn('"type": "heatmap"',html)

    def test_klqe_result_has_native_state_measurement_and_innovation_panels(self):
        from estimation import KLQEResult
        time=np.linspace(0.,1.,5); estimate=TimeSeries(time,np.sin(time)[:,None])
        truth=TimeSeries(time,np.cos(time)[:,None]); measurement=TimeSeries(time,truth.values)
        predicted=TimeSeries(time,.9*truth.values); innovation=TimeSeries(time,.1*truth.values)
        result=KLQEResult(estimate,estimate,measurement,predicted,innovation,
                          {"estimator":{"stable":True,"observer_poles":[{"real":-1.,"imag":0.}]}},truth,
                          TimeSeries(time,truth.values-estimate.values))
        report=Viewer().add(PlotRequest("klqe",KLQESource(result))).render()
        self.assertEqual(len(report.figures),5)
        self.assertEqual([figure.metadata["role"] for figure in report.figures],
                         ["state_estimate","measurement","innovation","error_norm","observer_poles"])
        self.assertEqual(required("klqe"),("klqe",))

    def test_scalar_observable_time_series_and_phase_heatmap_are_native(self):
        system=load_system("duffing",preset="gibson_damped_stabilization_2022")
        time=np.linspace(0.,1.,5); states=np.column_stack((np.linspace(-1.,1.,5),np.zeros(5)))
        series=TimeSeries(time,states); source=ArchiveSource("unused.npz",label="duffing")
        metadata={"system":"duffing","integration_form":"dimensionless",
                  "system_parameters":serializable_system(system)}
        with patch("viewer.core._load_archive_source",return_value=[("duffing",series,metadata)]):
            time_figure=Viewer()._observable_figures(PlotRequest("observable_time",source,observable="total_energy"))[0]
            heatmap=Viewer()._observable_figures(PlotRequest("observable_heatmap",source,observable="potential_energy",domain=((-1.,1.),(-2.,2.)),domain_density=9))[0]
        self.assertEqual(time_figure.kind,"observable_time")
        self.assertEqual(time_figure.traces[0]["y"].shape,(5,))
        self.assertEqual(heatmap.kind,"observable_heatmap")
        self.assertEqual(heatmap.traces[0]["z"].shape,(9,9))
        self.assertTrue(heatmap.metadata["square"])
        with tempfile.TemporaryDirectory() as directory:
            report=ViewerReport([time_figure,heatmap],TilePlan("observables",2))
            plotly=Viewer(renderer="plotly").write_html(report,Path(directory)/"observables-plotly.html").read_text(encoding="utf-8")
            svg=Viewer(renderer="svg").write_html(report,Path(directory)/"observables-svg.html").read_text(encoding="utf-8")
        self.assertIn('"type": "heatmap"',plotly)
        self.assertIn("Potential energy",svg)

    def test_complex_scalar_observable_defaults_to_magnitude(self):
        system=load_system("duffing"); time=np.array((0.,1.)); states=np.array(((0.,0.),(.5,1.)))
        source=ArchiveSource("unused.npz"); metadata={"system":"duffing","integration_form":"dimensionless","system_parameters":serializable_system(system)}
        with patch("viewer.core._load_archive_source",return_value=[("duffing",TimeSeries(time,states),metadata)]):
            figure=Viewer()._observable_figures(PlotRequest("observable_time",source,observable="random_observation"))[0]
        self.assertEqual(figure.metadata["observable_component"],"magnitude")
        np.testing.assert_allclose(figure.traces[0]["y"],1/np.cosh(states[:,1]))

    def test_observable_plots_route_to_solution_archives(self):
        self.assertEqual(required("observable_time"),("solution",))
        self.assertEqual(required("observable_heatmap"),("solution",))
    def test_source_transforms_are_read_only_and_dimension_aware(self):
        series=TimeSeries(np.array([0.,1.]),np.array([[1.,2.],[3.,4.]]))
        time,values,context=display_series(series,{"system":"other"})
        self.assertTrue(np.array_equal(time,series.time)); self.assertTrue(np.array_equal(values,series.values)); self.assertEqual(context["time_label"],"t")
    def test_rpe_axis_context_reports_nondimensional_labels(self):
        context=axis_context({"system":"rayleigh-plesset","integration_form":"dimensionless","system_parameters":{"equilibrium_radius_m":1e-6,"ambient_pressure_pa":1e5,"preset":{"rho_kg_m3":1000.}}})
        self.assertEqual(context["time_label"],"τ")

    def test_dimensionless_derived_states_are_not_scaled_twice(self):
        source_metadata={"system":"rayleigh-plesset","integration_form":"dimensionless",
            "system_parameters":{"equilibrium_radius_m":1e-5,"ambient_pressure_pa":1e5,"preset":{"rho_kg_m3":998.}},
            "nondimensionalization":{"time_scale_s":2e-6}}
        series=TimeSeries(np.array((0.,2e-6)),np.array(((2.,0.),(1.,1.))))
        time,values,context=_display_analysis_series(series,source_metadata,{"preprocessing":{"coordinate_system":"dimensionless"}})
        np.testing.assert_allclose(time,np.array((0.,1.)))
        np.testing.assert_array_equal(values,series.values)
        self.assertEqual(context["state_labels"][0],"r = R/R₀")
    def test_routing_policy_is_independent_of_rendering(self):
        class ArchiveSource: pass
        self.assertEqual(classify(ArchiveSource()),"solution")
        self.assertEqual(required("eigenfunction"),("koopman_rom",))

    def test_routing_returns_only_compatible_attached_sources(self):
        class Entry:
            def __init__(self,data_type): self.data_type=data_type
        entries=[Entry("solution"),Entry("koopman_rom"),Entry("campaign")]
        self.assertEqual(len(matching(entries,"time")),2)
        self.assertEqual(len(matching(entries,"koopman_spectrum")),1)

    def test_unknown_source_classification_is_explicit(self):
        self.assertEqual(classify(object()),"unknown")

    def test_physical_rpe_axis_labels_are_engineering_units(self):
        context=axis_context({"system":"rayleigh-plesset","integration_form":"dimensional"})
        self.assertEqual(context["time_label"],"t (µs)")

    def test_time_plot_requires_solution_or_campaign_source(self):
        self.assertEqual(required("time"),("solution","campaign"))

    def test_time_plot_uses_distinct_default_colors_for_state_components(self):
        series=TimeSeries(np.array([0.,1.]),np.array([[1.,2.],[3.,4.]]))
        source=ArchiveSource("unused.npz",label="trajectory")
        with patch("viewer.core._load_archive_source",return_value=[("trajectory",series,{"system":"other"})]):
            figures=Viewer()._state_figures(PlotRequest("time",source))
        self.assertEqual(len(figures),2)
        self.assertEqual([len(figure.traces) for figure in figures],[1,1])
        self.assertNotEqual(figures[0].traces[0]["color"],figures[1].traces[0]["color"])

    def test_combined_time_series_remains_an_explicit_display_option(self):
        series=TimeSeries(np.array([0.,1.]),np.array([[1.,2.],[3.,4.]]))
        source=ArchiveSource("unused.npz",label="trajectory")
        with patch("viewer.core._load_archive_source",return_value=[("trajectory",series,{"system":"other"})]):
            figures=Viewer()._state_figures(PlotRequest("time",source,time_series_layout="combined"))
        self.assertEqual(len(figures),1)
        self.assertEqual(len(figures[0].traces),2)
        self.assertEqual(figures[0].metadata["time_series_layout"],"combined")

    def test_state_display_limit_preserves_sharp_control_events_without_changing_archives(self):
        time=np.linspace(0.,20.,20001); values=np.column_stack((np.sin(time),np.cos(time)))
        control_values=np.zeros((time.size,1)); control_values[1000,0]=900.
        series=TimeSeries(time,values); control=TimeSeries(time,control_values)
        source=ArchiveSource("unused.npz",label="controlled")
        metadata={"system":"other","controller":{"type":"test"},"control_history":{"present":True}}
        with patch("viewer.core._load_archive_source",return_value=[("controlled",series,metadata)]), \
             patch("viewer.core._sources.control_target_series",return_value=None), \
             patch("viewer.core._sources.control_series",return_value=control):
            figures=Viewer()._state_figures(PlotRequest("time",source,max_display_samples=250))
        self.assertTrue(all(max(len(trace["x"]) for trace in figure.traces)<=250 for figure in figures))
        control_figure=figures[-1]
        self.assertEqual(float(np.max(control_figure.traces[0]["y"])),900.)
        self.assertEqual(control_figure.metadata["display_downsampling"]["maximum_source_samples"],time.size)
        self.assertTrue(control_figure.metadata["display_downsampling"]["archive_data_unchanged"])
        self.assertEqual(series.time.size,time.size)

    def test_controlled_solution_displays_target_overlay_and_control_history(self):
        time=np.array([0.,1.,2.]); series=TimeSeries(time,np.column_stack((time,time**2)))
        target=TimeSeries(time,np.column_stack((1.+time,2.*time)))
        control=TimeSeries(time,np.column_stack((-.5*time,.25*time)))
        source=ArchiveSource("unused.npz",label="controlled")
        metadata={"system":"other","controller":{"type":"test"},"control_target_history":{"present":True},"control_history":{"present":True}}
        with patch("viewer.core._load_archive_source",return_value=[("controlled",series,metadata)]), \
             patch("viewer.core._sources.control_target_series",return_value=target), \
             patch("viewer.core._sources.control_series",return_value=control):
            time_figures=Viewer()._state_figures(PlotRequest("time",source))
            phase=Viewer()._state_figures(PlotRequest("phase",source))[0]
        self.assertEqual(len(time_figures),4)
        self.assertTrue(all(figure.metadata["control_target_overlay"] for figure in time_figures[:2]))
        self.assertTrue(all(figure.metadata["control_history"] for figure in time_figures[2:]))
        overlay=phase.traces[-1]
        self.assertEqual((overlay["color"],overlay["line_style"],overlay["role"]),("#d62728","dash","control_target"))

    def test_rpe_control_history_is_displayed_as_dimensionless_pressure(self):
        time=np.array([0.,1.]); series=TimeSeries(time,np.ones((2,2))); control=TimeSeries(time,np.ones((2,1)))
        source=ArchiveSource("unused.npz",label="controlled")
        metadata={"system":"rayleigh-plesset","controller":{"type":"test"},"control_history":{"present":True},
                  "nondimensionalization":{"pressure_scale_pa":4.}}
        with patch("viewer.core._load_archive_source",return_value=[("controlled",series,metadata)]), \
             patch("viewer.core._sources.control_target_series",return_value=None), \
             patch("viewer.core._sources.control_series",return_value=control):
            control_figure=Viewer()._state_figures(PlotRequest("time",source))[-1]
        self.assertEqual(control_figure.metadata["axes"]["y"],"dimensionless pressure control p_c* = p_c/P0")
        self.assertIn("dimensionless pressure control",control_figure.title)
        np.testing.assert_allclose(control_figure.traces[0]["y"],np.full(2,.25))

    def test_rpe_tagged_dimensionless_control_is_not_scaled_twice(self):
        time=np.array([0.,1.]); series=TimeSeries(time,np.ones((2,2))); control=TimeSeries(time,np.full((2,1),.25))
        source=ArchiveSource("unused.npz",label="controlled")
        metadata={"system":"rayleigh-plesset","controller":{"type":"test"},
                  "control_history":{"present":True,"coordinate_system":"dimensionless","units":"p_control/p0"},
                  "nondimensionalization":{"pressure_scale_pa":4.}}
        with patch("viewer.core._load_archive_source",return_value=[("controlled",series,metadata)]), \
             patch("viewer.core._sources.control_target_series",return_value=None), \
             patch("viewer.core._sources.control_series",return_value=control):
            control_figure=Viewer()._state_figures(PlotRequest("time",source))[-1]
        self.assertEqual(control_figure.metadata["axes"]["y"],"dimensionless pressure control u* = p_control/p0")
        np.testing.assert_allclose(control_figure.traces[0]["y"],np.full(2,.25))

    def test_dashed_target_style_reaches_plotly_and_svg_renderers(self):
        figure=ViewerFigure("phase","Controlled",[{"x":np.array([0.,1.]),"y":np.array([0.,1.]),"color":"#d62728","line_style":"dash"}],{"square":True})
        report=ViewerReport([figure],TilePlan.preset("single"))
        with tempfile.TemporaryDirectory() as directory:
            plotly=Viewer(renderer="plotly").write_html(report,Path(directory)/"plotly.html").read_text(encoding="utf-8")
            svg=Viewer(renderer="svg").write_html(report,Path(directory)/"svg.html").read_text(encoding="utf-8")
        self.assertIn('"dash": "dash"',plotly)
        self.assertIn("stroke-dasharray='8 6'",svg)

    def test_time_phase_layout_spans_phase_above_separate_component_panels(self):
        series=TimeSeries(np.array([0.,1.]),np.array([[1.,2.],[3.,4.]]))
        source=ArchiveSource("unused.npz",label="trajectory")
        viewer=Viewer(tile_plan=TilePlan.preset("time_phase")).add(PlotRequest("phase",source)).add(PlotRequest("time",source))
        with patch("viewer.core._load_archive_source",return_value=[("trajectory",series,{"system":"other"})]):
            report=viewer.render()
        self.assertEqual([figure.kind for figure in report.figures],["phase","time","time"])
        self.assertEqual([figure.metadata.get("grid_span",1) for figure in report.figures],[2,1,1])

    def test_time_phase_layout_spans_control_history_across_full_width(self):
        time=np.array([0.,1.]); series=TimeSeries(time,np.ones((2,2))); control=TimeSeries(time,np.ones((2,1)))
        source=ArchiveSource("unused.npz",label="controlled")
        metadata={"system":"other","control_history":{"present":True}}
        viewer=Viewer(tile_plan=TilePlan.preset("time_phase")).add(PlotRequest("time",source))
        with patch("viewer.core._load_archive_source",return_value=[("controlled",series,metadata)]), \
             patch("viewer.core._sources.control_target_series",return_value=None), \
             patch("viewer.core._sources.control_series",return_value=control):
            report=viewer.render()
        self.assertTrue(report.figures[-1].metadata["control_history"])
        self.assertEqual(report.figures[-1].metadata["grid_span"],2)

    def test_krom_dashboard_routes_to_koopman_sources(self):
        self.assertEqual(required("krom"),("koopman_rom",))

    def test_control_error_archive_routes_and_renders_native_heatmaps(self):
        x=np.linspace(-1.,1.,5); y=np.linspace(-2.,2.,7); values=np.add.outer(y*y,x*x)
        field=ControlApproximationErrorField(x,y,values,values+.1,(values+.1)/np.maximum(values,1e-3),
            "constant",{"norm":"frobenius"})
        with tempfile.TemporaryDirectory() as directory:
            archive=save_control_approximation_error_fields(Path(directory)/"errors.npz",{"constant":field},
                metadata={"axes":{"x":"x","y":"xdot"}})
            source=ControlErrorSource(archive,label="control errors")
            viewer=Viewer(renderer="plotly").add(PlotRequest("control_error",source))
            report=viewer.render()
            self.assertEqual(required("control_error"),("control_error",))
            self.assertEqual(report.figures[0].kind,"control_error")
            self.assertEqual(report.figures[0].metadata["axes"],{"x":"x","y":"xdot"})
            plotly=viewer.write_html(report,Path(directory)/"plotly.html").read_text(encoding="utf-8")
            svg=Viewer(renderer="svg").write_html(report,Path(directory)/"svg.html").read_text(encoding="utf-8")
        self.assertIn('"type": "heatmap"',plotly)
        self.assertNotIn('"scaleanchor": "x"',plotly)
        self.assertIn('"square": true',plotly)
        self.assertIn("width:min(100%,500px);aspect-ratio:1/1",plotly)
        self.assertIn("relative Frobenius error over phase space",svg)

    def test_control_evaluation_archive_renders_each_retained_fit(self):
        x=np.linspace(-1.,1.,4); y=np.linspace(-1.,1.,3); values=np.ones((3,4))
        field=ControlApproximationErrorField(x,y,values,values*.2,values*.2,
            "constant",{"axes":{"x":"x","y":"v"}})
        evaluation=ControlApproximationEvaluation("fit-1","controlled set A","constant",
            ConstantControlMatrix(np.asarray(((1.+0j,),))),
            {"fingerprint_sha256":"abc"},field)
        with tempfile.TemporaryDirectory() as directory:
            archive=save_control_approximation_evaluations(Path(directory)/"evaluations.npz",(evaluation,))
            report=Viewer().add(PlotRequest("control_error",ControlErrorSource(archive))).render()
        self.assertEqual(len(report.figures),1)
        self.assertEqual(report.figures[0].title,"controlled set A: relative Frobenius error")

    def test_krom_dashboard_selection_never_changes_model_retention(self):
        retained=(0,2,4,6,8,10,12)
        self.assertEqual(_krom_display_indices(PlotRequest("krom"),retained),retained[:5])
        self.assertEqual(_krom_display_indices(PlotRequest("krom",krom_mode_limit=None),retained),retained)
        self.assertEqual(_krom_display_indices(PlotRequest("krom",krom_mode_indices=(8,2)),retained),(8,2))
        with self.assertRaisesRegex(ValueError,"not retained"):
            _krom_display_indices(PlotRequest("krom",krom_mode_indices=(1,)),retained)
        colors=_krom_colors(PlotStyle(),12)
        self.assertEqual(len(colors),len(set(colors)))

    def test_archived_krom_reconstructs_only_default_display_columns(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/"sindy.npz"; states=np.column_stack((np.arange(5.),np.arange(5.)**2))
            retained=(0,2,4,6,8,10,12)
            metadata={"source_archive":"source.npz","analysis_metadata":{
                "krom":{"retained_indices":list(retained)},"preprocessing":{"coordinate_system":"dimensional"}}}
            np.savez(path,training_states=states,metadata_json=json.dumps(metadata))
            displayed=retained[:5]; selected_values=np.ones((5,5),dtype=complex)
            with patch("viewer.core._sources.derived_array",return_value=selected_values) as loader, \
                 patch.object(Viewer,"_archived_eigenvalues",return_value=np.linspace(.8,1.,13).astype(complex)), \
                 patch("viewer.core._solution_record",return_value=(TimeSeries(np.arange(5.),states),{})):
                figures=Viewer()._krom_figures(PlotRequest("krom",KoopmanSource(path)))
            loader.assert_called_once_with(path.resolve(),"eigenfunction_values",columns=displayed)
            self.assertEqual(len(figures),26)
            spectrum=figures[0]
            self.assertEqual(spectrum.metadata["retained_mode_count"],7)
            self.assertEqual(spectrum.metadata["displayed_mode_indices"],list(displayed))
            self.assertEqual(spectrum.metadata["omitted_dashboard_mode_count"],2)
            self.assertFalse(spectrum.metadata["display_selection_is_model_truncation"])
            self.assertEqual(spectrum.traces[2]["x"].size,7)

    def test_viewer_attaches_named_koopman_source_and_auto_resolves_it(self):
        viewer=Viewer(); source=KoopmanSource(object(),label="spectrum")
        viewer.attach_source(source)
        resolved=viewer._resolve_request_source(PlotRequest(kind="koopman_spectrum"))
        self.assertIs(resolved.source,source)

    def test_viewer_requires_explicit_source_when_multiple_match(self):
        viewer=Viewer(); viewer.attach_source(KoopmanSource(object(),label="one")); viewer.attach_source(KoopmanSource(object(),label="two"))
        with self.assertRaises(ValueError): viewer._resolve_request_source(PlotRequest(kind="koopman_spectrum"))

    def test_viewer_rejects_duplicate_source_names(self):
        viewer=Viewer(); viewer.attach_source(KoopmanSource(object()),name="same")
        with self.assertRaises(ValueError): viewer.attach_source(KoopmanSource(object()),name="same")

    def test_tile_plan_presets_have_positive_columns(self):
        for name in ("single","grid_2x2","time_spectrum","dashboard_2x2"):
            self.assertGreater(TilePlan.preset(name).columns,0)

    def test_viewer_configure_rejects_unknown_renderer(self):
        with self.assertRaises(ValueError): Viewer().configure(renderer="unknown")

    def test_plotly_line_trace_omits_null_marker_and_declares_utf8(self):
        figure=ViewerFigure("phase","RPE ensemble",[
            {"label":"trajectory","x":np.array([1.,2.]),"y":np.array([0.,1.]),"color":"#123456","width":1.5}
        ],{"axes":{"x":"R (µm)","y":"dR/dt (m/s)"}})
        report=ViewerReport([figure],TilePlan.preset("single"))
        with tempfile.TemporaryDirectory() as directory:
            path=Viewer(renderer="plotly").write_html(report,Path(directory)/"phase.html")
            html=path.read_text(encoding="utf-8")
        self.assertIn("<meta charset='utf-8'>",html)
        self.assertIn("R (µm)",html)
        self.assertNotIn('"marker": null',html)
        self.assertIn('"orientation": "h"',html)
        self.assertIn("section{width:100%;min-width:0",html)
        self.assertNotIn("overflow:hidden",html)

    def test_plotly_payload_deduplicates_repeated_large_coordinate_arrays(self):
        shared=np.arange(1000.,1100.)
        figures=[ViewerFigure("time",f"panel-{index}",[{"x":shared,"y":np.full(100,index,dtype=float)}])
                 for index in range(2)]
        with tempfile.TemporaryDirectory() as directory:
            path=Viewer(renderer="plotly").write_html(ViewerReport(figures,TilePlan.preset("vertical")),
                                                       Path(directory)/"shared.html")
            html=path.read_text(encoding="utf-8")
        self.assertIn("dynsyslab.viewer-plotly/v2",html)
        self.assertEqual(html.count("1000.0, 1001.0, 1002.0"),1)
        self.assertGreaterEqual(html.count('"__dynsyslab_array__": 0'),2)

    def test_time_detail_inset_detects_a_large_transient_but_not_an_ordinary_signal(self):
        time=np.linspace(0.,10.,1001); detail=12.*np.sin(2.*np.pi*time)
        transient=detail.copy(); transient[0]=-63000.
        inset=_automatic_time_detail_inset([{"x":time,"y":transient}])
        self.assertIsNotNone(inset)
        self.assertEqual(inset["reason"],"transient-scale-disparity")
        self.assertGreater(inset["scale_ratio"],100.)
        self.assertGreater(inset["ylim"][0],-100.)
        self.assertIsNone(_automatic_time_detail_inset([{"x":time,"y":detail}]))

    def test_time_detail_inset_reaches_plotly_and_svg_renderers(self):
        time=np.linspace(0.,10.,101); values=np.sin(time); values[0]=-1000.
        inset=_automatic_time_detail_inset([{"x":time,"y":values,"color":"#123456"}])
        figure=ViewerFigure("time","Transient",[{"x":time,"y":values,"color":"#123456"}],
                            {"axes":{"x":"t","y":"u"},"detail_inset":inset})
        report=ViewerReport([figure],TilePlan.preset("single"))
        with tempfile.TemporaryDirectory() as directory:
            plotly=Viewer(renderer="plotly").write_html(report,Path(directory)/"plotly.html").read_text(encoding="utf-8")
            svg=Viewer(renderer="svg").write_html(report,Path(directory)/"svg.html").read_text(encoding="utf-8")
        self.assertIn('"xaxis2"',plotly)
        self.assertIn('"yaxis": "y2"',plotly)
        self.assertIn("data-viewer-detail-inset='true'",svg)

    def test_plotly_koopman_spectrum_locks_complex_plane_aspect_ratio(self):
        figure=ViewerFigure("koopman_spectrum","Spectrum",[
            {"label":"eigenvalues","x":np.array([0.,1.]),"y":np.array([1.,0.]),"marker":True}
        ],{"axes":{"x":"Re(λ)","y":"Im(λ)"},"title_color":"#a020f0","height":420,"grid_span":2})
        report=ViewerReport([figure],TilePlan.preset("single"))
        with tempfile.TemporaryDirectory() as directory:
            path=Viewer(renderer="plotly").write_html(report,Path(directory)/"spectrum.html")
            html=path.read_text(encoding="utf-8")
        self.assertIn('"scaleanchor": "x"',html)
        self.assertIn('"scaleratio": 1.0',html)
        self.assertIn('"color": "#a020f0"',html)
        self.assertIn("height:420px;grid-column:span 2",html)

    def test_archived_spectrum_can_use_an_explicit_eigenvalue_source(self):
        source=KoopmanSource("sindy.npz",eigenvalue_source="hdmd.npz")
        with patch.object(Viewer,"_archived_eigenvalues",return_value=np.array((1.+0j,.9+.1j))) as loader:
            figure=Viewer()._koopman_figure(PlotRequest("koopman_spectrum",source))
        loader.assert_called_once_with(Path("hdmd.npz"))
        self.assertEqual(figure.traces[-1]["x"].size,2)
        self.assertEqual(figure.metadata["domain"],"discrete-time")
        self.assertEqual(figure.metadata["axes"],{"x":"Re(μ) (unitless)","y":"Im(μ) (unitless)"})

    def test_archived_spectrum_follows_exact_lineage_path(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory); upstream=root/"hdmd.npz"; derived=root/"sindy.npz"
            np.savez(upstream,eigenvalues=np.array((1.+0j,.9+.1j)),metadata_json=json.dumps({}))
            np.savez(derived,metadata_json=json.dumps({"analysis_metadata":{"eigenvalue_source_archive":str(upstream)}}))
            np.testing.assert_array_equal(Viewer._archived_eigenvalues(derived),np.array((1.+0j,.9+.1j)))

    def test_plotly_phase_portrait_uses_a_square_canvas(self):
        figure=ViewerFigure("phase","Phase portrait",[
            {"label":"orbit","x":np.array([-1.,1.]),"y":np.array([-1.,1.])}
        ],{"axes":{"x":"x","y":"v"},"height":520,"square":True})
        report=ViewerReport([figure],TilePlan.preset("single"))
        with tempfile.TemporaryDirectory() as directory:
            path=Viewer(renderer="plotly").write_html(report,Path(directory)/"phase.html")
            html=path.read_text(encoding="utf-8")
        self.assertNotIn('"scaleanchor": "x"',html)
        self.assertIn('"square": true',html)
        self.assertIn("width:min(100%,520px);aspect-ratio:1/1",html)

    def test_canonical_krom_dashboard_has_five_panels_per_retained_mode(self):
        eigenvalues=np.array([.8+.1j,.8-.1j,.9+.2j,.9-.2j])
        states=np.column_stack((np.linspace(1.,2.,5),np.linspace(-1.,1.,5)))
        values=np.column_stack([np.exp(1j*np.linspace(0.,1.,5)*(index+1)) for index in range(4)])
        colors=("#112233","#445566")
        figures=Viewer._krom_dashboard_figures(eigenvalues,(0,2),states,np.linspace(0.,4.,5),values,
            {"time_label":"t","state_labels":["x","v"]},colors)
        self.assertEqual(len(figures),11)
        self.assertEqual(figures[0].kind,"koopman_spectrum")
        self.assertEqual([trace["color"] for trace in figures[0].traces[-2:]],list(colors))
        self.assertEqual(figures[1].metadata["title_color"],colors[0])
        self.assertEqual([figure.title for figure in figures[2:6]],
                         ["ψ0 · Phase","ψ0 · Magnitude","ψ0 · Real part","ψ0 · Imaginary part"])

    def test_canonical_krom_dashboard_accepts_preselected_value_columns(self):
        eigenvalues=np.array([.8+.1j,.8-.1j,.9+.2j,.9-.2j])
        states=np.column_stack((np.linspace(1.,2.,5),np.linspace(-1.,1.,5)))
        full=np.column_stack([np.exp(1j*np.linspace(0.,1.,5)*(index+1)) for index in range(4)])
        figures=Viewer._krom_dashboard_figures(eigenvalues,(0,2),states,np.linspace(0.,4.,5),full[:,[0,2]],
            {"time_label":"t","state_labels":["x","v"]},("#112233","#445566"))
        np.testing.assert_allclose(figures[2].traces[0]["y"],np.angle(full[:,0]))
        np.testing.assert_allclose(figures[7].traces[0]["y"],np.angle(full[:,2]))

    def test_krom_spectrum_distinguishes_retained_and_displayed_modes(self):
        eigenvalues=np.array([.8+.1j,.8-.1j,.9+.2j,.9-.2j])
        states=np.column_stack((np.linspace(1.,2.,5),np.linspace(-1.,1.,5)))
        values=np.column_stack((np.ones(5),1j*np.ones(5)))
        figures=Viewer._krom_dashboard_figures(eigenvalues,(0,2),states,np.arange(5.),values,
            {"time_label":"t","state_labels":["x","v"]},("#112233","#445566"),retained_indices=(0,1,2,3))
        spectrum=figures[0]
        self.assertEqual(spectrum.metadata["full_spectrum_count"],4)
        self.assertEqual(spectrum.metadata["retained_mode_count"],4)
        self.assertEqual(spectrum.metadata["displayed_mode_indices"],[0,2])
        self.assertEqual(spectrum.traces[2]["label"],"KROM retained spectrum")
