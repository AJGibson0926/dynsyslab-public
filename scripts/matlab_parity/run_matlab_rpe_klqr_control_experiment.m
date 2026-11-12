function [outputFile, summaryFile, figureFile] = run_matlab_rpe_klqr_control_experiment(options)
%RUN_MATLAB_RPE_KLQR_CONTROL_EXPERIMENT Reproduce one MATLAB KLQR pathway.
%
% This is a diagnostic wrapper around the existing MATLAB DynSysLab suite.
% It does not edit the suite and is not a specification for the Python/C++
% implementation. It reproduces the RPE pathway used by control_dynamics.m:
%
%   free RPE -> radius H-DMD -> SINDy eigenfunctions
%   -> first complex eigenfunction/conjugate pair -> KLQR
%   -> equilibrium-start closed-loop ode15s solve
%
% The output records both uniformly sampled results and the raw sequence of
% adaptive-solver controller calls, including the recomputed gains and the
% state-dependent lifted input matrix B_psi.

arguments
    options.SuiteRoot (1,1) string = ...
        "~/Documents/MATLAB/dynsyslab (5-15-25) (cleanup begins)"
    options.OutputFile (1,1) string = ""
    options.Seed (1,1) double = 2141442
    options.TrainingInitialState (1,2) double = [2.3333865986017348, 0]
    options.TrainingHorizonTau (1,1) double {mustBePositive} = 25
    options.TrainingDTau (1,1) double {mustBePositive} = 0.001
    options.HankelRows (1,1) double {mustBeInteger,mustBePositive} = 20000
    options.HankelColumns (1,1) double {mustBeInteger,mustBePositive} = 2000
    options.ControlInitialState (1,2) double = [1, 0]
    options.ControlHorizonTau (1,1) double {mustBePositive} = 10
    options.ControlDTau (1,1) double {mustBePositive} = 0.001
    options.TargetOffset (1,1) double = 4
    options.TargetAmplitudes (1,:) double = [1.5, 1]
    options.TargetAngularFrequencies (1,:) double = [4, 3]
    options.TargetPhases (1,:) double = [0, 0]
    options.RelTol (1,1) double {mustBePositive} = 1e-9
    options.AbsTol (1,1) double {mustBePositive} = 1e-9
    options.SaveRawControlCalls (1,1) logical = true
    options.MakeFigure (1,1) logical = true
end

validate_target_options(options);
suiteRoot = char(options.SuiteRoot);
assert(isfolder(suiteRoot), "KLQRExperiment:MissingSuite", ...
    "SuiteRoot does not exist: %s", suiteRoot);

if options.OutputFile == ""
    stamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
    outputFile = fullfile(pwd, "matlab_rpe_klqr_control_" + stamp + ".mat");
else
    outputFile = char(options.OutputFile);
end
outputFile = char(java.io.File(outputFile).getCanonicalPath());
[outputDir, outputStem] = fileparts(outputFile);
if isempty(outputDir), outputDir = pwd; end
if ~isfolder(outputDir), mkdir(outputDir); end
summaryFile = char(fullfile(outputDir, outputStem + ".json"));
figureFile = char(fullfile(outputDir, outputStem + "_results.png"));

assert(~is_path_inside(outputFile, suiteRoot), "KLQRExperiment:UnsafeOutput", ...
    "OutputFile must be outside SuiteRoot.");

originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(genpath(suiteRoot));
assert_required_symbols();

rngBefore = rng;
rng(options.Seed, "twister");
totalClock = tic;
timings = struct("training_trajectory_seconds",NaN, "hdmd_seconds",NaN, ...
    "sindy_seconds",NaN, "control_seconds",NaN, ...
    "diagnostic_reconstruction_seconds",NaN, "total_seconds",NaN);

%% Training trajectory and KROM pathway
trainingSystem = make_rpe_system(options.TrainingInitialState, ...
    options.TrainingHorizonTau, options.TrainingDTau, options.RelTol, options.AbsTol);
trainingBatch = dynbatch(systems={trainingSystem});
clock = tic;
trainingBatch.compute("soft");
timings.training_trajectory_seconds = toc(clock);

assert(numel(trainingSystem.sols.pool{1}.t) >= ...
    options.HankelRows + options.HankelColumns, ...
    "KLQRExperiment:ShortTrainingTrajectory", ...
    "Training trajectory is too short for the requested Hankel window.");

koop = Koopman;
radiusObservable = trainingSystem.dyn.obsv(1);
clock = tic;
romHDMD = koop.process(data=trainingSystem, ...
    eigenvalue_method='HDMD', eigenfunction_method='HDMD', ...
    observable=radiusObservable, rows=options.HankelRows, ...
    cols=options.HankelColumns, stack_method='horizontal', sortmode='Energy');
timings.hdmd_seconds = toc(clock);

clock = tic;
romSINDy = koop.process(data=trainingSystem, ...
    eigenvalue_method='provided', spectrum=romHDMD.spec, ...
    eigenfunction_method='SINDy', data_type='trajectory');
timings.sindy_seconds = toc(clock);

% This exact selection appears in generate_training.m and supplies the
% plant used by control_dynamics.m.
controlKROM = romSINDy.extract(rank=[0,1], method="real");
assert(controlKROM.rank == 1 && controlKROM.dim == 2, ...
    "KLQRExperiment:UnexpectedKROM", ...
    "Expected one complex eigenfunction represented as a 2-state real pair.");

%% Closed-loop MATLAB suite solve
controlSystem = make_rpe_system(options.ControlInitialState, ...
    options.ControlHorizonTau, options.ControlDTau, options.RelTol, options.AbsTol);
ctrl = controlSystem.CTRL;
ctrl.train(ctrl='KLQR', plant=controlKROM);

offset = options.TargetOffset;
amplitudes = options.TargetAmplitudes;
frequencies = options.TargetAngularFrequencies;
phases = options.TargetPhases;
ctrl.ref = @(x,t) target_state(t, offset, amplitudes, frequencies, phases);
ctrl.interp = 1;
ctrl.time = [0, options.ControlHorizonTau];

controlBatch = dynbatch(systems={controlSystem});
clock = tic;
controlBatch.compute("soft");
timings.control_seconds = toc(clock);

%% Extract uniform and raw controller evidence
uniform = extract_uniform_results(controlSystem, options);
clock = tic;
raw = extract_raw_control_calls(ctrl, controlKROM, controlSystem.dyn, options);
timings.diagnostic_reconstruction_seconds = toc(clock);
metrics = tracking_metrics(uniform);
model = extract_model_evidence(romHDMD, controlKROM, ctrl, options);

provenance = struct;
provenance.created_utc = char(datetime("now", "TimeZone","UTC", ...
    "Format","yyyy-MM-dd'T'HH:mm:ss.SSSXXX"));
provenance.wrapper_file = mfilename("fullpath");
provenance.suite_root = suiteRoot;
provenance.matlab_version = version;
provenance.matlab_release = version('-release');
provenance.computer = computer;
provenance.seed = options.Seed;
provenance.rng_before = rngBefore;
provenance.rng_after = rng;
provenance.resolved_files = struct( ...
    "control_dynamics", which("control_dynamics"), ...
    "Bubble", which("Bubble"), "CTRL", which("CTRL"), ...
    "plant", which("plant"), "Koopman", which("Koopman"), ...
    "DMD", which("DMD"), "SINDy", which("SINDy"));

experiment = struct;
experiment.schema = "dynsyslab.matlab-rpe-klqr-control-experiment";
experiment.schema_version = 1;
experiment.scope = [ ...
    "Diagnostic reproduction of one MATLAB Koopman-control pathway."; ...
    "Not a Python/C++ feature specification and not authority to copy MATLAB code." ...
    ];
experiment.provenance = provenance;
experiment.physical_parameters = scalar_table_to_struct(controlSystem.dyn.prms);
experiment.dimensionless_groups = scalar_table_to_struct(controlSystem.dyn.dimgrps);
experiment.training = struct( ...
    "initial_state", options.TrainingInitialState, ...
    "horizon_tau", options.TrainingHorizonTau, ...
    "dtau", options.TrainingDTau, ...
    "sample_count", numel(trainingSystem.sols.pool{1}.t));
experiment.model = model;
experiment.target = struct( ...
    "offset",offset, "amplitudes",amplitudes, ...
    "angular_frequencies",frequencies, "phases",phases, ...
    "radius_formula","offset+sum(amplitude_i*sin(omega_i*tau+phase_i))", ...
    "velocity_formula","sum(amplitude_i*omega_i*cos(omega_i*tau+phase_i))");
experiment.control = struct( ...
    "initial_state",options.ControlInitialState, ...
    "horizon_tau",options.ControlHorizonTau, ...
    "dtau",options.ControlDTau, ...
    "solver","ode15s through dynbatch.compute(""soft"")", ...
    "rel_tol",options.RelTol, "abs_tol",options.AbsTol, ...
    "input_units","dimensionless ambient-pressure perturbation p_control/p0", ...
    "physical_input_formula","p_control_Pa = dimensionless_input * p0");
experiment.uniform = uniform;
if options.SaveRawControlCalls, experiment.raw_control_calls = raw; end
experiment.metrics = metrics;
experiment.timings = timings;

timings.total_seconds = toc(totalClock);
experiment.timings = timings;
save(outputFile, "experiment", "-v7.3");
write_json(summaryFile, make_summary(experiment, outputFile, figureFile));

if options.MakeFigure
    make_result_figure(experiment, figureFile);
else
    figureFile = "";
end

fprintf("\nMATLAB KLQR experiment complete.\n");
fprintf("  MAT:    %s\n  JSON:   %s\n", outputFile, summaryFile);
if options.MakeFigure, fprintf("  Figure: %s\n", figureFile); end
fprintf("  Radius RMSE: %.8g (second half %.8g)\n", ...
    metrics.radius_rmse, metrics.radius_rmse_second_half);
fprintf("  Velocity RMSE: %.8g (second half %.8g)\n", ...
    metrics.velocity_rmse, metrics.velocity_rmse_second_half);
end


function validate_target_options(options)
n = numel(options.TargetAmplitudes);
assert(numel(options.TargetAngularFrequencies)==n && numel(options.TargetPhases)==n, ...
    "KLQRExperiment:TargetShape", ...
    "Target amplitudes, angular frequencies, and phases must have equal lengths.");
assert(options.TrainingDTau == options.ControlDTau, ...
    "KLQRExperiment:DTauMismatch", ...
    "The MATLAB spectrum conversion assumes the training and control dtau agree.");
end


function assert_required_symbols()
symbols = ["Bubble","dynsys","dynbatch","Koopman","CTRL","plant","lqr"];
for symbol = symbols
    assert(~isempty(which(symbol)), "KLQRExperiment:MissingSymbol", ...
        "Required MATLAB symbol '%s' was not found.", symbol);
end
end


function sys = make_rpe_system(initialState, horizon, dtau, relTol, absTol)
sys = dynsys(dynamics=Bubble, trajectories=1);
sys.dyn.prms.CL = 0;
sys.dyn.prms.CG = 0;
sys.load_inits(initialState);
sys.load_time({[0,horizon; 0,horizon],dtau});
sys.options = odeset(sys.options,"RelTol",relTol,"AbsTol",absTol);
end


function value = target_state(t, offset, amplitudes, frequencies, phases)
t = reshape(t,1,[]);
angle = frequencies(:)*t + phases(:);
radius = offset + sum(amplitudes(:).*sin(angle),1);
velocity = sum((amplitudes(:).*frequencies(:)).*cos(angle),1);
value = [radius; velocity];
end


function uniform = extract_uniform_results(sys, options)
sol = sys.sols.pool{1};
tau = double(sol.t(:));
states = double(sol.ts.');
target = target_state(tau.', options.TargetOffset, options.TargetAmplitudes, ...
    options.TargetAngularFrequencies, options.TargetPhases).';

inst = sys.clst{1}.inst{1};
ctrlTime = double(inst.ctrl_hist.t(:));
ctrlValues = double(inst.ctrl_hist.ts.');

uniform = struct;
uniform.tau = tau;
uniform.states = states;
uniform.target_states = target;
uniform.state_error = states-target;
uniform.control_tau = ctrlTime;
uniform.control_history = ctrlValues;
uniform.control_history_layout = [ ...
    "column 1: dimensionless controller output u=p_control/p0"; ...
    "column 2: filtered/plant-applied scalar Bu recorded by MATLAB CTRL" ...
    ];
uniform.control_pressure_Pa = ctrlValues(:,1) * sys.dyn.prms.p0;
end


function raw = extract_raw_control_calls(ctrl, krom, dyn, options)
N = numel(ctrl.tape.t);
raw = struct("available",false,"call_count",N);
if ~options.SaveRawControlCalls || N==0, return; end

tau = zeros(N,1);
state = zeros(N,2);
input = zeros(N,1);
applied = zeros(N,1);
liftedError = zeros(N,krom.dim);
cost = zeros(N,1);
gain = zeros(1,krom.dim,N);
liftedB = zeros(krom.dim,1,N);
closedLoopPoles = zeros(krom.dim,N);
synthetic = false(N,1);
A = ctrl.plant.dyn.A;
o = @(r) r(krom.dims);

for n = 1:N
    datum = ctrl.tape.data{n};
    tau(n) = ctrl.tape.t{n};
    input(n) = datum{1};
    applied(n) = datum{2};
    liftedError(n,:) = reshape(datum{3},1,[]);
    cost(n) = datum{4};
    state(n,:) = reshape(datum{5},1,[]);
    K = reshape(datum{9},1,[]);
    gain(1,:,n) = K;
    synthetic(n) = all(state(n,:)==0) && all(K==0);
    if ~synthetic(n)
        Bpsi = ctrl.plant.dyn.B(state(n,:).',tau(n),o);
        liftedB(:,:,n) = Bpsi;
        closedLoopPoles(:,n) = eig(A-Bpsi*K);
    else
        liftedB(:,:,n) = NaN;
        closedLoopPoles(:,n) = NaN;
    end
end

raw.available = true;
raw.tau = tau;
raw.state = state;
raw.dimensionless_input = input;
raw.control_pressure_Pa = input*dyn.prms.p0;
raw.applied_scalar_record = applied;
raw.lifted_error = liftedError;
raw.instantaneous_cost = cost;
raw.gain_K = gain;
raw.lifted_input_Bpsi = liftedB;
raw.closed_loop_poles = closedLoopPoles;
raw.synthetic_initial_record = synthetic;
raw.note = "Calls include accepted/rejected/internal ode15s evaluations and may repeat or move backward in time.";
end


function metrics = tracking_metrics(uniform)
e = uniform.state_error;
n = size(e,1);
second = max(1,floor(n/2)):n;
metrics = struct;
metrics.radius_rmse = sqrt(mean(e(:,1).^2));
metrics.velocity_rmse = sqrt(mean(e(:,2).^2));
metrics.radius_mae = mean(abs(e(:,1)));
metrics.velocity_mae = mean(abs(e(:,2)));
metrics.radius_max_abs_error = max(abs(e(:,1)));
metrics.velocity_max_abs_error = max(abs(e(:,2)));
metrics.radius_rmse_second_half = sqrt(mean(e(second,1).^2));
metrics.velocity_rmse_second_half = sqrt(mean(e(second,2).^2));
metrics.radius_target_span = max(uniform.target_states(:,1))-min(uniform.target_states(:,1));
metrics.radius_nrmse_by_target_span = metrics.radius_rmse/metrics.radius_target_span;
metrics.radius_nrmse_second_half_by_target_span = ...
    metrics.radius_rmse_second_half/metrics.radius_target_span;
metrics.control_rms_dimensionless = sqrt(mean(uniform.control_history(:,1).^2));
metrics.control_peak_abs_dimensionless = max(abs(uniform.control_history(:,1)));
end


function model = extract_model_evidence(hdmd, krom, ctrl, options)
psi = krom.psi(1);
model = struct;
model.hdmd_settings = struct("rows",options.HankelRows, ...
    "columns",options.HankelColumns,"observable","dimensionless radius", ...
    "stack_method","horizontal","sort_mode","Energy", ...
    "svd_absolute_tolerance",1e-10);
model.hdmd_retained_dimension = numel(hdmd.spec.get_spec("Discrete","Full"));
model.sindy_settings = struct("polynomial_order",psi.SINDy.polyorder, ...
    "coefficient_threshold",psi.SINDy.tol,"uses_sines",logical(psi.SINDy.usesine), ...
    "includes_constant",false,"normalization","none");
model.selection = "romSINDy.extract(rank=[0,1],method='real'): first complex eigenfunction plus implicit conjugate dual";
model.selected_discrete_eigenvalue = psi.lambda.d;
model.selected_continuous_eigenvalue = psi.lambda.c;
model.selected_dual_discrete_eigenvalue = psi.dual.lambda.d;
model.selected_dual_continuous_eigenvalue = psi.dual.lambda.c;
model.selected_sindy_coefficients = psi.psi(:);
model.selected_dual_sindy_coefficients = psi.dual.psi(:);
model.A_continuous_real_pair = ctrl.plant.dyn.A;
model.Q = ctrl.Q;
model.R = ctrl.R;
model.B_physical_state_formula = "B(x)=[0; -Eu/radius]";
model.Bpsi_formula = "Bpsi(x)=gradient(Psi(x))*B(x), recomputed inside every controller call";
model.control_law = "K(x)=lqr(A,Bpsi(x),Q,R); u=-K(x)*(Psi(x)-Psi(reference))";
end


function summary = make_summary(experiment, outputFile, figureFile)
summary = struct;
summary.schema = experiment.schema;
summary.schema_version = experiment.schema_version;
summary.mat_file = outputFile;
summary.figure_file = figureFile;
summary.scope = experiment.scope;
summary.provenance = rmfield(experiment.provenance,["rng_before","rng_after"]);
summary.physical_parameters = experiment.physical_parameters;
summary.dimensionless_groups = experiment.dimensionless_groups;
summary.training = experiment.training;
summary.model = rmfield(experiment.model, ...
    ["selected_sindy_coefficients","selected_dual_sindy_coefficients"]);
summary.model.selected_discrete_eigenvalue = complex_pair( ...
    summary.model.selected_discrete_eigenvalue);
summary.model.selected_continuous_eigenvalue = complex_pair( ...
    summary.model.selected_continuous_eigenvalue);
summary.model.selected_dual_discrete_eigenvalue = complex_pair( ...
    summary.model.selected_dual_discrete_eigenvalue);
summary.model.selected_dual_continuous_eigenvalue = complex_pair( ...
    summary.model.selected_dual_continuous_eigenvalue);
summary.target = experiment.target;
summary.control = experiment.control;
summary.metrics = experiment.metrics;
summary.timings = experiment.timings;
if isfield(experiment,"raw_control_calls")
    summary.raw_control_call_count = experiment.raw_control_calls.call_count;
end
end


function value = complex_pair(z)
value = struct("real",real(z),"imaginary",imag(z));
end


function make_result_figure(experiment, filename)
u = experiment.uniform;
fig = figure("Visible","off","Color","white","Position",[100,100,1200,950]);
layout = tiledlayout(fig,3,2,"TileSpacing","compact","Padding","compact");
title(layout,"MATLAB RPE KLQR diagnostic reproduction");

nexttile(layout,[1,2]);
plot(u.states(:,1),u.states(:,2),"b-","LineWidth",1.2); hold on;
plot(u.target_states(:,1),u.target_states(:,2),"r--","LineWidth",1.1);
xlabel("Dimensionless radius"); ylabel("Dimensionless radial velocity");
legend("Controlled solution","Target","Location","best"); grid on; axis square;

nexttile(layout);
plot(u.tau,u.states(:,1),"b-","LineWidth",1.1); hold on;
plot(u.tau,u.target_states(:,1),"r--","LineWidth",1.1);
xlabel("\tau"); ylabel("Dimensionless radius"); grid on;

nexttile(layout);
plot(u.tau,u.states(:,2),"Color",[0.1,0.55,0.15],"LineWidth",1.1); hold on;
plot(u.tau,u.target_states(:,2),"r--","LineWidth",1.1);
xlabel("\tau"); ylabel("Dimensionless radial velocity"); grid on;

nexttile(layout,[1,2]);
plot(u.control_tau,u.control_history(:,1),"k-","LineWidth",1);
xlabel("\tau"); ylabel("u/p_0 (dimensionless)"); grid on;

exportgraphics(fig,filename,"Resolution",160);
close(fig);
end


function value = scalar_table_to_struct(tbl)
value = table2struct(tbl,"ToScalar",true);
end


function write_json(filename,value)
encoded = jsonencode(value,"PrettyPrint",true);
fid = fopen(filename,"w");
assert(fid>=0,"KLQRExperiment:CannotWriteJSON","Cannot write %s",filename);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid,"%s\n",encoded);
end


function tf = is_path_inside(candidate,root)
candidate = lower(char(java.io.File(candidate).getCanonicalPath()));
root = lower(char(java.io.File(root).getCanonicalPath()));
if ~endsWith(root,filesep), root = [root filesep]; end
tf = startsWith(candidate,root);
end
