function [outputFile, summaryFile] = export_dynsyslab_rpe_hdmd_sindy(options)
%EXPORT_DYNSYSLAB_RPE_HDMD_SINDY Export a MATLAB DynSysLab parity bundle.
%
% This wrapper is intentionally external to the MATLAB DynSysLab suite. It
% calls the suite through its existing public interfaces and saves only
% plain numeric/struct data. It never edits source files in SuiteRoot.
%
% The default run reproduces the canonical comparison chain:
%   dimensionless RPE trajectory (ode15s/"soft") -> H-DMD -> SINDy
%
% The MAT file is the authoritative bundle. The JSON sidecar is a compact,
% human-readable manifest and deliberately omits large/complex arrays.

arguments
    options.SuiteRoot (1,1) string = ...
        fullfile(userpath, "dynsyslab (5-15-25) (cleanup begins)")
    options.OutputFile (1,1) string = ""
    options.ExistingSystem = []
    options.Seed (1,1) double = 2141442
    options.InitialState (1,2) double = [2.3333865986017348, 0]
    options.HorizonTau (1,1) double {mustBePositive} = 25
    options.DTau (1,1) double {mustBePositive} = 0.001
    options.RelTol (1,1) double {mustBePositive} = 1e-9
    options.AbsTol (1,1) double {mustBePositive} = 1e-9
    options.HankelRows (1,1) double {mustBeInteger,mustBePositive} = 20000
    options.HankelColumns (1,1) double {mustBeInteger,mustBePositive} = 2000
    options.HDMDTolerance (1,1) double {mustBePositive} = 1e-10
    options.SaveHDMDBasis (1,1) logical = false
    options.SaveSINDyEigenfunctionValues (1,1) logical = true
end

suiteRoot = char(options.SuiteRoot);
assert(isfolder(suiteRoot), "ParityExport:MissingSuite", ...
    "MATLAB DynSysLab SuiteRoot does not exist: %s", suiteRoot);

if options.OutputFile == ""
    stamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
    outputFile = fullfile(pwd, "matlab_rpe_hdmd_sindy_parity_" + stamp + ".mat");
else
    outputFile = char(options.OutputFile);
end
outputFile = char(java.io.File(outputFile).getCanonicalPath());
[outputDir, outputStem] = fileparts(outputFile);
if isempty(outputDir), outputDir = pwd; end
if ~isfolder(outputDir), mkdir(outputDir); end
summaryFile = fullfile(outputDir, outputStem + ".json");

assert(~is_path_inside(outputFile, suiteRoot), "ParityExport:UnsafeOutput", ...
    "OutputFile must be outside SuiteRoot so the source suite remains untouched.");

originalPath = path;
pathCleanup = onCleanup(@() path(originalPath));
addpath(genpath(suiteRoot));

requiredSymbols = ["Bubble", "dynsys", "dynbatch", "Koopman", "DMD", "SINDy"];
for symbol = requiredSymbols
    assert(~isempty(which(symbol)), "ParityExport:MissingSymbol", ...
        "Required suite symbol '%s' was not found after adding SuiteRoot.", symbol);
end

rngStateBefore = rng;
rng(options.Seed, "twister");

timings = struct("trajectory_seconds", NaN, "hankel_rebuild_seconds", NaN, ...
    "hdmd_seconds", NaN, "sindy_seconds", NaN, "export_seconds", NaN);
exportClock = tic;

if isempty(options.ExistingSystem)
    [sys, timings.trajectory_seconds] = make_canonical_system(options);
    sourceMode = "fresh_canonical_run";
else
    sys = options.ExistingSystem;
    assert(isa(sys, "dynsys"), "ParityExport:BadExistingSystem", ...
        "ExistingSystem must be a computed MATLAB DynSysLab dynsys object.");
    assert(~isempty(sys.sols) && ~isempty(sys.sols.pool), ...
        "ParityExport:UncomputedSystem", "ExistingSystem has no computed trajectories.");
    sourceMode = "existing_system";
end

trajectory = extract_trajectory(sys);
assert(numel(trajectory.tau) >= options.HankelRows + options.HankelColumns, ...
    "ParityExport:ShortTrajectory", ...
    ["The trajectory has %d samples, but the requested H-DMD window needs " ...
     "at least rows + columns = %d samples."], ...
    numel(trajectory.tau), options.HankelRows + options.HankelColumns);

physicalParameters = scalar_table_to_struct(sys.dyn.prms);
dimensionlessGroups = scalar_table_to_struct(sys.dyn.dimgrps);
solver = extract_solver(sys, options);

% Rebuild the exact H-DMD input matrices used by the suite. The large X/Y
% arrays are not saved by default; fingerprints and boundary slices are.
hankelClock = tic;
[X, Y, hankelInfo] = rebuild_hdmd_matrices(sys, options);
timings.hankel_rebuild_seconds = toc(hankelClock);

koop = Koopman;
radiusObservable = sys.dyn.obsv(1);

% Koopman.process contains a handler for a tol option but omits tol from
% its arguments block. The suite therefore exposes only DMD.HDMD's 1e-10
% default through this public pipeline.
assert(options.HDMDTolerance == 1e-10, "ParityExport:UnsupportedHDMDTolerance", ...
    ["This MATLAB suite's public Koopman.process API cannot pass tol. " ...
     "Use the implemented/default H-DMD tolerance 1e-10 for parity."]);

hdmdClock = tic;
romHDMD = koop.process(data=sys, ...
    eigenvalue_method='HDMD', ...
    eigenfunction_method='HDMD', ...
    observable=radiusObservable, ...
    rows=options.HankelRows, ...
    cols=options.HankelColumns, ...
    stack_method='horizontal', ...
    sortmode='Energy');
timings.hdmd_seconds = toc(hdmdClock);

sindyClock = tic;
romSINDy = koop.process(data=sys, ...
    eigenvalue_method='provided', ...
    spectrum=romHDMD.spec, ...
    eigenfunction_method='SINDy', ...
    data_type='trajectory');
timings.sindy_seconds = toc(sindyClock);

hdmd = extract_hdmd(romHDMD, X, Y, hankelInfo, options);
sindy = extract_sindy(romSINDy, romHDMD, sys, options);

provenance = struct;
provenance.created_utc = char(datetime("now", "TimeZone", "UTC", ...
    "Format", "yyyy-MM-dd'T'HH:mm:ss.SSSXXX"));
provenance.wrapper_file = mfilename("fullpath");
provenance.suite_root = suiteRoot;
provenance.source_mode = char(sourceMode);
provenance.matlab_version = version;
provenance.matlab_release = version('-release');
provenance.computer = computer;
provenance.seed = options.Seed;
provenance.rng_before = rngStateBefore;
provenance.rng_after = rng;
provenance.resolved_files = struct( ...
    "Bubble", which("Bubble"), "dynsys", which("dynsys"), ...
    "dynbatch", which("dynbatch"), "Koopman", which("Koopman"), ...
    "DMD", which("DMD"), "SINDy", which("SINDy"));

parity = struct;
parity.schema = "dynsyslab.matlab-rpe-hdmd-sindy-parity";
parity.schema_version = 1;
parity.provenance = provenance;
parity.physical_parameters = physicalParameters;
parity.dimensionless_groups = dimensionlessGroups;
parity.time_scaling = extract_time_scaling(sys);
parity.solver = solver;
parity.trajectory = trajectory;
parity.hdmd = hdmd;
parity.sindy = sindy;
parity.timings = timings;
parity.notes = [ ...
    "States and time are saved in the suite's dimensionless representation."; ...
    "H-DMD eigenvalues are saved in discrete and continuous forms."; ...
    "SINDy is fit against the continuous-time spectrum."; ...
    "The authoritative data are in this MAT file; JSON is summary-only." ...
    ];

if options.SaveHDMDBasis
    parity.hdmd.input_X = X;
    parity.hdmd.input_Y = Y;
end

timings.export_seconds = toc(exportClock);
parity.timings = timings;
save(outputFile, "parity", "-v7.3");

summary = make_summary(parity, outputFile);
write_json(summaryFile, summary);

fprintf("\nParity export complete.\n  MAT:  %s\n  JSON: %s\n", ...
    outputFile, summaryFile);
end


function [sys, elapsed] = make_canonical_system(options)
dyn = Bubble;
sys = dynsys(dynamics=dyn, trajectories=1);

% Bubble's parameter table is initialized by dynsys.load_dynamics, not by
% its bare constructor. Verify rather than overwrite the canonical suite
% defaults; this makes an unexpected suite revision fail loudly.
expected = struct("R0",1e-5, "K",4/3, "p0",1e5, "pV",2330, ...
    "rhoL",998, "sigma",0.0725, "muL",0.001);
names = fieldnames(expected);
for k = 1:numel(names)
    name = names{k};
    assert(isequal(sys.dyn.prms.(name), expected.(name)), ...
        "ParityExport:UnexpectedParameter", ...
        "Bubble default %s is %.17g, expected %.17g.", ...
        name, sys.dyn.prms.(name), expected.(name));
end

% These assignments match generate_training.m and are intentionally made
% after dynsys initialization. They are inactive for this RPE comparison.
sys.dyn.prms.CL = 0;
sys.dyn.prms.CG = 0;
sys.load_time({[0, options.HorizonTau; 0, options.HorizonTau], options.DTau});
sys.load_inits(options.InitialState);
sys.options = odeset(sys.options, "RelTol", options.RelTol, "AbsTol", options.AbsTol);

batch = dynbatch(systems={sys});
clock = tic;
batch.compute("soft");
elapsed = toc(clock);
end


function trajectory = extract_trajectory(sys)
ts = sys.sols.pool{1};
trajectory = struct;
trajectory.tau = double(ts.t(:));
trajectory.states = double(ts.ts.');
trajectory.state_layout = "rows=samples, columns=[dimensionless radius, dimensionless radial velocity]";
trajectory.initial_state = trajectory.states(1,:);
trajectory.final_state = trajectory.states(end,:);
trajectory.sample_count = size(trajectory.states,1);
trajectory.state_dimension = size(trajectory.states,2);
if numel(trajectory.tau) > 1
    trajectory.dtau_median = median(diff(trajectory.tau));
    trajectory.dtau_max_deviation = max(abs(diff(trajectory.tau) - trajectory.dtau_median));
else
    trajectory.dtau_median = NaN;
    trajectory.dtau_max_deviation = NaN;
end
end


function solver = extract_solver(sys, options)
solver = struct;
solver.requested_method = "ode15s through dynbatch.compute(""soft"")";
solver.actual_rel_tol = odeget(sys.options, "RelTol");
solver.actual_abs_tol = odeget(sys.options, "AbsTol");
solver.requested_rel_tol = options.RelTol;
solver.requested_abs_tol = options.AbsTol;
solver.time_spread = char(sys.time_spread);
solver.requested_dtau = options.DTau;
end


function scaling = extract_time_scaling(sys)
p = sys.dyn.prms;
pG0 = p.p0 - p.pV + 2*p.sigma/p.R0;
f0 = (1/(2*pi*p.R0*sqrt(p.rhoL))) * sqrt(3*p.K*pG0 - 2*p.sigma/p.R0);
scaling = struct;
scaling.pG0_Pa = pG0;
scaling.f0_Hz = f0;
scaling.characteristic_time_seconds = 1/f0;
scaling.formula = "f0=(1/(2*pi*R0*sqrt(rhoL)))*sqrt(3*K*pG0-2*sigma/R0); tau=t*f0";
scaling.dimensionless_time_symbol = "tau";
end


function [X, Y, info] = rebuild_hdmd_matrices(sys, options)
selected = sys.tindex(2);
sols = sys.sols.pool{1}.ts;
signal = sys.dyn.obsv(1).g(sols);

n = options.HankelColumns;
m = options.HankelRows;
N = n + m - 1;
LS = length(sols) - 1;
M = min(N, LS);
idx1 = 1:n;
idx2 = n:M;

C1 = signal(:,selected(idx1));
R1 = signal(:,selected(idx2)).';
C2 = signal(:,selected(idx1+1));
R2 = signal(:,selected(idx2+1)).';
X = hankel(C1(1,:), R1(:,1)).';
Y = hankel(C2(1,:), R2(:,1)).';

info = struct;
info.observable = "radius (sys.dyn.obsv(1))";
info.stack_method = "horizontal";
info.sort_mode = "Energy";
info.requested_rows = m;
info.requested_columns = n;
info.N_rows_plus_columns_minus_one = N;
info.effective_M = M;
info.selected_time_indices = double(selected(:));
info.signal_used = double(signal(:,selected).');
info.X_size = size(X);
info.Y_size = size(Y);
info.X_frobenius_norm = norm(X, "fro");
info.Y_frobenius_norm = norm(Y, "fro");
info.X_first_block = X(1:min(5,end), 1:min(5,end));
info.X_last_block = X(max(1,end-4):end, max(1,end-4):end);
info.Y_first_block = Y(1:min(5,end), 1:min(5,end));
info.Y_last_block = Y(max(1,end-4):end, max(1,end-4):end);
info.X_first_column = X(:,1);
info.X_last_column = X(:,end);
info.Y_first_column = Y(:,1);
info.Y_last_column = Y(:,end);
info.construction = [ ...
    "N=columns+rows-1; M=min(N,length(solution)-1)"; ...
    "X=hankel(signal(1:columns),signal(columns:M))'"; ...
    "Y=hankel(signal(2:columns+1),signal(columns+1:M+1))'" ...
    ];
end


function hdmd = extract_hdmd(rom, X, Y, hankelInfo, options)
[psi, lambdaDFromObjects, lambdaCFromObjects, nu] = expand_krom(rom);
lambdaD = reshape(rom.spec.get_spec("Discrete", "Full"), [], 1);
lambdaC = reshape(rom.spec.get_spec("Continuous", "Full"), [], 1);

assert(numel(lambdaD) == size(psi,2), "ParityExport:HDMDShape", ...
    "Expanded H-DMD spectrum and eigenfunction counts differ.");

energy = vecnorm(nu, 2, 1).';
[energySorted, energyOrder] = sort(energy, "descend");

singularValuesEnergyOrder = reshape(double(rom.aux{2}), [], 1);
singularValuesDescending = sort(singularValuesEnergyOrder, "descend");

hdmd = struct;
hdmd.settings = struct( ...
    "rows", options.HankelRows, ...
    "columns", options.HankelColumns, ...
    "svd_absolute_tolerance", options.HDMDTolerance, ...
    "stack_method", "horizontal", ...
    "sort_mode", "Energy", ...
    "observable", "dimensionless radius");
hdmd.hankel = hankelInfo;
hdmd.retained_rank = numel(lambdaD);
hdmd.full_spectrum_retained = true;
hdmd.discrete_eigenvalues_public_order = lambdaD;
hdmd.continuous_eigenvalues_public_order = lambdaC;
hdmd.discrete_eigenvalues_object_order = lambdaDFromObjects;
hdmd.continuous_eigenvalues_object_order = lambdaCFromObjects;
hdmd.projected_eigenfunctions_public_order = psi;
hdmd.terminal_amplitudes_public_order = nu;
hdmd.energy_public_order = energy;
hdmd.energy_sort_indices_from_public_order = energyOrder;
hdmd.energy_sorted_values = energySorted;
hdmd.discrete_eigenvalues_energy_order = lambdaD(energyOrder);
hdmd.continuous_eigenvalues_energy_order = lambdaC(energyOrder);
hdmd.projected_eigenfunctions_energy_order = psi(:,energyOrder);
hdmd.terminal_amplitudes_energy_order = nu(:,energyOrder);
hdmd.singular_values_stored_energy_order = singularValuesEnergyOrder;
hdmd.singular_values_descending = singularValuesDescending;
hdmd.singular_values_are_complete_for_econ_svd = ...
    numel(singularValuesDescending) == min(size(X));
hdmd.svd_formula = "[U,S,V]=svd(X,'econ'); retain diag(S)>tol";
hdmd.reduced_operator_formula = "Atilde=((U'*Y)*V)*diag(1./diag(S))";
hdmd.continuous_conversion = "lambda_continuous=log(lambda_discrete)/dtau";
hdmd.X_frobenius_norm_check = norm(X, "fro");
hdmd.Y_frobenius_norm_check = norm(Y, "fro");

if options.SaveHDMDBasis
    hdmd.U_stored_energy_order = rom.aux{1};
    hdmd.V_stored_energy_order = rom.aux{3};
end
end


function sindy = extract_sindy(rom, sourceRom, sys, options)
[coefficients, lambdaDObjects, lambdaCObjects, ~] = expand_krom(rom);
lambdaD = reshape(sourceRom.spec.get_spec("Discrete", "Full"), [], 1);
lambdaC = reshape(sourceRom.spec.get_spec("Continuous", "Full"), [], 1);

sind = rom.psi(1).SINDy;
serialized = sys.sols.serialize;
stateData = serialized.ts.';
dstateData = sys.dyn.eqn(0, stateData.', 0).';
Theta = sind.buildTheta(stateData);
Gamma = sind.buildGamma(stateData, dstateData);
exponents = monomial_exponents(sind.nvars, sind.polyorder);

assert(size(coefficients,1) == size(Theta,2), "ParityExport:SINDyShape", ...
    "SINDy coefficient count does not match the reconstructed library.");
assert(size(coefficients,2) == numel(lambdaC), "ParityExport:SINDyModes", ...
    "SINDy eigenfunction count does not match the supplied spectrum.");

residualNorm = zeros(numel(lambdaC),1);
relativeResidual = zeros(numel(lambdaC),1);
for k = 1:numel(lambdaC)
    operator = lambdaC(k)*Theta - Gamma;
    residualNorm(k) = norm(operator*coefficients(:,k));
    denominator = norm(operator,"fro") * norm(coefficients(:,k));
    relativeResidual(k) = residualNorm(k) / max(denominator, eps);
end

sindy = struct;
sindy.settings = struct( ...
    "polynomial_order", sind.polyorder, ...
    "coefficient_threshold", sind.tol, ...
    "uses_sines", logical(sind.usesine), ...
    "state_variable_count", sind.nvars, ...
    "includes_constant_term", false, ...
    "library_normalization", "none", ...
    "threshold_policy", "one pass: coefficients with abs(value)<tol set to zero");
sindy.input_discrete_eigenvalues = lambdaD;
sindy.input_continuous_eigenvalues = lambdaC;
sindy.object_discrete_eigenvalues = lambdaDObjects;
sindy.object_continuous_eigenvalues = lambdaCObjects;
sindy.discrete_to_continuous_formula = "lambda=log(mu)/dtau";
sindy.state_data = stateData;
sindy.state_derivatives_from_rpe_equation = dstateData;
sindy.library_theta = Theta;
sindy.library_directional_derivative_gamma = Gamma;
sindy.monomial_exponents = exponents;
sindy.monomial_description = "Each row is [power(x1), power(x2), ...], degree 1 through polynomial_order; no constant.";
sindy.coefficients = coefficients;
sindy.absolute_residual_norm = residualNorm;
sindy.relative_residual_norm = relativeResidual;
sindy.solve_formula = "v=last right singular vector of (lambda*Theta-Gamma); then abs(v)<tol -> 0";

if options.SaveSINDyEigenfunctionValues
    sindy.eigenfunction_values_on_training_data = Theta * coefficients;
end
end


function [psi, lambdaD, lambdaC, nu] = expand_krom(rom)
psiCells = {};
lambdaD = [];
lambdaC = [];
nuCells = {};

for k = 1:numel(rom.psi)
    item = rom.psi(k);
    psiCells{end+1} = double(item.psi(:)); %#ok<AGROW>
    lambdaD(end+1,1) = item.lambda.d; %#ok<AGROW>
    lambdaC(end+1,1) = item.lambda.c; %#ok<AGROW>
    nuCells{end+1} = reshape(double(item.obsv.nu), [], 1); %#ok<AGROW>
    if ~isempty(item.dual)
        dual = item.dual;
        psiCells{end+1} = double(dual.psi(:)); %#ok<AGROW>
        lambdaD(end+1,1) = dual.lambda.d; %#ok<AGROW>
        lambdaC(end+1,1) = dual.lambda.c; %#ok<AGROW>
        nuCells{end+1} = reshape(double(dual.obsv.nu), [], 1); %#ok<AGROW>
    end
end

psi = [psiCells{:}];
nu = [nuCells{:}];
end


function exponents = monomial_exponents(nvars, maxDegree)
exponents = zeros(0,nvars);
for degree = 1:maxDegree
    combos = nondecreasing_indices(nvars, degree, 1);
    for row = 1:size(combos,1)
        powers = zeros(1,nvars);
        for column = 1:degree
            powers(combos(row,column)) = powers(combos(row,column)) + 1;
        end
        exponents(end+1,:) = powers; %#ok<AGROW>
    end
end
end


function combos = nondecreasing_indices(nvars, degree, minimum)
if degree == 1
    combos = (minimum:nvars).';
    return
end
combos = zeros(0,degree);
for first = minimum:nvars
    tails = nondecreasing_indices(nvars, degree-1, first);
    combos = [combos; [repmat(first,size(tails,1),1), tails]]; %#ok<AGROW>
end
end


function value = scalar_table_to_struct(tbl)
if istable(tbl)
    value = table2struct(tbl, "ToScalar", true);
elseif isstruct(tbl)
    value = tbl;
else
    value = struct;
    names = properties(tbl);
    for k = 1:numel(names), value.(names{k}) = tbl.(names{k}); end
end
end


function summary = make_summary(parity, outputFile)
summary = struct;
summary.schema = parity.schema;
summary.schema_version = parity.schema_version;
summary.mat_file = outputFile;
summary.provenance = rmfield(parity.provenance, ["rng_before", "rng_after"]);
summary.physical_parameters = parity.physical_parameters;
summary.dimensionless_groups = parity.dimensionless_groups;
summary.time_scaling = parity.time_scaling;
summary.solver = parity.solver;
summary.trajectory = struct( ...
    "sample_count", parity.trajectory.sample_count, ...
    "initial_state", parity.trajectory.initial_state, ...
    "final_state", parity.trajectory.final_state, ...
    "dtau_median", parity.trajectory.dtau_median, ...
    "dtau_max_deviation", parity.trajectory.dtau_max_deviation);
summary.hdmd = struct( ...
    "settings", parity.hdmd.settings, ...
    "X_size", parity.hdmd.hankel.X_size, ...
    "Y_size", parity.hdmd.hankel.Y_size, ...
    "X_frobenius_norm", parity.hdmd.hankel.X_frobenius_norm, ...
    "Y_frobenius_norm", parity.hdmd.hankel.Y_frobenius_norm, ...
    "retained_rank", parity.hdmd.retained_rank, ...
    "singular_values_are_complete_for_econ_svd", ...
        parity.hdmd.singular_values_are_complete_for_econ_svd);
summary.sindy = struct( ...
    "settings", parity.sindy.settings, ...
    "library_size", size(parity.sindy.library_theta), ...
    "coefficient_matrix_size", size(parity.sindy.coefficients), ...
    "maximum_absolute_residual_norm", max(parity.sindy.absolute_residual_norm), ...
    "maximum_relative_residual_norm", max(parity.sindy.relative_residual_norm));
summary.timings = parity.timings;
summary.notes = parity.notes;
end


function write_json(filename, value)
encoded = jsonencode(value, "PrettyPrint", true);
fid = fopen(filename, "w");
assert(fid >= 0, "ParityExport:CannotWriteJSON", ...
    "Could not open JSON output for writing: %s", filename);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s\n", encoded);
end


function tf = is_path_inside(candidate, root)
candidate = lower(char(java.io.File(candidate).getCanonicalPath()));
root = lower(char(java.io.File(root).getCanonicalPath()));
if ~endsWith(root, filesep), root = [root filesep]; end
tf = startsWith(candidate, root);
end
