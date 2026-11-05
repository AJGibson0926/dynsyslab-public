function summary = export_playground_duffing_lambdaopt_sindy_regression(outputFile, options)
%EXPORT_PLAYGROUND_DUFFING_LAMBDAOPT_SINDY_REGRESSION Reproduce the
% cleanup-suite playground.m Duffing identification and B-regression path.
%
% This wrapper is external to, and does not modify, the MATLAB suite.  It
% reconstructs noisyDuffData.csv directly with readcell and the documented
% dynbatch/dynsys/dynclst/inst/time/sol/ctrl/u record grammar.  It never calls
% data_node.load, uigetfile, or any other file-selection UI.
%
% The free data use the same 30 retained ICs, parameters, time mesh, and
% the suite's own bat.compute("hard") execution path.  The suite's own
% Koopman.process implementation performs lambda optimization and SINDy.
% get_B_parts is expanded here statement-for-statement because that method
% contains an unconditional keyboard breakpoint in the read-only source.

arguments
    outputFile (1,1) string
    options.SuiteRoot (1,1) string = ...
        "~/Documents/MATLAB/dynsyslab (5-15-25) (cleanup begins)"
    options.InitialConditionFile (1,1) string = ...
        "~/Documents/MATLAB/dynsyslab (5-15-25) (cleanup begins)\scripts\dynsyslib\Duffing\DuffingICs.csv"
    options.ControlCSVFile (1,1) string = ...
        "~/Documents/MATLAB/dynsyslab (5-15-25) (cleanup begins)\noisyDuffData.csv"
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 2141442
    options.FreeTrajectoryCount (1,1) double {mustBeInteger,mustBePositive} = 30
    options.FreeDeltaT (1,1) double {mustBePositive} = 0.01
    options.FreeFinalTime (1,1) double {mustBePositive} = 8
    options.OptimizationOrder (1,1) double {mustBeInteger,mustBePositive} = 50
    options.SaveRegressionParts (1,1) logical = true
end

startedAt = datetime("now", "TimeZone", "UTC");
suiteRoot = canonical_path(options.SuiteRoot);
initialConditionFile = canonical_path(options.InitialConditionFile);
controlCSVFile = canonical_path(options.ControlCSVFile);
outputFile = canonical_path(outputFile);
assert(isfolder(suiteRoot), "MATLAB suite root not found: %s", suiteRoot);
assert(isfile(initialConditionFile), "Duffing IC file not found: %s", initialConditionFile);
assert(isfile(controlCSVFile), "Controlled Duffing CSV not found: %s", controlCSVFile);
assert(~startsWith(lower(outputFile), lower(suiteRoot + filesep)), ...
    "Output must remain outside the read-only MATLAB source suite.");
outputDirectory = fileparts(outputFile);
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end

playgroundFile = fullfile(suiteRoot, "playground.m");
koopmanFile = fullfile(suiteRoot, "scripts", "comp_modules", "Koopman.m");
dmdFile = fullfile(suiteRoot, "scripts", "comp_modules", "DMD.m");
sindyFile = fullfile(suiteRoot, "scripts", "comp_modules", "SINDy.m");
kromFile = fullfile(suiteRoot, "scripts", "data_structs", "Koopman_tree", "KROM.m");
loaderFile = fullfile(suiteRoot, "scripts", "data_structs", "data_tree", "data_node.m");
duffingFile = fullfile(suiteRoot, "scripts", "dynsyslib", "Duffing", "Duffing.m");
computeFile = fullfile(suiteRoot, "scripts", "comp_modules", "compute_node.m");
sourceFiles = [playgroundFile; koopmanFile; dmdFile; sindyFile; kromFile; ...
    loaderFile; duffingFile; computeFile; initialConditionFile; controlCSVFile];
sourceHashesBefore = strings(size(sourceFiles));
for fileIndex = 1:numel(sourceFiles)
    sourceHashesBefore(fileIndex) = sha256_file(sourceFiles(fileIndex));
end

originalFolder = string(pwd);
folderCleanup = onCleanup(@() cd(originalFolder));
addpath(genpath(suiteRoot));
cd(suiteRoot);

timings = struct;
stageStart = tic;
rng(options.Seed, "twister");
initialConditionsAll = readmatrix(initialConditionFile).';
assert(size(initialConditionsAll, 1) >= options.FreeTrajectoryCount, ...
    "Duffing IC file does not contain the requested number of trajectories.");
free_initial_conditions = initialConditionsAll(1:options.FreeTrajectoryCount, :);
assert(size(free_initial_conditions, 2) == 2, "Expected two-state Duffing ICs.");

duff = dynsys(dynamics=Duffing, trajectories=[options.FreeTrajectoryCount]);
duff.dyn.prms.alpha = 4;
duff.dyn.prms.beta = -1;
duff.dyn.prms.gamma = 1/2;
duff.load_inits(free_initial_conditions);
duff.load_time({[0 options.FreeFinalTime; 0 options.FreeFinalTime; ...
    0 options.FreeFinalTime], options.FreeDeltaT});
freeBatch = dynbatch(systems={duff});
freeBatch.compute("hard");
free_time = duff.clst{1}.inst{1}.sols.t;
free_trajectories = zeros(2, numel(free_time), options.FreeTrajectoryCount);
for trajectoryIndex = 1:options.FreeTrajectoryCount
    instance = duff.clst{1}.inst{trajectoryIndex};
    assert(isequal(instance.sols.t, free_time), ...
        "Suite hard-compute time mesh changed at free trajectory %d.", trajectoryIndex);
    free_trajectories(:, :, trajectoryIndex) = instance.sols.ts;
end
timings.free_training_seconds = toc(stageStart);
assert(numel(duff.sols.pool) == options.FreeTrajectoryCount, ...
    "Free-data hierarchy did not pool all trajectories.");
assert(length(duff.sols.ts) == numel(free_time), ...
    "The lambda optimizer's trajectory-length contract is not satisfied.");

stageStart = tic;
noisyduff = load_legacy_dynbatch_csv_direct(controlCSVFile);
timings.control_csv_decode_seconds = toc(stageStart);
controlledSystemCount = numel(noisyduff.sys);
assert(controlledSystemCount == 100, ...
    "Expected 100 controlled systems in noisyDuffData.csv, decoded %d.", controlledSystemCount);
[controlled_states, controlled_controls_sample_aligned, controlled_time] = ...
    collect_controlled_data(noisyduff);
controlled_interval_controls = controlled_controls_sample_aligned(:, 2:end, :);
controlled_dynamics_parameter_names = string( ...
    noisyduff.sys{end}.dyn.prms.Properties.VariableNames);
controlled_dynamics_parameter_values = cell2mat(table2cell( ...
    noisyduff.sys{end}.dyn.prms));

Koop = Koopman;
optstrct = struct;
optstrct.order = options.OptimizationOrder;
stageStart = tic;
duffrom = Koop.process(data=duff, eigenvalue_method='lambda-opt', ...
    eigenfunction_method='g-function', aux=optstrct);
timings.lambda_optimization_and_boundary_functions_seconds = toc(stageStart);
optimized_eigenvalues_discrete_full = duffrom.spec.get_spec("Discrete", "Full").';
optimized_eigenvalues_continuous_full = duffrom.spec.get_spec("Continuous", "Full").';
optimized_eigenvalues_discrete_compact = duffrom.spec.get_spec("Discrete", "Compact").';
optimized_eigenvalues_continuous_compact = duffrom.spec.get_spec("Continuous", "Compact").';

stageStart = tic;
rom = Koop.process(data=duff, eigenvalue_method='provided', spectrum=duffrom.spec, ...
    eigenfunction_method='SINDy', data_type='domain');
timings.sindy_eigenfunctions_seconds = toc(stageStart);
R = rom.rank;
sindy_eigenvalues_discrete_full = rom.spec.get_spec("Discrete", "Full").';
sindy_eigenvalues_continuous_full = rom.spec.get_spec("Continuous", "Full").';
sindy_eigenvalues_discrete_compact = rom.spec.get_spec("Discrete", "Compact").';
sindy_eigenvalues_continuous_compact = rom.spec.get_spec("Continuous", "Compact").';
assert(numel(sindy_eigenvalues_continuous_compact) == R, ...
    "Compact SINDy spectrum and no-dual KROM rank do not align.");
sindy_coefficients = complex(zeros(numel(rom.psi(1).psi), R));
for modeIndex = 1:R
    sindy_coefficients(:, modeIndex) = rom.psi(modeIndex).psi;
end
sindy_library_order = rom.mthd.aux.polyorder;
sindy_library_state_dimension = rom.mthd.aux.nvars;
sindy_threshold = rom.mthd.aux.tol;
sindy_domain = duff.sols.dom;
sindy_domain_points = sindy_domain.x;
sindy_domain_boundary = sindy_domain.boundary;
sindy_domain_density = sindy_domain.density;

stageStart = tic;
[A_continuous, Z_controlled, DZ_controlled, Y_controlled, U_controlled] = ...
    get_B_parts_without_breakpoint(rom, noisyduff);
lifted_control_residual = DZ_controlled - A_continuous * Z_controlled;
bilinear_coefficients = lifted_control_residual * pinv([Y_controlled; U_controlled]);
B1_bilinear = bilinear_coefficients(:, 1:R);
B2_bilinear = bilinear_coefficients(:, end);
B_constant = lifted_control_residual * pinv(U_controlled);
bilinear_fit = B1_bilinear * Y_controlled + B2_bilinear * U_controlled;
constant_fit = B_constant * U_controlled;
bilinear_relative_residual = norm(lifted_control_residual - bilinear_fit, "fro") / ...
    norm(lifted_control_residual, "fro");
constant_relative_residual = norm(lifted_control_residual - constant_fit, "fro") / ...
    norm(lifted_control_residual, "fro");
timings.control_regressions_seconds = toc(stageStart);

stageStart = tic;
free_state_samples = duff.sols.serialize().ts;
Z_free = rom.state(free_state_samples, form="complex");
C_projection = free_state_samples * pinv(Z_free);
C_native_before_projection = rom.get_C(form="complex", use_duals=0);
free_reconstruction_projection = C_projection * Z_free;
free_reconstruction_rmse_projection = sqrt(mean(abs( ...
    free_reconstruction_projection - free_state_samples).^2, 2));
timings.observable_projection_seconds = toc(stageStart);

sourceHashesAfter = strings(size(sourceFiles));
for fileIndex = 1:numel(sourceFiles)
    sourceHashesAfter(fileIndex) = sha256_file(sourceFiles(fileIndex));
end
assert(isequal(sourceHashesBefore, sourceHashesAfter), ...
    "A read-only MATLAB source/input hash changed during the wrapper run.");

completedAt = datetime("now", "TimeZone", "UTC");
metadata = struct( ...
    "schema", "dynsyslab.matlab-playground-duffing-lambdaopt-sindy-regression/v1", ...
    "candidate_only", true, ...
    "mutates_existing_matlab_suite", false, ...
    "mutates_existing_dynsyslab_framework", false, ...
    "data_node_load_called", false, ...
    "file_picker_called", false, ...
    "control_csv_decode", "direct readcell tagged-record state machine", ...
    "get_B_parts_execution", "statement-equivalent external expansion; source breakpoint skipped", ...
    "seed", options.Seed, ...
    "free_trajectory_count", options.FreeTrajectoryCount, ...
    "free_delta_t", options.FreeDeltaT, ...
    "free_final_time", options.FreeFinalTime, ...
    "free_solver", "suite dynbatch.compute hard mode", ...
    "duffing_parameters", struct("alpha", 4, "beta", -1, "gamma", 1/2), ...
    "controlled_derivative_parameters", struct( ...
        "names", controlled_dynamics_parameter_names, ...
        "values", controlled_dynamics_parameter_values, ...
        "source", "Duffing defaults reconstructed by legacy CSV loader"), ...
    "legacy_sindy_state_evaluation", ...
        "conjugates polynomial values through MATLAB apostrophe transpose", ...
    "optimization_order_requested", options.OptimizationOrder, ...
    "optimized_full_spectrum_count", numel(optimized_eigenvalues_discrete_full), ...
    "sindy_full_spectrum_count", numel(sindy_eigenvalues_discrete_full), ...
    "sindy_compact_krom_rank", R, ...
    "sindy_data_type", "domain", ...
    "sindy_polynomial_order", sindy_library_order, ...
    "sindy_threshold", sindy_threshold, ...
    "controlled_system_count", controlledSystemCount, ...
    "controlled_sample_count", size(controlled_states, 2), ...
    "control_alignment_sample", "CTRL.hist includes leading zero initializer", ...
    "control_alignment_interval", "CTRL.hist samples 2:end", ...
    "regression_form", "complex compact, use_duals=false", ...
    "source_files", sourceFiles, ...
    "source_hashes_before", sourceHashesBefore, ...
    "source_hashes_after", sourceHashesAfter, ...
    "source_files_unchanged", true, ...
    "wrapper_file", string(mfilename("fullpath") + ".m"), ...
    "wrapper_sha256", sha256_file(mfilename("fullpath") + ".m"), ...
    "started_at_utc", char(startedAt), ...
    "completed_at_utc", char(completedAt));

saveArguments = { ...
    "metadata", "timings", "free_time", "free_initial_conditions", ...
    "free_trajectories", "free_state_samples", "controlled_time", ...
    "controlled_states", "controlled_controls_sample_aligned", ...
    "controlled_interval_controls", "optimized_eigenvalues_discrete_full", ...
    "controlled_dynamics_parameter_names", ...
    "controlled_dynamics_parameter_values", ...
    "optimized_eigenvalues_continuous_full", ...
    "optimized_eigenvalues_discrete_compact", ...
    "optimized_eigenvalues_continuous_compact", ...
    "sindy_eigenvalues_discrete_full", "sindy_eigenvalues_continuous_full", ...
    "sindy_eigenvalues_discrete_compact", ...
    "sindy_eigenvalues_continuous_compact", "sindy_coefficients", ...
    "sindy_library_order", "sindy_library_state_dimension", "sindy_threshold", ...
    "sindy_domain_points", "sindy_domain_boundary", "sindy_domain_density", ...
    "A_continuous", "C_native_before_projection", "C_projection", ...
    "free_reconstruction_rmse_projection", "B1_bilinear", "B2_bilinear", ...
    "B_constant", "bilinear_relative_residual", "constant_relative_residual"};
if options.SaveRegressionParts
    saveArguments = [saveArguments, {"Z_free", "Z_controlled", "DZ_controlled", ...
        "Y_controlled", "U_controlled", "lifted_control_residual"}]; %#ok<AGROW>
end
save(outputFile, saveArguments{:}, "-v7");

matrixFile = replace(outputFile, ".mat", ".matrices.txt");
write_matrix_report(matrixFile, optimized_eigenvalues_discrete_full, ...
    optimized_eigenvalues_continuous_full, sindy_eigenvalues_discrete_compact, ...
    sindy_eigenvalues_continuous_compact, A_continuous, C_native_before_projection, ...
    C_projection, B1_bilinear, B2_bilinear, B_constant);

summary = struct( ...
    "schema", metadata.schema, ...
    "output_file", outputFile, ...
    "matrix_file", matrixFile, ...
    "output_sha256", sha256_file(outputFile), ...
    "matrix_sha256", sha256_file(matrixFile), ...
    "source_files_unchanged", true, ...
    "optimized_full_spectrum_count", numel(optimized_eigenvalues_discrete_full), ...
    "sindy_compact_krom_rank", R, ...
    "B1_shape", size(B1_bilinear), ...
    "B2_shape", size(B2_bilinear), ...
    "B_constant_shape", size(B_constant), ...
    "bilinear_relative_residual", bilinear_relative_residual, ...
    "constant_relative_residual", constant_relative_residual, ...
    "free_reconstruction_rmse_projection", free_reconstruction_rmse_projection.', ...
    "timings", timings);
jsonFile = replace(outputFile, ".mat", ".summary.json");
fid = fopen(jsonFile, "w");
assert(fid >= 0, "Could not create summary JSON: %s", jsonFile);
fileCleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s\n", jsonencode(summary, PrettyPrint=true));
clear fileCleanup
fprintf("Exported MATLAB playground Duffing parity bundle: %s\n", outputFile);
fprintf("Full-precision matrices: %s\n", matrixFile);
fprintf("Summary: %s\n", jsonFile);
clear folderCleanup
end


function batch = load_legacy_dynbatch_csv_direct(csvFile)
% Expand data_node.load's CSV branch without its uigetfile call.
tableCells = readcell(csvFile);
tableCells = tableCells(2:end, :);
batch = dynbatch();
system = [];
cluster = [];
instance = [];
controller = [];
trajectoryNumber = 0;
systemTime = [];
rowCount = size(tableCells, 1);
for rowIndex = 1:rowCount
    tag = tableCells{rowIndex, 1};
    if is_missing_cell(tag)
        continue
    end
    switch string(tag)
        case "dynbatch"
            % Header/count record only.
        case "dynsys"
            system = dynsys();
            systemNumber = tableCells{rowIndex, 2}; %#ok<NASGU>
            batch.sys{end + 1} = system;
            trajectoryNumber = 0;
            controller = [];
            dynamicsName = tableCells{rowIndex, 3};
            resolution = tableCells{rowIndex, 4};
            finalIndex = tableCells{rowIndex, 5};
            system.tres = resolution;
            system.time_struct = [[0 resolution * finalIndex]; [0 resolution * finalIndex]];
            system.tspan = 0:resolution:resolution * (finalIndex - 1);
            system.load_dynamics(eval(dynamicsName)); %#ok<EVLDOT>
            batch.children{end + 1} = system;
        case "dynclst"
            clusterNumber = tableCells{rowIndex, 2};
            if clusterNumber > numel(system.clst)
                system.add([[1]]);
            end
            cluster = system.clst{clusterNumber};
        case "inst"
            trajectoryNumber = tableCells{rowIndex, 2};
            if trajectoryNumber > numel(cluster.inst)
                cluster.add([[1]]);
            end
            instance = cluster.inst{trajectoryNumber};
        case "time"
            values = tableCells(rowIndex, 2:end);
            values = [values{:}];
            if trajectoryNumber == 0
                systemTime = values;
            else
                instance.sols.t = values;
            end
        case "sol"
            assert(trajectoryNumber > 0, "A sol record appeared outside an instance.");
            endRow = rowIndex;
            while endRow + 1 <= rowCount && is_missing_cell(tableCells{endRow + 1, 1})
                endRow = endRow + 1;
            end
            stateCells = tableCells(rowIndex:endRow, 2:end);
            instance.sols.ts = cell2mat(stateCells);
            if ~isempty(controller)
                instance.ctrl_hist.ts = controller.hist.ts;
                instance.ctrl_hist.t = controller.hist.t;
            end
        case "ctrl"
            controller = system.CTRL;
            controller.forcing = tableCells{rowIndex, 2};
            controller.subforcing = tableCells{rowIndex, 3};
        case "u"
            controller.hist.ts = [tableCells{rowIndex, 2:end}];
            controller.hist.t = systemTime;
        otherwise
            error("Unsupported legacy CSV tag '%s' at data row %d.", string(tag), rowIndex);
    end
end
batch.Pool();
end


function missing = is_missing_cell(value)
missing = isempty(value) || (isscalar(value) && ismissing(value));
end


function [states, controls, time] = collect_controlled_data(batch)
systemCount = numel(batch.sys);
sampleCount = length(batch.sys{1}.CTRL.hist.ts);
stateDimension = size(batch.sys{1}.sols.ts, 1);
inputDimension = size(batch.sys{1}.CTRL.hist.ts, 1);
states = zeros(stateDimension, sampleCount, systemCount);
controls = zeros(inputDimension, sampleCount, systemCount);
time = batch.sys{1}.CTRL.hist.t;
for systemIndex = 1:systemCount
    system = batch.sys{systemIndex};
    assert(size(system.sols.ts, 2) == sampleCount, ...
        "Controlled state sample count changed at system %d.", systemIndex);
    assert(isequal(system.CTRL.hist.t, time), ...
        "Controlled time mesh changed at system %d.", systemIndex);
    states(:, :, systemIndex) = system.sols.ts;
    controls(:, :, systemIndex) = system.CTRL.hist.ts;
end
end


function [A, Z, DZ, Y, U] = get_B_parts_without_breakpoint(rom, controlledBatch)
% Statement-equivalent expansion of KROM.get_B_parts for polynomial SINDy.
A = rom.get_A(mode='Continuous', form="complex", use_duals=0, sparse=1);
systemCount = numel(controlledBatch.sys);
stateCells = cell(1, systemCount);
controlCells = cell(1, systemCount);
for systemIndex = 1:systemCount
    system = controlledBatch.sys{systemIndex};
    stateCells{systemIndex} = system.sols.ts;
    controlCells{systemIndex} = system.CTRL.hist.ts;
end
X = [stateCells{:}];
U = [controlCells{:}];
sindy = rom.mthd.aux;
system = controlledBatch.sys{end};
DX = system.dyn.eqn(0, X, U);
DZZ = sindy.buildGamma(X.', DX.').';
Z = rom.state(X, form="complex");
R = size(Z, 1);
bilinearRows = cell(1, R);
derivativeRows = cell(1, R);
for modeIndex = 1:R
    bilinearRows{modeIndex} = Z(modeIndex, :) .* U;
    derivativeRows{modeIndex} = transpose(rom.psi(modeIndex).psi) * DZZ;
end
Y = vertcat(bilinearRows{:});
DZ = vertcat(derivativeRows{:});
end


function write_matrix_report(path, optimizedDiscrete, optimizedContinuous, ...
        sindyDiscrete, sindyContinuous, A, Cnative, Cprojection, B1, B2, Bconstant)
fid = fopen(path, "w");
assert(fid >= 0, "Could not create matrix report: %s", path);
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%% Full-precision MATLAB literals from the playground Duffing parity wrapper.\n");
fprintf(fid, "%% B1_bilinear and B2_bilinear satisfy B_psi(z) ~= B1*z + B2.\n");
fprintf(fid, "%% B_constant is the state-independent approximation.\n\n");
write_complex_matrix(fid, "optimized_eigenvalues_discrete_full", optimizedDiscrete);
write_complex_matrix(fid, "optimized_eigenvalues_continuous_full", optimizedContinuous);
write_complex_matrix(fid, "sindy_eigenvalues_discrete_compact", sindyDiscrete);
write_complex_matrix(fid, "sindy_eigenvalues_continuous_compact", sindyContinuous);
write_complex_matrix(fid, "A_continuous", full(A));
write_complex_matrix(fid, "C_native_before_projection", Cnative);
write_complex_matrix(fid, "C_projection", Cprojection);
write_complex_matrix(fid, "B1_bilinear", B1);
write_complex_matrix(fid, "B2_bilinear", B2);
write_complex_matrix(fid, "B_constant", Bconstant);
clear cleanup
end


function write_complex_matrix(fid, name, matrix)
fprintf(fid, "%s = [\n", name);
for rowIndex = 1:size(matrix, 1)
    for columnIndex = 1:size(matrix, 2)
        value = matrix(rowIndex, columnIndex);
        fprintf(fid, "%.17e%+.17ei", real(value), imag(value));
        if columnIndex < size(matrix, 2)
            fprintf(fid, ", ");
        end
    end
    if rowIndex < size(matrix, 1)
        fprintf(fid, ";");
    end
    fprintf(fid, "\n");
end
fprintf(fid, "];\n\n");
end


function path = canonical_path(path)
path = string(char(java.io.File(char(path)).getCanonicalPath()));
end


function digest = sha256_file(path)
fid = fopen(path, "rb");
assert(fid >= 0, "Could not open file for SHA-256: %s", path);
cleanup = onCleanup(@() fclose(fid));
bytes = fread(fid, Inf, "*uint8");
md = java.security.MessageDigest.getInstance("SHA-256");
md.update(typecast(bytes, "int8"));
raw = typecast(md.digest(), "uint8");
digest = lower(reshape(dec2hex(raw, 2).', 1, []));
clear cleanup
end
