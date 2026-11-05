function summary = export_korda_mezic_duffing_control_candidate(outputFile, options)
%EXPORT_KORDA_MEZIC_DUFFING_CONTROL_CANDIDATE Reproduce and export the
% precomputed-eigenvalue Duffing identification path from Korda and Mezic.
%
% This is an experiment wrapper.  It never edits the Korda--Mezic source
% tree and never edits DynSysLab.  The original script's nonstandard RK4
% fourth stage and its source-input scaling are reproduced intentionally so
% that the resulting bundle is a faithful comparison artifact.
%
% The wrapper stops after open-loop prediction and identification.  It does
% not enter the source script's interactive `keyboard` statements or run
% its MPC section.

arguments
    outputFile (1,1) string
    options.KordaRoot (1,1) string = "~/Documents/MATLAB/KordaMPC"
    options.Seed (1,1) double {mustBeInteger,mustBeNonnegative} = 2141442
    options.TrajectoryCount (1,1) double {mustBeInteger,mustBePositive} = 100
    options.TrajectoryLength (1,1) double {mustBeInteger,mustBeGreaterThan(options.TrajectoryLength,2)} = 800
    options.ControlTransitionCount (1,1) double {mustBeInteger,mustBePositive} = 200
    options.DeltaT (1,1) double {mustBePositive} = 0.01
    options.GridPoints (1,1) double {mustBeInteger,mustBeGreaterThan(options.GridPoints,4)} = 101
    options.SaveRegressionDesign (1,1) logical = true
end

startedAt = datetime("now", "TimeZone", "UTC");
kordaRoot = canonical_path(options.KordaRoot);
outputFile = canonical_path(outputFile);
sourceFile = fullfile(kordaRoot, "Duffing.m");
lambdaFile = fullfile(kordaRoot, "Functions", "lam_opt_duffing.mat");

assert(isfile(sourceFile), "Korda--Mezic Duffing.m was not found: %s", sourceFile);
assert(isfile(lambdaFile), "Precomputed eigenvalue file was not found: %s", lambdaFile);
assert(~startsWith(lower(outputFile), lower(kordaRoot + filesep)), ...
    "Output must be outside the read-only Korda--Mezic source tree.");

outputDir = fileparts(outputFile);
if ~isfolder(outputDir)
    mkdir(outputDir);
end

loaded = load(lambdaFile, "LAM_OPT");
LAM_OPT = loaded.LAM_OPT;
stateDimension = 2;
inputDimension = 1;
modeCount = sum(cellfun(@numel, LAM_OPT));
assert(modeCount == 20, "Expected the published 20-eigenfunction Duffing model.");
assert(numel(LAM_OPT) == stateDimension, ...
    "Expected one optimized eigenvalue block per physical state.");

eigenvalues_discrete = reshape(cell2mat(LAM_OPT), [], 1);
A_discrete = diag(eigenvalues_discrete);
A_continuous = diag(log(eigenvalues_discrete) ./ options.DeltaT);
C_official = zeros(stateDimension, modeCount);
offset = 0;
for stateIndex = 1:stateDimension
    count = numel(LAM_OPT{stateIndex});
    C_official(stateIndex, offset + (1:count)) = 1;
    offset = offset + count;
end

% Match the published random stream.  With precomputed lambdas enabled the
% source does not call rng(31), so random forcing follows IC generation.
rng(options.Seed, "twister");
free_trajectories = zeros(stateDimension, options.TrajectoryLength, options.TrajectoryCount);
free_initial_conditions = zeros(stateDimension, options.TrajectoryCount);
for trajectoryIndex = 1:options.TrajectoryCount
    x0 = randn(stateDimension, 1);
    x0 = x0 ./ norm(x0);
    free_initial_conditions(:, trajectoryIndex) = x0;
    free_trajectories(:, 1, trajectoryIndex) = x0;
    for timeIndex = 1:(options.TrajectoryLength - 1)
        free_trajectories(:, timeIndex + 1, trajectoryIndex) = ...
            source_step_bug_compatible( ...
                free_trajectories(:, timeIndex, trajectoryIndex), 0, options.DeltaT);
    end
end

% Reproduce the boundary-function regression and interpolation used by the
% source.  Ordering is trajectory-major, matching F_vec and bdiag.
X_free = reshape(free_trajectories, stateDimension, []);
Psi_free_regression = zeros(modeCount, size(X_free, 2));
boundary_values = cell(stateDimension, 1);
free_fit_residual_percent = zeros(stateDimension, 1);
modeOffset = 0;
for stateIndex = 1:stateDimension
    F = reshape(free_trajectories(stateIndex, :, :), [], 1);
    eigenvaluesForState = LAM_OPT{stateIndex}(:);
    modeCountForState = numel(eigenvaluesForState);
    blocks = cell(1, modeCountForState);
    for modeIndex = 1:modeCountForState
        temporalPowers = eigenvaluesForState(modeIndex) .^ (0:(options.TrajectoryLength - 1));
        blocks{modeIndex} = kron(speye(options.TrajectoryCount), temporalPowers(:));
    end
    Lbig = horzcat(blocks{:});
    normalMatrix = Lbig' * Lbig;
    normalRhs = Lbig' * F;
    g00 = normalMatrix \ normalRhs;
    boundary_values{stateIndex} = reshape(g00, options.TrajectoryCount, modeCountForState);
    fit = Lbig * g00;
    free_fit_residual_percent(stateIndex) = 100 * norm(fit - F) / norm(F);
    for modeIndex = 1:modeCountForState
        rows = (modeIndex - 1) * options.TrajectoryCount + (1:options.TrajectoryCount);
        Psi_free_regression(modeOffset + modeIndex, :) = (blocks{modeIndex} * g00(rows)).';
    end
    modeOffset = modeOffset + modeCountForState;
end

interpolants = cell(modeCount, 1);
for modeIndex = 1:modeCount
    interpolants{modeIndex} = scatteredInterpolant( ...
        X_free(1, :).', X_free(2, :).', Psi_free_regression(modeIndex, :).', ...
        "linear", "linear");
end
Psi_free_evaluated = lift_points(interpolants, X_free);
C_projection = X_free * pinv(Psi_free_evaluated);
free_reconstruction_official = real(C_official * Psi_free_evaluated);
free_reconstruction_projection = real(C_projection * Psi_free_evaluated);
free_reconstruction_rmse_official = sqrt(mean((free_reconstruction_official - X_free).^2, 2));
free_reconstruction_rmse_projection = sqrt(mean((free_reconstruction_projection - X_free).^2, 2));

% Reproduce the random-forced controlled dataset.  U_source is exactly the
% input appearing in the source script.  U_canonical = 0.5*U_source is the
% additive acceleration used by DynSysLab's canonical Duffing convention.
U_source_all = 2 * rand(options.ControlTransitionCount, options.TrajectoryCount) - 1;
X_controlled_all = zeros( ...
    stateDimension, options.ControlTransitionCount + 1, options.TrajectoryCount);
X_controlled_all(:, 1, :) = reshape(free_initial_conditions, stateDimension, 1, []);
for timeIndex = 1:options.ControlTransitionCount
    current = reshape(X_controlled_all(:, timeIndex, :), stateDimension, []);
    X_controlled_all(:, timeIndex + 1, :) = reshape( ...
        source_step_bug_compatible(current, U_source_all(timeIndex, :), options.DeltaT), ...
        stateDimension, 1, []);
end
accepted_mask = squeeze(all(all(abs(X_controlled_all) <= 3 & ~isnan(X_controlled_all), 1), 2)).';
accepted_indices = find(accepted_mask);
X_controlled = X_controlled_all(:, :, accepted_mask);
U_source = permute(U_source_all(:, accepted_mask), [3, 1, 2]);
U_canonical = 0.5 .* U_source;
acceptedTrajectoryCount = size(X_controlled, 3);
assert(acceptedTrajectoryCount > 0, "No controlled trajectories passed the source acceptance test.");

Psi_controlled = zeros(modeCount, options.ControlTransitionCount + 1, acceptedTrajectoryCount);
for trajectoryIndex = 1:acceptedTrajectoryCount
    Psi_controlled(:, :, trajectoryIndex) = lift_points( ...
        interpolants, X_controlled(:, :, trajectoryIndex));
end

% Compute chain-rule lifted derivatives at controlled interval starts.  This
% avoids contaminating bilinear regression with temporal differencing across
% the discontinuities of the interval-held random forcing.
controlled_points = reshape( ...
    X_controlled(:, 1:options.ControlTransitionCount, :), stateDimension, []);
controlled_inputs_source = reshape(U_source, 1, []);
positionStepControlled = max(1e-7, 1e-5 * range(X_free(1, :)));
velocityStepControlled = max(1e-7, 1e-5 * range(X_free(2, :)));
positionPlus = controlled_points; positionMinus = controlled_points;
velocityPlus = controlled_points; velocityMinus = controlled_points;
positionPlus(1, :) = positionPlus(1, :) + positionStepControlled;
positionMinus(1, :) = positionMinus(1, :) - positionStepControlled;
velocityPlus(2, :) = velocityPlus(2, :) + velocityStepControlled;
velocityMinus(2, :) = velocityMinus(2, :) - velocityStepControlled;
dPsi_dposition_controlled = (lift_points(interpolants, positionPlus) - ...
    lift_points(interpolants, positionMinus)) ./ (2 * positionStepControlled);
dPsi_dvelocity_controlled = (lift_points(interpolants, velocityPlus) - ...
    lift_points(interpolants, velocityMinus)) ./ (2 * velocityStepControlled);
controlled_flow_source = source_rhs(controlled_points, controlled_inputs_source);
lifted_derivatives_controlled = ...
    dPsi_dposition_controlled .* controlled_flow_source(1, :) + ...
    dPsi_dvelocity_controlled .* controlled_flow_source(2, :);
lifted_derivatives_controlled = reshape( ...
    lifted_derivatives_controlled, modeCount, options.ControlTransitionCount, []);
controlled_exact_Bpsi_source = reshape( ...
    0.5 .* dPsi_dvelocity_controlled, modeCount, options.ControlTransitionCount, []);
controlled_exact_Bpsi_canonical = 2 .* controlled_exact_Bpsi_source;

% Published multistep observed-output regression for the discrete B matrix.
CA_powers = cell(options.ControlTransitionCount + 1, 1);
Apower = eye(modeCount);
for powerIndex = 0:options.ControlTransitionCount
    CA_powers{powerIndex + 1} = C_official * Apower;
    Apower = Apower * A_discrete;
end
rowCount = acceptedTrajectoryCount * options.ControlTransitionCount * stateDimension;
regression_design = complex(zeros(rowCount, modeCount * inputDimension));
regression_target = complex(zeros(rowCount, 1));
rowOffset = 0;
for trajectoryIndex = 1:acceptedTrajectoryCount
    z0 = Psi_controlled(:, 1, trajectoryIndex);
    for horizonIndex = 1:options.ControlTransitionCount
        rows = rowOffset + (1:stateDimension);
        autonomous = CA_powers{horizonIndex + 1} * z0;
        regression_target(rows) = X_controlled(:, horizonIndex + 1, trajectoryIndex) - autonomous;
        accumulated = complex(zeros(stateDimension, modeCount));
        for inputIndex = 0:(horizonIndex - 1)
            remainingPower = horizonIndex - inputIndex - 1;
            accumulated = accumulated + ...
                U_source(1, inputIndex + 1, trajectoryIndex) .* CA_powers{remainingPower + 1};
        end
        regression_design(rows, :) = accumulated;
        rowOffset = rowOffset + stateDimension;
    end
end
B_discrete_matlab = reshape(pinv(regression_design) * regression_target, modeCount, inputDimension);
regression_residual = regression_design * B_discrete_matlab(:) - regression_target;
regression_relative_residual = norm(regression_residual) / norm(regression_target);

% Convert Bd to a continuous-time source-input matrix through the exact
% zero-order-hold integral.  Since A is diagonal this is mode-wise and does
% not require a branch-sensitive matrix logarithm beyond log(mu)/dt.
lambda_continuous = log(eigenvalues_discrete) ./ options.DeltaT;
integration_factor = zeros(modeCount, 1);
nearZero = abs(lambda_continuous) < 1e-12;
integration_factor(nearZero) = options.DeltaT;
integration_factor(~nearZero) = ...
    (exp(lambda_continuous(~nearZero) .* options.DeltaT) - 1) ./ lambda_continuous(~nearZero);
B_continuous_source = B_discrete_matlab ./ integration_factor;
B_continuous_canonical = 2 .* B_continuous_source;

% One-step and rolled-output diagnostics on the retained controlled data.
one_step_squared_error = zeros(stateDimension, 1);
rollout_squared_error = zeros(stateDimension, 1);
sampleCount = acceptedTrajectoryCount * options.ControlTransitionCount;
for trajectoryIndex = 1:acceptedTrajectoryCount
    zRoll = Psi_controlled(:, 1, trajectoryIndex);
    for timeIndex = 1:options.ControlTransitionCount
        zOne = A_discrete * Psi_controlled(:, timeIndex, trajectoryIndex) + ...
            B_discrete_matlab * U_source(1, timeIndex, trajectoryIndex);
        oneError = real(C_official * zOne) - X_controlled(:, timeIndex + 1, trajectoryIndex);
        one_step_squared_error = one_step_squared_error + oneError.^2;
        zRoll = A_discrete * zRoll + B_discrete_matlab * U_source(1, timeIndex, trajectoryIndex);
        rollError = real(C_official * zRoll) - X_controlled(:, timeIndex + 1, trajectoryIndex);
        rollout_squared_error = rollout_squared_error + rollError.^2;
    end
end
one_step_output_rmse = sqrt(one_step_squared_error ./ sampleCount);
rollout_output_rmse = sqrt(rollout_squared_error ./ sampleCount);

% Published open-loop validation signal and IC.
predictionTransitionCount = round(2 / options.DeltaT);
prediction_input_source = (-1) .^ round((1:predictionTransitionCount) ./ 30);
prediction_state_true = zeros(stateDimension, predictionTransitionCount + 1);
prediction_state_true(:, 1) = [-0.9685; -0.2689];
prediction_lifted = zeros(modeCount, predictionTransitionCount + 1);
prediction_lifted(:, 1) = lift_points(interpolants, prediction_state_true(:, 1));
for timeIndex = 1:predictionTransitionCount
    prediction_state_true(:, timeIndex + 1) = source_step_bug_compatible( ...
        prediction_state_true(:, timeIndex), prediction_input_source(timeIndex), options.DeltaT);
    prediction_lifted(:, timeIndex + 1) = A_discrete * prediction_lifted(:, timeIndex) + ...
        B_discrete_matlab * prediction_input_source(timeIndex);
end
prediction_state_rom = real(C_official * prediction_lifted);
prediction_output_rmse = sqrt(mean((prediction_state_rom - prediction_state_true).^2, 2));

% Exact continuous B_psi for these interpolated eigenfunctions, evaluated by
% centered finite differences.  Source convention has dx2/dt += 0.5*u;
% canonical DynSysLab convention has dx2/dt += u.
x_values = linspace(min(X_free(1, :)), max(X_free(1, :)), options.GridPoints);
v_values = linspace(min(X_free(2, :)), max(X_free(2, :)), options.GridPoints);
[grid_x, grid_v] = meshgrid(x_values, v_values);
grid_points_flat = [grid_x(:).'; grid_v(:).'];
velocityStep = max(1e-7, 1e-5 * range(v_values));
plusPoints = grid_points_flat;
minusPoints = grid_points_flat;
plusPoints(2, :) = plusPoints(2, :) + velocityStep;
minusPoints(2, :) = minusPoints(2, :) - velocityStep;
Psi_grid = lift_points(interpolants, grid_points_flat);
dPsi_dvelocity = (lift_points(interpolants, plusPoints) - ...
    lift_points(interpolants, minusPoints)) ./ (2 * velocityStep);
grid_exact_Bpsi_source = reshape( ...
    0.5 .* dPsi_dvelocity, modeCount, options.GridPoints, options.GridPoints);
grid_exact_Bpsi_canonical = reshape( ...
    dPsi_dvelocity, modeCount, options.GridPoints, options.GridPoints);
[hull_indices, hull_area] = convhull(X_free(1, :).', X_free(2, :).');
training_hull = X_free(:, hull_indices);
grid_inside_training_hull = inpolygon( ...
    grid_x, grid_v, training_hull(1, :), training_hull(2, :));
constant_error_source = squeeze(vecnorm( ...
    grid_exact_Bpsi_source - reshape(B_continuous_source, modeCount, 1, 1), 2, 1));

source_sha256 = sha256_file(sourceFile);
lambda_sha256 = sha256_file(lambdaFile);
wrapper_sha256 = sha256_file(mfilename("fullpath") + ".m");
completedAt = datetime("now", "TimeZone", "UTC");
metadata = struct( ...
    "schema", "dynsyslab.korda-mezic-duffing-control-candidate/v1", ...
    "candidate_only", true, ...
    "mutates_existing_framework", false, ...
    "source_script", sourceFile, ...
    "source_script_sha256", source_sha256, ...
    "lambda_file", lambdaFile, ...
    "lambda_file_sha256", lambda_sha256, ...
    "wrapper_sha256", wrapper_sha256, ...
    "seed", options.Seed, ...
    "precomputed_lambdas", true, ...
    "source_rk4_k4_uses_k1", true, ...
    "input_convention_source_acceleration_gain", 0.5, ...
    "input_convention_canonical_acceleration_gain", 1.0, ...
    "discrete_time_step", options.DeltaT, ...
    "trajectory_count_requested", options.TrajectoryCount, ...
    "trajectory_count_accepted", acceptedTrajectoryCount, ...
    "trajectory_length", options.TrajectoryLength, ...
    "control_transition_count", options.ControlTransitionCount, ...
    "grid_points_per_axis", options.GridPoints, ...
    "started_at_utc", char(startedAt), ...
    "completed_at_utc", char(completedAt));

saveArguments = { ...
    "metadata", "A_discrete", "A_continuous", "C_official", "C_projection", ...
    "eigenvalues_discrete", "lambda_continuous", "integration_factor", ...
    "free_initial_conditions", "free_trajectories", "X_free", ...
    "Psi_free_regression", "Psi_free_evaluated", "boundary_values", ...
    "free_fit_residual_percent", "free_reconstruction_rmse_official", ...
    "free_reconstruction_rmse_projection", "accepted_mask", "accepted_indices", ...
    "X_controlled", "U_source", "U_canonical", "Psi_controlled", ...
    "lifted_derivatives_controlled", "controlled_exact_Bpsi_source", ...
    "controlled_exact_Bpsi_canonical", "positionStepControlled", ...
    "velocityStepControlled", ...
    "regression_target", "B_discrete_matlab", "B_continuous_source", ...
    "B_continuous_canonical", "regression_relative_residual", ...
    "one_step_output_rmse", "rollout_output_rmse", ...
    "prediction_input_source", "prediction_state_true", "prediction_lifted", ...
    "prediction_state_rom", "prediction_output_rmse", ...
    "x_values", "v_values", "grid_x", "grid_v", "Psi_grid", ...
    "grid_exact_Bpsi_source", "grid_exact_Bpsi_canonical", ...
    "constant_error_source", "grid_inside_training_hull", "training_hull", ...
    "hull_area", "velocityStep"};
if options.SaveRegressionDesign
    saveArguments = [saveArguments, {"regression_design"}]; %#ok<AGROW>
end
save(outputFile, saveArguments{:}, "-v7");

summary = struct( ...
    "schema", metadata.schema, ...
    "output_file", outputFile, ...
    "candidate_only", true, ...
    "trajectory_count_accepted", acceptedTrajectoryCount, ...
    "free_fit_residual_percent", free_fit_residual_percent.', ...
    "free_reconstruction_rmse_official", free_reconstruction_rmse_official.', ...
    "free_reconstruction_rmse_projection", free_reconstruction_rmse_projection.', ...
    "regression_relative_residual", regression_relative_residual, ...
    "one_step_output_rmse", one_step_output_rmse.', ...
    "rollout_output_rmse", rollout_output_rmse.', ...
    "prediction_output_rmse", prediction_output_rmse.', ...
    "source_script_sha256", source_sha256, ...
    "lambda_file_sha256", lambda_sha256, ...
    "wrapper_sha256", wrapper_sha256);
jsonFile = replace(outputFile, ".mat", ".summary.json");
fid = fopen(jsonFile, "w");
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, "%s\n", jsonencode(summary, PrettyPrint=true));
clear cleanup
fprintf("Exported candidate parity bundle: %s\n", outputFile);
fprintf("Summary: %s\n", jsonFile);
end


function derivative = source_rhs(state, input)
% Algebraic simplification of Duffing.m's f_u anonymous function.
derivative = [state(2, :); -0.5 .* state(2, :) + state(1, :) - ...
    4 .* state(1, :).^3 + 0.5 .* input];
end


function next = source_step_bug_compatible(state, input, deltaT)
% Reproduce Duffing.m exactly: its fourth RK stage advances from k1.
k1 = source_rhs(state, input);
k2 = source_rhs(state + k1 .* deltaT ./ 2, input);
k3 = source_rhs(state + k2 .* deltaT ./ 2, input);
k4 = source_rhs(state + k1 .* deltaT, input);
next = state + (deltaT ./ 6) .* (k1 + 2 .* k2 + 2 .* k3 + k4);
end


function lifted = lift_points(interpolants, points)
modeCount = numel(interpolants);
lifted = complex(zeros(modeCount, size(points, 2)));
for modeIndex = 1:modeCount
    lifted(modeIndex, :) = interpolants{modeIndex}(points(1, :).', points(2, :).').';
end
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
