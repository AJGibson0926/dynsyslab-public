function [outputFile, summaryFile] = export_legacy_random_forced_rpe(options)
%EXPORT_LEGACY_RANDOM_FORCED_RPE Flatten a saved legacy Bubble/RPE object.
%
% This removable interoperability wrapper is intentionally external to the
% MATLAB DynSysLab source tree.  It reads one saved object archive and writes
% only plain numeric arrays and structs to an output directory outside that
% tree.  The source MAT file and MATLAB suite are never modified.

arguments
    options.SuiteRoot (1,1) string = ...
        fullfile(userpath, "dynsyslab (9-21-24) (working copy)")
    options.InputFile (1,1) string = ...
        fullfile(userpath, "dynsyslab (9-21-24) (working copy)", "dynsyslib", "Bubble", "Bubble.mat")
    options.OutputFile (1,1) string = ""
end

suiteRoot = canonical_path(options.SuiteRoot);
inputFile = canonical_path(options.InputFile);
assert(isfolder(suiteRoot), "LegacyRPEExport:MissingSuite", ...
    "MATLAB suite root does not exist: %s", suiteRoot);
assert(isfile(inputFile), "LegacyRPEExport:MissingInput", ...
    "Saved MATLAB object does not exist: %s", inputFile);

if options.OutputFile == ""
    stamp = string(datetime("now", "Format", "yyyyMMdd_HHmmss"));
    outputFile = fullfile(pwd, "legacy_random_forced_rpe_" + stamp + ".mat");
else
    outputFile = canonical_path(options.OutputFile);
end
[outputDir, outputStem] = fileparts(outputFile);
if isempty(outputDir), outputDir = pwd; end
if ~isfolder(outputDir), mkdir(outputDir); end
summaryFile = fullfile(outputDir, outputStem + ".json");

assert(~is_path_inside(outputFile, suiteRoot), "LegacyRPEExport:UnsafeOutput", ...
    "OutputFile must be outside SuiteRoot so the MATLAB suite remains read-only.");
assert(~strcmpi(outputFile, inputFile), "LegacyRPEExport:OverwriteInput", ...
    "OutputFile must not overwrite the source MAT file.");

originalPath = path;
pathCleanup = onCleanup(@() path(originalPath)); %#ok<NASGU>
addpath(genpath(suiteRoot), "-begin");

loaded = load(inputFile, "dyn");
assert(isfield(loaded, "dyn"), "LegacyRPEExport:MissingDyn", ...
    "The source MAT file has no top-level 'dyn' object.");
dyn = loaded.dyn;
assert(isa(dyn, "Bubble") || isa(dyn, "RPE"), ...
    "LegacyRPEExport:UnexpectedDynamics", ...
    "Expected a Bubble or RPE object, found %s.", class(dyn));
assert(~isempty(dyn.parent), "LegacyRPEExport:MissingSystem", ...
    "The saved dynamics object has no parent system.");
sys = dyn.parent;

clusterCells = as_cells(sys.clst);
trajectoryRecords = struct([]);
recordIndex = 0;
for clusterIndex = 1:numel(clusterCells)
    cluster = clusterCells{clusterIndex};
    instanceCells = as_cells(cluster.inst);
    for trajectoryIndex = 1:numel(instanceCells)
        instance = instanceCells{trajectoryIndex};
        if isempty(instance.sols), continue; end
        if ~isempty(instance.sols.noise)
            error("LegacyRPEExport:NoisyGetter", ...
                ["Trajectory %d/%d has a nonempty time_series.noise setting. " ...
                 "Reading .ts would generate new noise rather than recover stored data."], ...
                clusterIndex, trajectoryIndex);
        end

        stateTime = double(instance.sols.t(:));
        stateValues = orient_samples(double(instance.sols.ts), numel(stateTime), "state");
        if isempty(stateTime) || size(stateValues,1) ~= numel(stateTime), continue; end

        recordIndex = recordIndex + 1;
        trajectoryRecords(recordIndex).cluster_index_one_based = clusterIndex; %#ok<AGROW>
        trajectoryRecords(recordIndex).trajectory_index_one_based = trajectoryIndex;
        trajectoryRecords(recordIndex).initial_state_configured = double(instance.inits(:).');
        trajectoryRecords(recordIndex).state_time = stateTime;
        trajectoryRecords(recordIndex).states = stateValues;
        trajectoryRecords(recordIndex).state_dimension = size(stateValues,2);
        trajectoryRecords(recordIndex).state_sample_count = size(stateValues,1);

        [controlTime, controlValues] = extract_history(instance.ctrl_hist);
        trajectoryRecords(recordIndex).control_call_time = controlTime;
        trajectoryRecords(recordIndex).control_call_values = controlValues;
        trajectoryRecords(recordIndex).control_call_count = size(controlValues,1);
        trajectoryRecords(recordIndex).control_columns = ...
            ["command_u", "filtered_or_applied_u"];
        trajectoryRecords(recordIndex).control_time_has_repeats = ...
            numel(unique(controlTime)) ~= numel(controlTime);
        trajectoryRecords(recordIndex).control_time_is_monotone_nondecreasing = ...
            isempty(controlTime) || all(diff(controlTime) >= 0);
    end
end
assert(~isempty(trajectoryRecords), "LegacyRPEExport:NoTrajectories", ...
    "The saved object contains no nonempty instance-level trajectories.");

controller = sys.CTRL;
legacy = struct;
legacy.schema = "dynsyslab.matlab-legacy-random-forced-rpe/v1";
legacy.source = struct( ...
    "mat_file", inputFile, ...
    "suite_root", suiteRoot, ...
    "top_level_variable", "dyn", ...
    "loaded_class", class(dyn), ...
    "matlab_version", version, ...
    "exported_utc", char(datetime("now", "TimeZone", "UTC", ...
        "Format", "yyyy-MM-dd'T'HH:mm:ss.SSSXXX")));
legacy.classification = struct( ...
    "system_family", "rayleigh-plesset-equation", ...
    "legacy_class", class(dyn), ...
    "model_title", safe_text(dyn.model.title), ...
    "coordinate_system", ternary(logical(dyn.reduced), "dimensionless", "dimensional"), ...
    "forcing_type", safe_text(controller.forcing), ...
    "forcing_subtype", safe_text(controller.subforcing), ...
    "record_kind", "saved randomly/externally forced trajectory ensemble candidate");
legacy.parameters = table_to_struct(dyn.prms);
legacy.dimensionless_groups = table_to_struct(dyn.dimgrps);
legacy.model_parameters = table_to_struct(dyn.model.prms);
legacy.model_dimensionless_groups = table_to_struct(dyn.model.dimgrps);
legacy.system = struct( ...
    "time_resolution", double(sys.tres), ...
    "time_span", double(sys.tspan(:).'), ...
    "time_structure", double(sys.time_struct), ...
    "time_spread", safe_text(sys.time_spread), ...
    "cycle_count", double(sys.cyc_comp), ...
    "solver_options", sys.options);
legacy.controller = struct( ...
    "forcing", safe_text(controller.forcing), ...
    "subforcing", safe_text(controller.subforcing), ...
    "forcing_parameters", double(controller.forcprms), ...
    "active_time", double(controller.time(:).'), ...
    "channel_count", double(controller.chn), ...
    "hard_cap_enabled", logical(controller.hardcap), ...
    "interpolated_history", logical(controller.interp), ...
    "controller_function", safe_function_text(controller.controller));
legacy.trajectories = trajectoryRecords;
legacy.trajectory_count = numel(trajectoryRecords);
legacy.notes = [ ...
    "State arrays are transposed to rows=samples, columns=state components."; ...
    "Controller call records preserve repeated/non-monotone adaptive-solver evaluations."; ...
    "No interpolation or deduplication is performed by this MATLAB wrapper."; ...
    "The second legacy control column is CTRL.B_filt(x,t)*u(end), named conservatively until validated." ...
    ];

save(outputFile, "legacy", "-v7");

summary = make_summary(legacy, outputFile);
write_json(summaryFile, summary);
fprintf("Legacy RPE export complete.\n  MAT:  %s\n  JSON: %s\n", ...
    outputFile, summaryFile);
end


function cells = as_cells(value)
if iscell(value)
    cells = value;
else
    cells = num2cell(value);
end
end


function values = orient_samples(values, sampleCount, label)
if isempty(values)
    values = zeros(sampleCount,0);
elseif size(values,1) == sampleCount
    % Already rows=samples.
elseif size(values,2) == sampleCount
    values = values.';
else
    error("LegacyRPEExport:ShapeMismatch", ...
        "%s array shape %s does not match %d time samples.", ...
        label, mat2str(size(values)), sampleCount);
end
end


function [time, values] = extract_history(historyObject)
if isempty(historyObject)
    time = zeros(0,1);
    values = zeros(0,2);
    return
end
if ~isempty(historyObject.noise)
    error("LegacyRPEExport:NoisyControlGetter", ...
        "A control history has nonempty noise metadata; refusing nondeterministic extraction.");
end
time = double(historyObject.t(:));
values = orient_samples(double(historyObject.ts), numel(time), "control history");
end


function output = table_to_struct(value)
output = struct;
if isempty(value), return; end
assert(istable(value), "LegacyRPEExport:ExpectedTable", ...
    "Expected a parameter table, found %s.", class(value));
assert(height(value) == 1, "LegacyRPEExport:ExpectedScalarTable", ...
    "Expected a one-row parameter table, found %d rows.", height(value));
names = value.Properties.VariableNames;
for index = 1:numel(names)
    item = value.(names{index});
    if iscell(item) && isscalar(item), item = item{1}; end
    if isnumeric(item) || islogical(item)
        output.(names{index}) = double(item);
    else
        output.(names{index}) = safe_text(item);
    end
end
end


function text = safe_text(value)
if isempty(value)
    text = "";
elseif iscell(value)
    text = strjoin(string(value), " / ");
else
    text = strjoin(string(value(:).'), " / ");
end
text = char(text);
end


function text = safe_function_text(value)
if isempty(value)
    text = "";
elseif isa(value, "function_handle")
    text = func2str(value);
else
    text = class(value);
end
end


function output = ternary(condition, whenTrue, whenFalse)
if condition, output = whenTrue; else, output = whenFalse; end
end


function summary = make_summary(legacy, outputFile)
records = legacy.trajectories;
stateCounts = [records.state_sample_count];
controlCounts = [records.control_call_count];
summary = struct;
summary.schema = legacy.schema;
summary.output_file = outputFile;
summary.source = legacy.source;
summary.classification = legacy.classification;
summary.parameters = legacy.parameters;
summary.dimensionless_groups = legacy.dimensionless_groups;
summary.model_parameters = legacy.model_parameters;
summary.model_dimensionless_groups = legacy.model_dimensionless_groups;
summary.system = legacy.system;
summary.controller = legacy.controller;
summary.trajectory_count = legacy.trajectory_count;
summary.state_sample_counts = stateCounts;
summary.control_call_counts = controlCounts;
summary.state_time_spans = arrayfun(@(r) [r.state_time(1), r.state_time(end)], ...
    records, "UniformOutput", false);
summary.control_time_repeat_flags = [records.control_time_has_repeats];
summary.control_time_monotone_flags = [records.control_time_is_monotone_nondecreasing];
summary.notes = legacy.notes;
end


function write_json(filename, value)
encoded = jsonencode(value, "PrettyPrint", true);
fid = fopen(filename, "w");
assert(fid >= 0, "LegacyRPEExport:CannotWriteJSON", ...
    "Could not open JSON output: %s", filename);
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, "%s\n", encoded);
end


function value = canonical_path(pathValue)
value = char(java.io.File(char(pathValue)).getCanonicalPath());
end


function tf = is_path_inside(candidate, root)
candidate = lower(canonical_path(candidate));
root = lower(canonical_path(root));
if ~endsWith(root, filesep), root = [root filesep]; end
tf = startsWith(candidate, root);
end
