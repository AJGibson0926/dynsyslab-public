function metrics = compare_rigged_reference(inputFile, outputFile, mpedmdFolder, riggedFolder)
%COMPARE_RIGGED_REFERENCE Compare a native fit with pinned upstream algorithms.
% Test adapter only; not a production fitter. Input MAT: PX, PY, W, B,
% g_coeffs, TH, TH2, native_action, native_packets, optional epsilon/order.
% PX and PY already use the same frozen input whitening basis B. Upstream
% mpEDMD recomputes its own operator. Output files are immutable.
% Sources: Colbrook, mpEDMD (2023); Colbrook, Drysdale & Horning,
% Rigged DMD, https://arxiv.org/abs/2405.00782.
if isfile(outputFile)
    error('DynSysLab:ExistingReference', 'Refusing to overwrite reference output');
end
previousPath = path;
cleanupPath = onCleanup(@() path(previousPath)); %#ok<NASGU>
addpath(mpedmdFolder);
addpath(riggedFolder);
D = load(inputFile);
epsilon = .05;
order = 6;
if isfield(D, 'epsilon'), epsilon = double(D.epsilon); end
if isfield(D, 'order'), order = double(D.order); end
started = tic;
[gTH, xi] = riggedDMD(D.PX, D.PY, D.W, epsilon, D.TH, [], ...
    'order', order, 'PoleType', 'equi', 'g_coeffs', D.g_coeffs, 'TH2', D.TH2);
coefficients = D.B * gTH;
[mpK, ~, ~] = mpEDMDqr(D.PX, D.PY, D.W);
metrics.operator_relative_error = norm(mpK-D.native_action, 'fro') / norm(D.native_action, 'fro');
metrics.packet_relative_error = norm(D.PX*(gTH-D.B\D.native_packets), 'fro') / norm(D.PX*gTH, 'fro');
% Audit the supplied whitening independently; rewhitening can amplify tiny
% singular directions even when both operator implementations are correct.
gram = D.PX' * (D.W .* D.PX);
metrics.whitening_relative_frobenius_error = norm(gram-eye(size(gram)), 'fro') / sqrt(size(gram,1));
gramEigenvalues = eig((gram+gram')/2);
metrics.whitening_minimum_eigenvalue = min(real(gramEigenvalues));
metrics.whitening_maximum_eigenvalue = max(real(gramEigenvalues));
metrics.seconds = toc(started);
metrics.reference_files = {which('riggedDMD'), which('mpEDMDqr')};
save(outputFile, 'coefficients', 'xi', 'metrics', '-v7');
end
