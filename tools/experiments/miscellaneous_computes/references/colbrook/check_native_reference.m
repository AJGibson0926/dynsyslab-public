% Run only the pinned mpEDMD and Rigged DMD routines on a synthetic rotation.
% The native parity test compares operators, densities and wave packets.
root = fileparts(fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))))));
refs = fullfile(root,'build','reference_sources','colbrook');
savedPath = path;
cleanupPath = onCleanup(@() path(savedPath));
addpath(fullfile(refs,'MColbrook--Measure-preserving-Extended-Dynamic-Mode-Decomposition--e473dee7b32d'));
addpath(fullfile(refs,'MColbrook--Rigged-Dynamic-Mode-Decomposition--013f5a429868','main_routines'));
theta = (0:127)'*2*pi/128;
PX = [cos(theta),sin(theta)];
PY = [cos(theta+0.31),sin(theta+0.31)];
W = ones(128,1)/128;
TH = linspace(-pi,pi,257);
[mpK,~,~] = mpEDMDqr(PX,PY,W);
[g1,xi1] = riggedDMD(PX,PY,W,0.08,TH,PX(:,1),'order',1,'PoleType','equi','g_coeffs',[1;0],'TH2',TH);
[g4,xi4] = riggedDMD(PX,PY,W,0.08,TH,PX(:,1),'order',4,'PoleType','equi','g_coeffs',[1;0],'TH2',TH);
save(fullfile(root,'build','colbrook-reference-parity.mat'),'PX','PY','W','TH','mpK','g1','xi1','g4','xi4','-v7');
disp('Saved pinned upstream mpEDMD / Rigged DMD numerical reference.');
