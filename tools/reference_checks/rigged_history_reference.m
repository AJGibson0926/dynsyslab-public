function rigged_history_reference(manifestPath, outputDirectory, referenceRoot)
% Reference-only adapter for the pinned, unmodified complex mpEDMD/Rigged DMD.
% Production fitting remains native. Inputs are sampled histories, not equations.
% Manifest: files (.mat G, rows x times); weight_files (.mat W, rows x plans);
% cases (name, source, successor, target, weight_set), angles, density_angles,
% bandwidth, order. Indices are ONE-BASED, target indexes the source dictionary.
% Block QR retains joint weighted products, never normal-equation fitting.
% Sources: Colbrook SINUM (2023), doi:10.1137/22M1521407;
% Colbrook/Drysdale/Horning arXiv:2405.00782; Demmel et al.,
% doi:10.1137/080731992. No upstream source is modified.
q=jsondecode(fileread(manifestPath));
assert(numel(q.files)==numel(q.weight_files),'File and weight block counts differ');
assert(~isfolder(outputDirectory),'Use a new immutable output directory');
addpath(fullfile(referenceRoot,'MColbrook--Rigged-Dynamic-Mode-Decomposition--013f5a429868','main_routines'));
addpath(fullfile(referenceRoot,'MColbrook--Measure-preserving-Extended-Dynamic-Mode-Decomposition--e473dee7b32d'));
plans=max([q.cases.weight_set]);R=cell(plans,1);probes=[];probeNorm=[];samples=0;weightMass=zeros(plans,1);started=tic;
for block=1:numel(q.files)
    data=load(q.files{block},'G');weights=load(q.weight_files{block},'W');G=data.G;W=weights.W;
    assert(size(G,1)==size(W,1)&&size(W,2)>=plans,'Weight dimensions differ');
    assert(all(isfinite(G),'all')&&all(isfinite(W),'all')&&all(W>=0,'all'),'Nonfinite data or negative weights');
    if isempty(probes)
        k=(1:size(G,2)).';probes=[ones(size(k)),sin(.37*k)+1i*cos(.23*k),exp(1i*.71*k)];probeNorm=zeros(plans,size(probes,2));
    end
    values=G*probes;
    for plan=1:plans
        probeNorm(plan,:)=probeNorm(plan,:)+sum(W(:,plan).*abs(values).^2,1);
        weightMass(plan)=weightMass(plan)+sum(W(:,plan));
        R{plan}=qr([R{plan};sqrt(W(:,plan)).*G],0);
    end
    samples=samples+size(G,1);fprintf('Compressed block %d/%d, %d original rows\n',block,numel(q.files),samples);
end
assert(all(weightMass>0),'Positive total weights required');
compressionProbeError=zeros(plans,1);
for plan=1:plans
    actual=sum(abs(R{plan}*probes).^2,1);
    compressionProbeError(plan)=max(abs(actual-probeNorm(plan,:))./max(probeNorm(plan,:),realmin));
    assert(compressionProbeError(plan)<1e-10,'Compression probe failed');
end
mkdir(outputDirectory);
save(fullfile(outputDirectory,'compression.mat'),'R','weightMass','samples','compressionProbeError','-v7');
for j=1:numel(q.cases)
    spec=q.cases(j);assert(~isempty(regexp(spec.name,'^[A-Za-z0-9_-]+$','once')),'Unsafe case identifier');
    px=R{spec.weight_set}(:,spec.source);py=R{spec.weight_set}(:,spec.successor);N=size(px,2);
    assert(isequal(size(px),size(py))&&spec.target>=1&&spec.target<=N,'Case dimensions differ');
    g=zeros(N,1);g(spec.target)=1;singular=svd(px);
    [~,sourceR]=qr(px,0);whitenedY=py/sourceR;
    gramShiftRelative=norm(py'*py-px'*px,'fro')/norm(px'*px,'fro');
    whitenedGramShift=norm(whitenedY'*whitenedY-eye(N),'fro')/sqrt(N);
    fprintf('Fitting %s (%d complex features)\n',spec.name,N);
    [coefficients,density]=riggedDMD(px,py,ones(size(px,1),1),q.bandwidth,q.angles,[],...
        'order',q.order,'g_coeffs',g,'TH2',q.density_angles);
    seconds=toc(started);source=spec.source;successor=spec.successor;angles=q.angles;
    density_angles=q.density_angles;manifest=q;
    save(fullfile(outputDirectory,[spec.name '.mat']),'coefficients','density','source','successor',...
        'angles','density_angles','singular','samples','weightMass','compressionProbeError','seconds','manifest',...
        'gramShiftRelative','whitenedGramShift','-v7');
    fprintf('Saved %s at %.1f s, condition %.4g\n',spec.name,seconds,singular(1)/singular(end));
end
end
