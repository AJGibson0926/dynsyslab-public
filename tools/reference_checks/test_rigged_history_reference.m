function test_rigged_history_reference(referenceRoot, workDirectory)
% Synthetic reference-adapter regression; writes only in a new caller directory.
assert(~isfolder(workDirectory),'Use a fresh verification directory');mkdir(workDirectory);
previousRandomState=rng;restoreRandomState=onCleanup(@()rng(previousRandomState)); %#ok<NASGU>
rng(712);omega=linspace(-pi,pi,96).';phase=rand(96,1)*2*pi;
G=exp(1i*(phase+omega*(0:16)));W=.2+rand(96,1);
files=cell(2,1);weight_files=cell(2,1);
for j=1:2
    rows=(j-1)*48+(1:48);part=struct('G',G(rows,:));weights=struct('W',W(rows));
    files{j}=fullfile(workDirectory,sprintf('data%d.mat',j));weight_files{j}=fullfile(workDirectory,sprintf('weights%d.mat',j));
    save(files{j},'-struct','part');save(weight_files{j},'-struct','weights');
end
q=struct('files',{files},'weight_files',{weight_files},'cases',struct('name','synthetic','source',1:16,'successor',2:17,'target',9,'weight_set',1),...
    'angles',[.3,.7],'density_angles',linspace(-pi,pi,30),'bandwidth',.08,'order',4);
manifest=fullfile(workDirectory,'manifest.json');f=fopen(manifest,'w');fprintf(f,'%s',jsonencode(q));fclose(f);
rigged_history_reference(manifest,fullfile(workDirectory,'output'),referenceRoot);
g=zeros(16,1);g(9)=1;
[coeff,density]=riggedDMD(G(:,1:16),G(:,2:17),W,.08,q.angles,[],'order',4,'g_coeffs',g,'TH2',q.density_angles);
saved=load(fullfile(workDirectory,'output','synthetic.mat'));
fieldError=norm(G(:,1:16)*(saved.coefficients-coeff),'fro')/norm(G(:,1:16)*coeff,'fro');
densityError=norm(saved.density-density)/norm(density);
assert(fieldError<1e-10&&densityError<1e-10,'Compressed reference differs from raw reference');
fprintf('Reference adapter passed: field %.3g, density %.3g\n',fieldError,densityError);
end
