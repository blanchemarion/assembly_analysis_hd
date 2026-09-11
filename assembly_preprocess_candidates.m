function P=assembly_preprocess_candidates(candidatePath,rasterPath)
% Exact candidate preprocessing used by the non-z-scored fitting path.
C=load(candidatePath);
required={'calcium_traces','candidate_cell_ids','network_phase_rad','phi_bin_centers_rad','tuning_curves_phi_all_cells'};
for k=1:numel(required),assert(isfield(C,required{k}),'Candidate file missing %s.',required{k});end
ids=double(C.candidate_cell_ids(:));assert(numel(unique(ids))==numel(ids),'candidate_cell_ids must be unique.');
Y=double(C.calcium_traces);
if size(Y,2)==numel(ids),elseif size(Y,1)==numel(ids),Y=Y';else,error('calcium_traces dimensions do not match candidate_cell_ids.');end
phase=wrapPi(double(C.network_phase_rad(:)));assert(size(Y,1)==numel(phase),'Candidate trace/phase lengths differ.');
R=load(rasterPath,'deltaFoF','raster');assert(isfield(R,'deltaFoF')&&isfield(R,'raster'),'RASTER must contain deltaFoF and raster.');
D=orientRaster(R.deltaFoF,size(Y,1));B=orientRaster(R.raster,size(Y,1));
assert(isequal(size(D),size(B))&&size(D,1)>=size(Y,1),'RASTER fields are incompatible with candidate traces.');
D=D(1:size(Y,1),:);B=B(1:size(Y,1),:);
[roi,assigned,margin]=cachedMapping(candidatePath,rasterPath,Y,ids);
if isempty(roi)
 allR=corr(Y,D,'rows','pairwise');allR(~isfinite(allR))=-Inf;[~,roi]=max(allR,[],2);
 if numel(unique(roi))<numel(roi)&&exist('matchpairs','file')==2
  pairs=matchpairs(-allR,1e6,'min');assert(size(pairs,1)==size(Y,2),'Candidate/RASTER assignment is incomplete.');[~,order]=sort(pairs(:,1));roi=pairs(order,2);
 end
 assert(numel(unique(roi))==numel(roi),'Candidate/RASTER mapping is not one-to-one.');
 assigned=allR(sub2ind(size(allR),(1:size(allR,1))',roi));margin=nan(size(assigned));
 for k=1:numel(assigned),a=allR(k,:);a(roi(k))=-Inf;margin(k)=assigned(k)-max(a);end
end
assert(all(assigned>=0.9999),'Candidate/RASTER full-trace correlation is below 0.9999.');
assert(all(margin>=1e-5),'Candidate/RASTER full-trace correlation margin is below 1e-5.');
mapping=table((1:numel(ids))',ids,roi,assigned,margin,'VariableNames',{'candidateOrder','candidateCellID','originalRasterRoiIndex','fullTraceCorrelation','correlationMargin'});
phaseOk=getPhaseOk(C,numel(phase));tuning=selectTuned(Y,phase,phaseOk,ids);[pref,harmonic]=harmonicPreferredPhase(C,ids);
selected=tuning.selectedCandidateOrder;
assert(isequal(ids(selected),tuning.selectedCandidateIDs)&&isequal(mapping.candidateCellID(selected),tuning.selectedCandidateIDs),'Selection alignment mismatch.');
P=struct('nTotalRois',size(D,2),'mappedDeltaFoF',D(:,roi(selected)),'mappedRaster',B(:,roi(selected)),'candidateIDs',ids, ...
 'originalRoiIndices',roi,'mapping',mapping,'tuning',tuning,'harmonic',harmonic, ...
 'selectedCandidateOrder',selected,'selectedCandidateIDs',ids(selected), ...
 'selectedOriginalRoiIndices',roi(selected),'selectedPreferredPhaseRad',pref(selected));
assert(size(P.mappedDeltaFoF,2)==numel(P.selectedCandidateIDs)&&isequal(mapping.originalRasterRoiIndex(selected),P.selectedOriginalRoiIndices),'Final alignment invariant failed.');
end
function [roi,r,margin]=cachedMapping(candidatePath,rasterPath,Y,ids)
roi=[];r=[];margin=[];assemblyDir=fileparts(mfilename('fullpath'));root=fileparts(assemblyDir);
[recordingDir,candidateName]=fileparts(candidatePath);[morphDir,recording]=fileparts(recordingDir);[~,morph]=fileparts(morphDir);
file=fullfile(root,'data_processed','candidate_raster_mappings',morph,recording,[candidateName '_raster_mapping.mat']);
if ~isfile(file),return,end;Q=load(file,'candidateRasterMapping');if ~isfield(Q,'candidateRasterMapping'),return,end;M=Q.candidateRasterMapping;
ci=dir(candidatePath);ri=dir(rasterPath);valid=isfield(M,'version')&&M.version==2&&strcmp(char(M.candidatePath),candidatePath)&&strcmp(char(M.rasterFile),rasterPath)&&M.candidateBytes==ci.bytes&&M.rasterBytes==ri.bytes&&isequal(M.candidateTraceSize,size(Y))&&isequal(double(M.candidateIds(:)),ids(:));
if ~valid,return,end;roi=double(M.originalRoiIds(:));r=double(M.correlation(:));margin=double(M.margin(:));assert(numel(unique(roi))==numel(ids),'Cached mapping is not one-to-one.');fprintf('  Reused candidate/RASTER mapping: %s\n',file);
end
function X=orientRaster(X,n),X=double(X);assert(ismatrix(X));if size(X,1)<n&&size(X,2)>=n,X=X';end,end
function ok=getPhaseOk(C,n)
ok=true(n,1);names={'phase_ok','network_phase_ok','phase_valid_mask','network_phase_valid','manual_phase_ok'};
for k=1:numel(names),if isfield(C,names{k})&&numel(C.(names{k}))==n,v=C.(names{k});ok=logical(v(:));ok(~isfinite(double(v(:))))=false;return,end,end
if isfield(C,'hdMetrics')&&isstruct(C.hdMetrics)&&isfield(C.hdMetrics,'r1pi')&&isfield(C.hdMetrics.r1pi,'phase')
 for tag={'r1pi','all'},s=tag{1};if isfield(C.hdMetrics.r1pi.phase,s)&&isfield(C.hdMetrics.r1pi.phase.(s),'phase_ok')&&numel(C.hdMetrics.r1pi.phase.(s).phase_ok)==n,v=C.hdMetrics.r1pi.phase.(s).phase_ok;ok=logical(v(:));ok(~isfinite(double(v(:))))=false;return,end,end
end
end
function T=selectTuned(Y,phase,phaseOk,ids)
old=rng;cleanup=onCleanup(@()rng(old));rng(1,'twister'); %#ok<NASGU>
Z=zscore(double(Y),0,1);Z(~isfinite(Z))=0;lo=prctile(Z,2,2);hi=prctile(Z,98,2);Z=bsxfun(@min,Z,hi);Z=bsxfun(@max,Z,lo);Z=bsxfun(@minus,Z,mean(Z,2));
valid=isfinite(phase)&logical(phaseOk);phi=phase(valid);A=Z(valid,:);assert(numel(phi)>=20,'Fewer than 20 valid phase frames.');
nb=36;edges=linspace(-pi,pi,nb+1);centers=((edges(1:end-1)+edges(2:end))/2)';bins=discretize(phi,edges);vb=~isnan(bins);
counts=accumarray(bins(vb),1,[nb 1],@sum,0);occ=counts/max(sum(counts),1);nc=size(A,2);info=nan(nc,1);pref=nan(nc,1);strength=nan(nc,1);rates=nan(nb,nc);
for k=1:nc,rb=accumarray(bins(vb),max(A(vb,k),0),[nb 1],@mean,0);rb=smooth3(rb);[info(k),pref(k),strength(k)]=infoBins(rb,occ,centers);rates(:,k)=rb;end
null=nan(nc,500);minShift=max(1,floor(.1*numel(phi)));maxShift=max(1,floor(.9*numel(phi)));assert(maxShift>minShift,'Invalid circular-shift range.');
for k=1:nc,a=max(A(:,k),0);for s=1:500,sh=circshift(a,randi([minShift maxShift]));rb=accumarray(bins(vb),sh(vb),[nb 1],@mean,0);null(k,s)=infoBins(smooth3(rb),occ,centers);end,end
p=mean(null>=info,2,'omitnan');try,q=mafdr(p,'BHFDR',true);catch,q=bh(p);end;selected=isfinite(q)&q<.05;
T=struct('method','ORI_V15_STEP17_information_circular_shift_BH_FDR','candidateIDs',ids,'infoBits',info,'prefRad',pref,'vectorStrength',strength,'rateBins',rates,'occupancy',occ,'pShuffle',p,'qFDR',q,'null95',prctile(null,95,2),'isPhaseTuned',selected,'validFrameMask',valid,'nValidPhaseFrames',nnz(valid),'nCandidates',nc,'nSelected',nnz(selected),'selectedCandidateOrder',find(selected),'selectedCandidateIDs',ids(selected),'shuffleInfo',null,'params',struct('nPhaseBins',36,'smoothBins',3,'nShuffles',500,'minShiftFraction',.1,'maxShiftFraction',.9,'fdrAlpha',.05,'percentileLow',2,'percentileHigh',98,'randomSeed',1));
end
function y=smooth3(x),y=(circshift(x,-1)+x+circshift(x,1))/3;end
function [info,pref,strength]=infoBins(rate,occ,centers)
mu=sum(occ.*rate);valid=occ>0&rate>0&isfinite(rate);if mu>0,ratio=rate(valid)/mu;info=sum(occ(valid).*ratio.*log2(ratio));else,info=0;end
v=sum(occ.*rate.*exp(1i*centers));pref=angle(v);strength=abs(v)/max(sum(occ.*rate),eps);
end
function q=bh(p),q=nan(size(p));valid=isfinite(p);pv=p(valid);[sp,o]=sort(pv);m=numel(sp);if m==0,return,end;sq=sp.*m./(1:m)';for k=m-1:-1:1,sq(k)=min(sq(k),sq(k+1));end;u=nan(size(pv));u(o)=min(sq,1);q(valid)=u;end
function [pref,H]=harmonicPreferredPhase(C,ids)
phi=double(C.phi_bin_centers_rad(:));Q=double(C.tuning_curves_phi_all_cells);n=numel(ids);
if size(Q,1)==numel(phi)&&size(Q,2)==n,elseif size(Q,2)==numel(phi)&&size(Q,1)==n,Q=Q';else,error('Stored tuning curves do not align with candidates.');end
X=[ones(size(phi)),cos(phi),sin(phi)];pref=nan(n,1);coef=nan(n,3);
for k=1:n,v=isfinite(phi)&isfinite(Q(:,k));if nnz(v)>=3,b=X(v,:)\Q(v,k);coef(k,:)=b';pref(k)=wrapPi(atan2(b(3),b(2)));end,end
H=struct('candidateIDs',ids,'intercept',coef(:,1),'cosCoeff',coef(:,2),'sinCoeff',coef(:,3),'preferredPhaseRad',pref,'method','first circular harmonic of stored tuning_curves_phi_all_cells');
end
function x=wrapPi(x),x=atan2(sin(x),cos(x));x(x<=-pi)=pi;end
