function S=load_assembly_session(row,opts)
R=load(char(row.rasterPath));assert(isfield(R,opts.pcaField)&&isfield(R,opts.activityField),'Required raster fields missing.');
candidatePath=char(row.hdMetricsPath);if ismember('candidatePath',row.Properties.VariableNames),candidatePath=char(row.candidatePath);end
P=assembly_preprocess_candidates(candidatePath,char(row.rasterPath));Fall=P.mappedDeltaFoF;Aall=P.mappedRaster;assert(isequal(size(Fall),size(Aall)),'Mapped fields differ.');
fps=findFps(R);assert(isfinite(fps)&&fps>0,'Imaging fps missing or invalid.');roi=P.selectedOriginalRoiIndices(:);candidateIDs=P.selectedCandidateIDs(:);
n=min(size(Fall,1),round(opts.analysisDurationMin*60*fps));F=Fall(1:n,:);A=Aall(1:n,:);F(~isfinite(F))=0;A(~isfinite(A))=0;assert(all(std(F,0,1)>0),'A selected trace is genuinely constant in the assembly window.');assert(numel(roi)>=2,'Fewer than two highly phase-tuned cells.');
label=opts.analysisLabel;if ~isempty(opts.residualField),error('Residual mode is incompatible with strictly mapped candidate preprocessing.');end
S=struct('pca',F,'activity',A,'candidateRoiIdx',roi,'candidateCellIDs',candidateIDs,'preferredPhaseRad',P.selectedPreferredPhaseRad(:),'selectedCandidateOrder',P.selectedCandidateOrder(:),'candidateRasterMapping',P.mapping,'phaseTuning',P.tuning,'harmonicDiagnostics',P.harmonic,'excludedRoiIdx',P.originalRoiIndices(~P.tuning.isPhaseTuned),'nTotalRois',P.nTotalRois,'nFrames',n,'fps',fps,'durationSec',n/fps,'behavior',extractBehavior(R,n),'analysisLabel',label,'candidateSource','candidate-file order mapped by full-trace correlation');
assert(isequal(S.candidateRoiIdx(:),roi(:))&&isequal(S.candidateCellIDs(:),candidateIDs(:))&&size(S.pca,2)==numel(roi)&&numel(S.preferredPhaseRad)==numel(roi),'Candidate mapping invariant failed.');
end
function fps=findFps(R)
fps=NaN;names={'fps','fps_im','frameRate','framerate','samplingRate'};scopes={R};if isfield(R,'params')&&isstruct(R.params),scopes{end+1}=R.params;end
for s=1:numel(scopes),Q=scopes{s};for i=1:numel(names),if isfield(Q,names{i}),v=double(Q.(names{i}));if isscalar(v)&&v>0&&isfinite(v),fps=v;return,end,end,end,end
end
function B=extractBehavior(R,n)
B=struct('available',false,'ahv',[],'missingReason','AHV not found');names={'AHV','ahv','angularHeadVelocity'};
for i=1:numel(names),if isfield(R,names{i}),x=double(R.(names{i})(:));if numel(x)>=n,B.available=true;B.ahv=x(1:n);B.missingReason='';return,end,end,end
end
