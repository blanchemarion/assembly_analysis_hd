function report = compare_assembly_triggered_averages_across_morphs(batchResults,varargin)
%COMPARE_ASSEMBLY_TRIGGERED_AVERAGES_ACROSS_MORPHS Onset-aligned assembly signals.
%
% report = compare_assembly_triggered_averages_across_morphs(results)
% report = compare_assembly_triggered_averages_across_morphs(results, ...
%   'controlRandomSeed',2)
% Default offsets are -2:4; baseline uses negative offsets. Controls use uniform
% unique pseudo-onsets; windows may overlap and other assemblies may be active.
% Controls sample up to controlWindowCount (default 10000) per assembly,
% independently of retained onsets. No signal values guide selection.
% Dashed references use all outside frames exactly as the bar helper does.
% Fully inactive windows exclude short gaps and change frame weighting.
%
% Windows are averaged within assembly, assemblies within fish, then fish
% receive equal weight in morph means and SEMs. By default, onset windows
% exclude prior activation in the prewindow and subsequent onsets in the postwindow.
% Activation classification uses the same shared defaults and computation as
% compare_assembly_properties_across_morphs. Pass identical activation name/value
% overrides to both scripts when comparing nondefault runs. Complete-window
% filtering affects averaging only; episodeSelection retains all detected onsets.

if nargin<1||isempty(batchResults)
 scriptDir=fileparts(mfilename('fullpath'));
 savedBatch=fullfile(scriptDir,'data_processed','assemblies','session_manifest.mat');
 assert(isfile(savedBatch),'Saved assembly manifest not found: %s',savedBatch);
 saved=load(savedBatch,'manifest','opts');
 assert(isfield(saved,'manifest')&&isfield(saved,'opts'), ...
  'Saved manifest must contain manifest and opts.');
 batchResults=struct('manifest',saved.manifest,'options',saved.opts);
end

p=inputParser;
addRequired(p,'batchResults',@(x)isstruct(x)&&isfield(x,'manifest')&&isfield(x,'options'));
addParameter(p,'outputDir','',@(x)ischar(x)||isstring(x));
addParameter(p,'relativeFrames',-2:50,@(x)isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(diff(x)==1)&&x(1)<0&&x(end)>0&&all(x==round(x)));
addParameter(p,'isolatedWindows',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'controlRandomSeed',2,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0&&x<=2^32-1&&x==round(x));
addParameter(p,'controlWindowCount',10000,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1&&x==round(x));
activationDefaults=assembly_activation_defaults();
addParameter(p,'minActiveMembers',activationDefaults.minActiveMembers,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&x==round(x));
addParameter(p,'onRecruitmentFraction',activationDefaults.onRecruitmentFraction,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p,'offRecruitmentFraction',activationDefaults.offRecruitmentFraction,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p,'minMatchingIndex',activationDefaults.minMatchingIndex,@(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
addParameter(p,'activationNullAlpha',activationDefaults.activationNullAlpha,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<=1);
addParameter(p,'nNullShuffles',activationDefaults.nNullShuffles,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&x==round(x));
addParameter(p,'minDurationSeconds',activationDefaults.minDurationSeconds,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(p,'minNullShiftSeconds',activationDefaults.minNullShiftSeconds,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=0);
addParameter(p,'activationRandomSeed',activationDefaults.randomSeed,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x==round(x));
addParameter(p,'figureVisible','off',@(x)ischar(x)||isstring(x));
addParameter(p,'pngResolution',300,@(x)isnumeric(x)&&isscalar(x)&&x>=150&&x==round(x));
parse(p,batchResults,varargin{:}); q=p.Results;

opts=batchResults.options;
if strlength(string(q.outputDir))==0
 outputName=['assembly_triggered_averages_' char(opts.analysisLabel)];
 tag=sprintf('_frames%g_to%g_K%d_Ron%g_Roff%g_MI%g_a%g_n%d_dur%g_safe%g_seed%d_controlSeed%d', ...
  q.relativeFrames(1),q.relativeFrames(end),q.minActiveMembers,q.onRecruitmentFraction, ...
  q.offRecruitmentFraction,q.minMatchingIndex,q.activationNullAlpha, ...
  q.nNullShuffles,q.minDurationSeconds,q.minNullShiftSeconds,q.activationRandomSeed,q.controlRandomSeed);
 outputName=[outputName strrep(tag,'.','p') sprintf('_controls%d_isolated%d_fishWeighted',q.controlWindowCount,q.isolatedWindows)];
 outputDir=fullfile(opts.figureRoot,outputName);
else
 outputDir=char(q.outputDir);
end
if ~isfolder(outputDir),mkdir(outputDir);end

M=batchResults.manifest;
assert(all(ismember({'status','outputPath','morph','recording'},M.Properties.VariableNames)), ...
 'Manifest must contain status, outputPath, morph, and recording.');
use=M.status=="processed"|M.status=="existing";
M=M(use,:);
if isempty(M),error('assembly:NoBatchOutputs', ...
  'No processed/existing recordings in the batch manifest.');end

activationOptions=struct('minActiveMembers',q.minActiveMembers, ...
 'onRecruitmentFraction',q.onRecruitmentFraction, ...
 'offRecruitmentFraction',q.offRecruitmentFraction, ...
 'minMatchingIndex',q.minMatchingIndex,'activationNullAlpha',q.activationNullAlpha, ...
 'nNullShuffles',q.nNullShuffles,'minDurationSeconds',q.minDurationSeconds, ...
 'minNullShiftSeconds',q.minNullShiftSeconds,'randomSeed',q.activationRandomSeed);
assemblyRecords=struct('morph',{},'fish',{},'assemblyID',{},'file',{}, ...
 'fps',{},'nValidEpisodes',{},'relativeFrames',{},'assemblyDeltaFoF',{}, ...
 'activeMemberDeltaFoF',{},'inactiveMemberDeltaFoF',{}, ...
 'episodeBaselineDeltaFoF',{},'fractionMembersActive',{}, ...
 'baselineCorrectedFractionMembersActive',{});
controlRecords=assemblyRecords;
controlSelection=struct('morph',{},'fish',{},'file',{},'assemblyID',{}, ...
 'requestedWindowCount',{},'obtainedWindowCount',{},'eligiblePseudoOnsetCount',{}, ...
 'shortfall',{},'seed',{},'sampledPseudoOnsetIndices',{},'exclusionReason',{});
controlStream=RandStream('mt19937ar','Seed',q.controlRandomSeed);
relativeFrames=q.relativeFrames(:);
activationAudits=struct('morph',{},'fish',{},'file',{},'audit',{});
episodeSelection=struct('morph',{},'fish',{},'file',{},'assemblyID',{}, ...
 'nDetectedEpisodes',{},'nRetainedEpisodes',{},'nExcludedBoundaryEpisodes',{}, ...
 'detectedOnsetIndices',{},'retainedOnsetIndices',{},'excludedBoundaryOnsetIndices',{}, ...
 'nExcludedOverlapEpisodes',{},'excludedOverlapOnsetIndices',{});
exclusions=struct('morph',{},'fish',{},'file',{},'assemblyID',{}, ...
 'category',{},'reason',{});
runLengthValues=cell(3,1);runLengthMorphs=["surface","molino","pachon"]; 

fprintf('\n=== Assembly-triggered averages: %d files ===\n',height(M));
for i=1:height(M)
 file=char(M.outputPath(i));morph=lower(string(M.morph(i)));fish=string(M.recording(i));
 fprintf('[%d/%d] %s/%s\n',i,height(M),morph,fish);
 try
  if ~isfile(file),error('assembly:MissingFile','Output file missing: %s',file);end
  X=load(file,'model','S');
  requireFields(X,{'model','S'},'saved output');
  requireFields(X.model,{'membership'},'model');
  requireFields(X.S,{'pca','activity','fps'},'S');
  membership=logical(X.model.membership);
  requireFields(X.model,{'candidateRoiIdx'},'model');
  requireFields(X.S,{'candidateRoiIdx'},'S');
  modelCandidateIDs=double(X.model.candidateRoiIdx(:));
  sessionCandidateIDs=double(X.S.candidateRoiIdx(:));
  assert(isequal(modelCandidateIDs,sessionCandidateIDs)&&size(membership,1)==numel(modelCandidateIDs), ...
   'Assembly membership rows do not match the detected candidate-neuron population.');
  deltaFoFTraces=double(X.S.pca);
  activityRaster=double(X.S.activity);fps=double(X.S.fps);
  assert(ismatrix(membership)&&ismatrix(deltaFoFTraces)&&ismatrix(activityRaster), ...
   'model.membership, S.pca, and S.activity must be two-dimensional.');
  assert(isequal(size(deltaFoFTraces),size(activityRaster)), ...
   'S.pca and S.activity must have identical [frames x neurons] dimensions.');
  assert(size(activityRaster,2)==size(membership,1), ...
   'S activity columns must align one-to-one with model.membership rows.');
  assert(isscalar(fps)&&isfinite(fps)&&fps>0,'S.fps must be a finite positive scalar.');
  nAssemblies=size(membership,2);nFrames=size(activityRaster,1);
  if nAssemblies==0
   exclusions=addExclusion(exclusions,morph,fish,file,NaN,"zero assemblies", ...
    "zero assemblies detected");
   fprintf('  Excluded from comparison: zero assemblies detected.\n');
   continue
  end
  assert(isfield(X.model,'candidateCellIDs')&&isfield(X.S,'candidateCellIDs')&& ...
   isequal(double(X.model.candidateCellIDs(:)),double(X.S.candidateCellIDs(:)))&& ...
   numel(X.S.candidateCellIDs)==size(membership,1),'Candidate IDs are not aligned with membership rows.');
  % Apply the recording-level exclusion before contributing to any output,
  % including the neuron run-length comparison.
  runMorph=find(runLengthMorphs==morph,1);
  assert(~isempty(runMorph),'Unknown morph for run-length audit: %s',morph);
  runLengthValues{runMorph}=[runLengthValues{runMorph};binaryRunLengths(activityRaster>0)];
  currentActivation=compute_assembly_activations(membership,activityRaster,fps,activationOptions);
  requireFields(currentActivation,{'active'},'activation result');
  active=logical(currentActivation.active);
  activationAudits(end+1)=struct('morph',morph,'fish',fish, ... %#ok<AGROW>
   'file',string(file),'audit',currentActivation);
  assert(isequal(size(active),[nAssemblies nFrames]), ...
   'Activation matrix must be [assemblies x frames] and align with model.membership and S.activity.');
  % Match the properties comparison: rising edges of duration-filtered activity.
  onsets=false(size(active));
  if nFrames>0,onsets(:,1)=active(:,1);end
  if nFrames>1,onsets(:,2:end)=active(:,2:end)&~active(:,1:end-1);end
  neuronActive=activityRaster>0; % Exact definition used by apply_frozen_assembly_template.
  assert(isequal(size(neuronActive),[nFrames size(membership,1)]), ...
   'Neuron active-state matrix must align with S.activity and model.membership.');
  for a=1:nAssemblies
   memberIdx=find(membership(:,a));
   if isempty(memberIdx)
    exclusions=addExclusion(exclusions,morph,fish,file,a,"zero members", ...
     "assembly membership contains no members");
    continue
   end
   fractionSignal=mean(neuronActive(:,memberIdx),2);
   assert(numel(fractionSignal)==nFrames, ...
    'Framewise assembly signals do not align with the activation matrix.');
   detectedOnset=find(onsets(a,:));
   valid=detectedOnset+relativeFrames(1)>=1&detectedOnset+relativeFrames(end)<=nFrames;
   boundaryValid=valid;
   overlap=false(size(detectedOnset));
   if q.isolatedWindows
    for e=find(boundaryValid)
     t=detectedOnset(e);
     overlap(e)=any(active(a,t+relativeFrames(relativeFrames<0))) || ...
      any(onsets(a,t+relativeFrames(relativeFrames>0)));
    end
   end
   valid=boundaryValid & ~overlap;
   onset=detectedOnset(valid);
   episodeSelection(end+1)=struct('morph',morph,'fish',fish,'file',string(file), ... %#ok<AGROW>
    'assemblyID',a,'nDetectedEpisodes',numel(detectedOnset), ...
    'nRetainedEpisodes',nnz(valid),'nExcludedBoundaryEpisodes',nnz(~boundaryValid), ...
    'detectedOnsetIndices',detectedOnset,'retainedOnsetIndices',onset, ...
    'excludedBoundaryOnsetIndices',detectedOnset(~boundaryValid), ...
    'nExcludedOverlapEpisodes',nnz(overlap),'excludedOverlapOnsetIndices',detectedOnset(overlap));
   if isempty(onset)
    exclusions=addExclusion(exclusions,morph,fish,file,a,"no eligible onset windows", ...
     "no activation onset passes boundary and isolated-window selection");
   else
    r=windowRecord(deltaFoFTraces,neuronActive,memberIdx,onset,relativeFrames,morph,fish,a,file,fps);
    assemblyRecords(end+1)=r; %#ok<AGROW>
   end
   [pseudoOnsets,eligibleCount]=sampleControls(active(a,:),relativeFrames,q.controlWindowCount,controlStream);
   obtained=numel(pseudoOnsets);reason="";
   if obtained==0,reason="no eligible non-activation windows";end
   controlSelection(end+1)=struct('morph',morph,'fish',fish,'file',string(file), ...
    'assemblyID',a,'requestedWindowCount',q.controlWindowCount,'obtainedWindowCount',obtained, ...
    'eligiblePseudoOnsetCount',eligibleCount,'shortfall',q.controlWindowCount-obtained, ...
    'seed',q.controlRandomSeed,'sampledPseudoOnsetIndices',pseudoOnsets,'exclusionReason',reason); %#ok<AGROW>
   if obtained>0
    controlRecords(end+1)=windowRecord(deltaFoFTraces,neuronActive,memberIdx,pseudoOnsets,relativeFrames,morph,fish,a,file,fps); %#ok<AGROW>
   end
  end
 catch ME
  warning('assembly:TriggeredAverageFailed','%s/%s failed: %s',morph,fish,ME.message);
  category="recording failure";
  if contains(lower(string(ME.message)),"missing")||contains(lower(string(ME.message)),"must contain")
   category="missing fields";
  end
  exclusions=addExclusion(exclusions,morph,fish,file,NaN,category,string(ME.message));
 end
end

[longTable,morphTable,morphSummary]=recordTables(assemblyRecords,relativeFrames);
[controlLongTable,controlMorphTable,controlMorphSummary]=recordTables(controlRecords,relativeFrames);
controlLongTable.Properties.VariableNames{'nValidEpisodes'}='nWindows';
controlSelectionTable=struct2table(controlSelection,'AsArray',true);
controlSelectionTable.sampledPseudoOnsetIndices=string(arrayfun(@(r)jsonencode(r.sampledPseudoOnsetIndices),controlSelection,'UniformOutput',false))';
writetable(controlLongTable,fullfile(outputDir,'assembly_triggered_averages_outside_activation_long.csv'));
writetable(controlMorphTable,fullfile(outputDir,'assembly_triggered_averages_outside_activation_morph_summary.csv'));
writetable(controlSelectionTable,fullfile(outputDir,'assembly_triggered_averages_outside_activation_window_selection.csv'));

exclusionTable=exclusionStructToTable(exclusions);
writetable(longTable,fullfile(outputDir,'assembly_triggered_averages_long.csv'));
writetable(morphTable,fullfile(outputDir,'assembly_triggered_averages_morph_summary.csv'));
writetable(exclusionTable,fullfile(outputDir,'assembly_triggered_averages_exclusions.csv'));
selectionTable=struct2table(episodeSelection);
selectionTable=selectionTable(:,{'morph','fish','assemblyID','nDetectedEpisodes','nRetainedEpisodes','nExcludedBoundaryEpisodes','nExcludedOverlapEpisodes'});
writetable(selectionTable,fullfile(outputDir,'onset_window_selection.csv'));
[runLengthSummary,runLengthDistribution]=makeRunLengthFigure(runLengthValues,runLengthMorphs,outputDir,q.figureVisible,q.pngResolution);
writetable(runLengthSummary,fullfile(outputDir,'binary_raster_run_length_summary.csv'));
writetable(runLengthDistribution,fullfile(outputDir,'binary_raster_run_length_distribution.csv'));
save(fullfile(outputDir,'assembly_activation_audit.mat'),'activationAudits','activationOptions','episodeSelection','-v7.3');

settings=q;settings.activationOptions=activationOptions;
settings.relativeFrames=relativeFrames;settings.baselineRelativeFrames=relativeFrames(relativeFrames<0);
settings.controlDefinition="uniform unique pseudo-onsets; all selected frames outside the same assembly final active mask; windows may overlap; other assemblies may be active; up to controlWindowCount per assembly independently of retained activation onsets; full-window selection excludes short gaps and changes frame weighting";
settings.controlRandomGenerator="mt19937ar; separate stream consumed in manifest/assembly order";
settings.observationalUnit="fish; windows averaged within assembly, assemblies within fish";
settings.neuronActiveDefinition="S.activity > 0 without temporal extension";
settings.episodeSelectionDefinition="all rising edges detected; complete windows required; isolatedWindows excludes prewindow activity and subsequent postwindow onsets without changing activation detection";
settings.activationDefinition="empirical circular-shift coactivation null with recruitment/matching-index onset criteria and hysteretic offsets";
settings.assemblyDeltaFoFSource="S.pca (mappedDeltaFoF)";
settings.assemblyDeltaFoFDefinition= ...
 "per-episode unweighted arithmetic mean of member-neuron S.pca delta-F-over-F traces, with each window baseline-corrected by subtracting its mean over all negative frame offsets before averaging windows";
settings.fractionActiveDefinition="uncorrected fraction of assembly members with S.activity > 0; windows averaged within assembly, then assemblies averaged within fish and fish equally weighted within morph";
settings.binaryEventsPerNeuronPerFrameDefinition="sum of binary-positive member-neuron frames divided by the number of assembly members at each aligned frame; numerically identical to the uncorrected fraction of assembly members active";
settings.primaryActivityFigureDefinition="left panel: onset-conditioned frames -2 through 2; right panel: unconditioned all-outside-frame morph mean shown over frames 15 through the requested endpoint for visual comparison only; panels use separate y-scales and the right panel is not a continuation of the onset trace";
settings.baselineCorrectedFractionActiveDefinition="member active fraction minus its per-window mean over all negative frame offsets; retained in tables for reference, not plotted";
settings.conditionalDeltaFoFDefinition="per-neuron baseline-corrected delta-F-over-F averaged separately over members with S.activity > 0 and ~(S.activity > 0) at each aligned frame";
settings.runLengthDefinition="contiguous runs of S.activity > 0, pooled over candidate neurons and recordings within morph; no event extension";
settings.runLengthTenFrameCheck="reports exact 10-frame and inclusive 9-to-11-frame percentages";
report=struct('controlPerAssembly',{controlRecords},'controlMorphSummaries',{controlMorphSummary}, ...
 'controlLongTable',controlLongTable,'controlMorphSummaryTable',controlMorphTable, ...
 'controlWindowSelection',{controlSelection},'controlWindowSelectionTable',controlSelectionTable, ...
 'perAssembly',{assemblyRecords},'morphSummaries',{morphSummary}, ...
 'assemblyLongTable',longTable,'morphSummaryTable',morphTable, ...
 'runLengthSummary',runLengthSummary,'runLengthDistribution',runLengthDistribution, ...
 'settings',settings,'activationAudits',{activationAudits},'episodeSelection',{episodeSelection},'exclusions',exclusionTable,'outputDir',string(outputDir));
report.onsetDiagnostics=table;
if ~isempty(assemblyRecords)
 report.onsetDiagnostics=plot_assembly_onset_diagnostics(report,outputDir,q.figureVisible,q.pngResolution);
end
report.outsideActivity=plot_assembly_outside_activity(report,outputDir,q.figureVisible,q.pngResolution);
report.outsideEstimatorComparison=outsideEstimatorComparison(report);
writetable(report.outsideEstimatorComparison,fullfile(outputDir,'outside_activation_estimator_comparison.csv'));
disp(report.outsideEstimatorComparison);
makeFigure(morphSummary,outputDir,q.figureVisible,q.pngResolution,false,report.outsideActivity);
makeFigure(controlMorphSummary,outputDir,q.figureVisible,q.pngResolution,true,report.outsideActivity,report.outsideEstimatorComparison);
save(fullfile(outputDir,'assembly_triggered_averages_report.mat'),'report','-v7.3');
fprintf('Triggered-average comparison complete. Outputs: %s\n',outputDir);
end

function requireFields(S,names,label)
missing=names(~isfield(S,names));
assert(isempty(missing),'%s missing required field(s): %s.',label,strjoin(missing,', '));
end

function E=addExclusion(E,morph,fish,file,assemblyID,category,reason)
r=struct('morph',string(morph),'fish',string(fish),'file',string(file), ...
 'assemblyID',assemblyID,'category',string(category),'reason',string(reason));
if isempty(E),E=r;else,E(end+1)=r;end
end

function T=exclusionStructToTable(E)
if isempty(E)
 T=table(strings(0,1),strings(0,1),strings(0,1),nan(0,1),strings(0,1),strings(0,1), ...
  'VariableNames',{'morph','fish','file','assemblyID','category','reason'});
else
 T=struct2table(E);
end
end

function [longTable,morphTable,summaries]=makeOutputTables(records,time,deltaFoF,activeDeltaFoF,inactiveDeltaFoF,fractionRaw,fractionCorrected)
nRecords=numel(records);nTime=numel(time);
morph=repelem(reshape(string({records.morph}),[],1),nTime,1);
fish=repelem(reshape(string({records.fish}),[],1),nTime,1);
assemblyID=repelem(reshape([records.assemblyID],[],1),nTime,1);
nValidEpisodes=repelem(reshape([records.nValidEpisodes],[],1),nTime,1);
relativeFrame=repmat(time,nRecords,1);
assemblyDeltaFoF=reshape(deltaFoF',[],1);
activeMemberDeltaFoF=reshape(activeDeltaFoF',[],1);
inactiveMemberDeltaFoF=reshape(inactiveDeltaFoF',[],1);
fractionMemberNeuronsActive=reshape(fractionRaw',[],1);
baselineCorrectedFractionActive=reshape(fractionCorrected',[],1);
longTable=table(morph,fish,assemblyID,nValidEpisodes,relativeFrame, ...
 assemblyDeltaFoF,activeMemberDeltaFoF,inactiveMemberDeltaFoF, ...
 fractionMemberNeuronsActive,baselineCorrectedFractionActive);

order=["surface","molino","pachon"];
summaries=struct('morph',{},'nAssemblies',{},'nFish',{},'relativeFrame',{}, ...
 'deltaFoFMean',{},'deltaFoFSD',{},'deltaFoFSEM',{}, ...
 'activeDeltaFoFMean',{},'activeDeltaFoFSD',{},'activeDeltaFoFSEM',{}, ...
 'inactiveDeltaFoFMean',{},'inactiveDeltaFoFSD',{},'inactiveDeltaFoFSEM',{}, ...
 'fractionMean',{},'fractionSD',{},'fractionSEM',{});
rows=cell(numel(order),1);recordMorph=string({records.morph})';
for m=1:numel(order)
 pick=recordMorph==order(m);n=sum(pick);
 fishIDs=string({records.fish})';nFish=numel(unique(fishIDs(pick)));
 [dm,ds,dsem]=fishTraceSummary(deltaFoF,pick,fishIDs);
 [am,as,asem]=fishTraceSummary(activeDeltaFoF,pick,fishIDs);
 [im,is,isem]=fishTraceSummary(inactiveDeltaFoF,pick,fishIDs);
 [fm,fs,fsem]=fishTraceSummary(fractionRaw,pick,fishIDs);
 [fcm,fcs,fcsem]=fishTraceSummary(fractionCorrected,pick,fishIDs);
 summaries(end+1)=struct('morph',order(m),'nAssemblies',n,'nFish',nFish,'relativeFrame',time, ...
  'deltaFoFMean',dm','deltaFoFSD',ds','deltaFoFSEM',dsem', ...
  'activeDeltaFoFMean',am','activeDeltaFoFSD',as','activeDeltaFoFSEM',asem', ...
  'inactiveDeltaFoFMean',im','inactiveDeltaFoFSD',is','inactiveDeltaFoFSEM',isem', ...
  'fractionMean',fm','fractionSD',fs','fractionSEM',fsem'); %#ok<AGROW>
 rows{m}=table(repmat(order(m),nTime,1),time,repmat(n,nTime,1),repmat(nFish,nTime,1), ...
  dm',ds',dsem',am',as',asem',im',is',isem',fm',fs',fsem',fcm',fcs',fcsem', ...
  'VariableNames',{'morph','relativeFrame','nAssemblies','nFish', ...
  'deltaFoFMean','deltaFoFSD','deltaFoFSEM', ...
  'activeDeltaFoFMean','activeDeltaFoFSD','activeDeltaFoFSEM', ...
  'inactiveDeltaFoFMean','inactiveDeltaFoFSD','inactiveDeltaFoFSEM', ...
  'fractionMean','fractionSD','fractionSEM', ...
  'baselineCorrectedFractionMean','baselineCorrectedFractionSD','baselineCorrectedFractionSEM'});
end
morphTable=vertcat(rows{:});
end

function [mu,sd,sem]=fishTraceSummary(X,pick,fishIDs)
ids=unique(fishIDs(pick));Y=nan(numel(ids),size(X,2));
for f=1:numel(ids),Y(f,:)=mean(X(pick & fishIDs==ids(f),:),1,'omitmissing');end
[mu,sd,sem]=traceSummary(Y,true(size(Y,1),1));
end

function [mu,sd,sem]=traceSummary(X,pick)
nTime=size(X,2);Y=X(pick,:);
if isempty(Y),mu=nan(1,nTime);sd=mu;sem=mu;return,end
mu=mean(Y,1,'omitmissing');sd=std(Y,0,1,'omitmissing');
nEffective=sum(isfinite(Y),1);sem=sd./sqrt(nEffective);sem(nEffective==0)=NaN;
sd(nEffective==1)=0;sem(nEffective==1)=0;
end

function makeFigure(S,outDir,visibility,resolution,isControl,outside,comparison)
if nargin<6,outside=[];end
if nargin<7,comparison=[];end
if ~isControl && ~isempty(outside)
 makeStateComparisonFigure(S,outside,outDir,visibility,resolution,false);
 makeStateComparisonFigure(S,outside,outDir,visibility,resolution,true);
 return
end
makeSignalFigure(S,outDir,visibility,resolution,isControl,outside,comparison,false);
makeSignalFigure(S,outDir,visibility,resolution,isControl,outside,comparison,true);
end

function makeStateComparisonFigure(S,outside,outDir,visibility,resolution,showBinaryEventsPerFrame)
colors=assembly_figure_style().morphColors;
displayNames=["Surface","Molino","Pachon"];
fig=figure("Visible",char(visibility),"Color","w","Position",[100 100 1100 500]);
layout=tiledlayout(fig,1,2,"TileSpacing","compact","Padding","compact");
ax1=nexttile(layout);hold(ax1,"on");
frames=S(1).relativeFrame; onsetFrames=frames>=-2 & frames<=2;
for m=1:numel(S)
 if S(m).nAssemblies==0,continue,end
 label=sprintf("%s (%d fish; %d assemblies)",displayNames(m),S(m).nFish,S(m).nAssemblies);
 shadedTrace(ax1,frames(onsetFrames),S(m).fractionMean(onsetFrames),S(m).fractionSEM(onsetFrames),colors(m,:),label);
end
formatAxis(ax1,"Fraction of members active");xlim(ax1,[-2 2]);xticks(ax1,-2:2);
ylim(ax1,[0 signalUpperLimit(S)]);title(ax1,"Assembly-onset response");
legend(ax1,"Location","best","Interpreter","none");
ax2=nexttile(layout);hold(ax2,"on");lateStart=max(15,frames(1));lateFrames=frames(frames>=lateStart);
if isempty(lateFrames),lateFrames=frames;end
upper=0;
for m=1:numel(S)
 ref=outside.summaryTable(outside.summaryTable.morph==S(m).morph & outside.summaryTable.metric=="fractionActive",:);
 if isempty(ref),continue,end
 x=lateFrames(:);mu=repmat(ref.mean,size(x));sem=repmat(ref.SEM,size(x));
 shadedTrace(ax2,x,mu,sem,colors(m,:),sprintf("%s all outside frames (%d fish)",displayNames(m),ref.nFish));
 upper=max(upper,ref.mean+ref.SEM);
end
grid(ax2,"on");box(ax2,"off");xlabel(ax2,"Display interval (not onset-aligned)");
ylabel(ax2,"Fraction of members active");xlim(ax2,[lateFrames(1) lateFrames(end)]);
xticks(ax2,displayFrameTicks(lateFrames));ylim(ax2,[0 1.12*upper]);
title(ax2,"All outside-episode frames");legend(ax2,"Location","best","Interpreter","none");
if showBinaryEventsPerFrame
 ylabel(ax1,"Binary-positive frames per neuron per frame");
 ylabel(ax2,"Binary-positive frames per neuron per frame");
 stem="assembly_triggered_averages_across_morphs_binary_events_per_neuron_per_frame";
else
 stem="assembly_triggered_averages_across_morphs";
end
title(layout,"Assembly activation and outside-episode activity | equal fish weights | mean +/- SEM");
subtitle(layout,"Separate estimands and y-scales; the right panel is not an onset-conditioned continuation");
saveFigure(fig,outDir,char(stem),resolution);close(fig);
end

function makeSignalFigure(S,outDir,visibility,resolution,isControl,outside,comparison,showBinaryEventsPerFrame)
colors=assembly_figure_style().morphColors;
displayNames=["Surface","Molino","Pachón"];
fig=figure('Visible',char(visibility),'Color','w','Position',[100 100 850 550]);
layout=tiledlayout(fig,1,1,'TileSpacing','compact','Padding','compact');
ax=nexttile(layout);hold(ax,'on');
for m=1:numel(S)
 if S(m).nAssemblies==0,continue,end
 label=sprintf('%s (%d fish; %d assemblies)',displayNames(m),S(m).nFish,S(m).nAssemblies);
 if isControl
  label=sprintf('%s windows (%d fish; %d assemblies; %d windows)',displayNames(m),S(m).nFish,S(m).nAssemblies,comparison.nWindows(m));
 end
 shadedTrace(ax,S(m).relativeFrame,S(m).fractionMean,S(m).fractionSEM,colors(m,:),label);
end
if showBinaryEventsPerFrame
 yLabel='Binary events per neuron per frame';
else
 yLabel='Fraction of assembly members active';
end
formatAxis(ax,yLabel);
frames=S(1).relativeFrame;
xticks(ax,displayFrameTicks(frames));xlim(ax,[frames(1) frames(end)]);
% Do not clip the post-episode tail. The old fixed [0.3 0.7] limits hid
% precisely the outside-activation range (roughly 0.05--0.2) and made the
% binary-event plot look as if all morphs converged after the onset peak.
upper=signalUpperLimit(S);
ylim(ax,[0 upper]);
if isControl
 upper=0;
 for m=1:numel(S)
  values=S(m).fractionMean+S(m).fractionSEM;
  upper=max([upper;values(isfinite(values))]);
  ref=outside.summaryTable(outside.summaryTable.morph==S(m).morph & outside.summaryTable.metric=="fractionActive",:);
  if ~isempty(ref)&&isfinite(ref.mean)
   yline(ax,ref.mean,'--','Color',colors(m,:),'LineWidth',1.8, ...
    'DisplayName',sprintf('%s all outside frames (%d fish; %d assemblies)',displayNames(m),ref.nFish,comparison.nOutsideAssemblies(m)));
   upper=max(upper,ref.mean);
  end
 end
 if upper==0,upper=0.01;end % Nondegenerate axis for an entirely zero baseline.
 ylim(ax,[0 1.08*upper]);
 xlabel(ax,'Frames from random pseudo-onset');
 lines=findall(ax,'Type','constantline','Label','Onset');set(lines,'Label','Pseudo-onset');
end
legend(ax,'Location','best','Interpreter','none');
if isControl
 if showBinaryEventsPerFrame
  title(layout,'Binary activity in random non-activation windows | equal fish weights | mean +/- SEM');
 else
  title(layout,'Random non-activation windows | equal fish weights | mean +/- SEM');
 end
 notes={'Solid: entire window inactive; dashed: all outside-frame mean', ...
  'Full-window selection excludes short gaps and changes frame weighting'};
 stem='assembly_triggered_averages_outside_activation_across_morphs';
else
 if showBinaryEventsPerFrame
  title(layout,'Assembly onset-aligned binary activity | equal fish weights | mean +/- SEM');
 else
  title(layout,'Assembly onset-aligned activity | equal fish weights | mean +/- SEM');
 end
 notes={};
 stem='assembly_triggered_averages_across_morphs';
end
if showBinaryEventsPerFrame
 stem=[stem '_binary_events_per_neuron_per_frame'];
 notes=[notes,{'Each binary-positive member frame counts as one event; per-frame rate equals active fraction.'}];
end
if ~isempty(notes),subtitle(layout,notes);end
saveFigure(fig,outDir,stem,resolution);close(fig);
end

function ticks=displayFrameTicks(frames)
% Keep long (-2:50 by default) onset axes readable.
frames=frames(:)';
if numel(frames)<=15,ticks=frames;return,end
step=max(5,5*ceil((frames(end)-frames(1))/50));
ticks=unique([frames(1),0,step:step:frames(end),frames(end)]);
ticks=ticks(ticks>=frames(1)&ticks<=frames(end));
end

function upper=signalUpperLimit(S)
upper=0;
for m=1:numel(S)
 values=S(m).fractionMean+S(m).fractionSEM;
 values=values(isfinite(values));
 if ~isempty(values),upper=max(upper,max(values));end
end
if upper<=0,upper=0.01;else,upper=min(1,1.08*upper);end
end

function shadedTrace(ax,x,mu,sem,color,label)
x=x(:);mu=mu(:);sem=sem(:);
% Preserve conditional-category gaps in both curves and SEM patches.
edges=diff([false;isfinite(mu)&isfinite(sem);false]);
starts=find(edges==1);stops=find(edges==-1)-1;
for k=1:numel(starts)
 ix=starts(k):stops(k);xx=x(ix);lo=mu(ix)-sem(ix);hi=mu(ix)+sem(ix);
 fill(ax,[xx;flipud(xx)],[lo;flipud(hi)],color,'FaceAlpha',0.13, ...
  'EdgeColor','none','HandleVisibility','off');
end
plot(ax,x,mu,'Color',color,'LineWidth',2.5,'DisplayName',label);
end

function lengths=binaryRunLengths(B)
lengths=zeros(0,1);
for neuron=1:size(B,2)
 d=diff([false;logical(B(:,neuron));false]);
 starts=find(d==1);stops=find(d==-1)-1;
 lengths=[lengths;stops-starts+1]; %#ok<AGROW>
end
end

function [summary,distribution]=makeRunLengthFigure(values,morphs,outDir,visibility,resolution)
colors=assembly_figure_style().morphColors;
displayNames=["Surface","Molino","Pachón"];maxShown=30;
summaryRows=cell(numel(morphs),1);distributionRows=cell(numel(morphs),1);
fig=figure('Visible',char(visibility),'Color','w','Position',[100 100 760 520]);ax=axes(fig);hold(ax,'on');
for m=1:numel(morphs)
 L=double(values{m}(:));n=numel(L);
 if n==0
  stats=nan(1,7);prob=nan(maxShown,1);
 else
  stats=[mean(L),median(L),prctile(L,75),prctile(L,90),prctile(L,95),100*mean(L==10),100*mean(L>=9&L<=11)];
  prob=accumarray(min(L,maxShown+1),1,[maxShown+1 1])/n;prob=prob(1:maxShown);
  plot(ax,1:maxShown,100*prob,'Color',colors(m,:),'LineWidth',2.5,'DisplayName',sprintf('%s (n=%d runs)',displayNames(m),n));
 end
 summaryRows{m}=table(morphs(m),n,stats(1),stats(2),stats(3),stats(4),stats(5),stats(6),stats(7), ...
  'VariableNames',{'morph','nRuns','meanFrames','medianFrames','p75Frames','p90Frames','p95Frames','percentExactly10Frames','percent9To11Frames'});
 distributionRows{m}=table(repmat(morphs(m),maxShown,1),(1:maxShown)',prob, ...
  'VariableNames',{'morph','runLengthFrames','probability'});
end
xline(ax,10,'--k','10 frames','LabelVerticalAlignment','bottom','HandleVisibility','off');
xlabel(ax,'Binary-raster active-run length (frames)');ylabel(ax,'Percentage of all runs');
title(ax,'S.activity > 0 run-length distribution by morph');legend(ax,'Location','best');grid(ax,'on');box(ax,'off');xlim(ax,[0.5 maxShown+0.5]);
saveFigure(fig,outDir,'binary_raster_run_lengths_by_morph',resolution);close(fig);
summary=vertcat(summaryRows{:});distribution=vertcat(distributionRows{:});
end

function formatAxis(ax,yLabel)
xline(ax,0,'--k','Onset','LabelVerticalAlignment','bottom','HandleVisibility','off');
xlabel(ax,'Frames from assembly activation onset');ylabel(ax,yLabel);
grid(ax,'on');box(ax,'off');
end

function r=windowRecord(deltaFoFTraces,neuronActive,memberIdx,onset,relativeFrames,morph,fish,a,file,fps)
fractionSignal=mean(neuronActive(:,memberIdx),2);
   deltaFoFWindows=nan(numel(relativeFrames),numel(onset));
   activeMemberWindows=nan(numel(relativeFrames),numel(onset));
   inactiveMemberWindows=nan(numel(relativeFrames),numel(onset));
   fractionWindows=nan(numel(relativeFrames),numel(onset));
   correctedFractionWindows=nan(numel(relativeFrames),numel(onset));
   preOnsetMask=relativeFrames<0;
   assert(any(preOnsetMask),'Baseline correction requires at least one pre-onset frame.');
   episodeBaselineDeltaFoF=nan(1,numel(onset));
   for e=1:numel(onset)
    frameIdx=onset(e)+relativeFrames;
    memberDeltaFoF=deltaFoFTraces(frameIdx,memberIdx);
    memberBaseline=mean(memberDeltaFoF(preOnsetMask,:),1);
    correctedMemberDeltaFoF=memberDeltaFoF-memberBaseline;
    deltaFoFWindows(:,e)=mean(correctedMemberDeltaFoF,2);
    episodeBaselineDeltaFoF(e)=mean(memberBaseline);
    memberActive=neuronActive(frameIdx,memberIdx);
    activeValues=correctedMemberDeltaFoF;activeValues(~memberActive)=NaN;
    inactiveValues=correctedMemberDeltaFoF;inactiveValues(memberActive)=NaN;
    activeMemberWindows(:,e)=mean(activeValues,2,'omitmissing');
    inactiveMemberWindows(:,e)=mean(inactiveValues,2,'omitmissing');
    fractionWindows(:,e)=fractionSignal(frameIdx);
    correctedFractionWindows(:,e)=fractionWindows(:,e)-mean(fractionWindows(preOnsetMask,e));
   end
   r=struct('morph',morph,'fish',fish,'assemblyID',a,'file',string(file), ...
    'fps',fps,'nValidEpisodes',numel(onset),'relativeFrames',relativeFrames, ...
    'assemblyDeltaFoF',mean(deltaFoFWindows,2), ...
    'activeMemberDeltaFoF',mean(activeMemberWindows,2,'omitmissing'), ...
    'inactiveMemberDeltaFoF',mean(inactiveMemberWindows,2,'omitmissing'), ...
    'episodeBaselineDeltaFoF',episodeBaselineDeltaFoF, ...
    'fractionMembersActive',mean(fractionWindows,2), ...
    'baselineCorrectedFractionMembersActive',mean(correctedFractionWindows,2));
end

function [onsets,nEligible]=sampleControls(active,relativeFrames,requested,stream)
% Every offset must be in bounds and inactive for this assembly only.
candidates=(1-relativeFrames(1)):(numel(active)-relativeFrames(end));
indices=relativeFrames(:)+candidates;
eligible=candidates(~any(reshape(active(indices),size(indices)),1));
nEligible=numel(eligible);
onsets=eligible(randperm(stream,nEligible,min(requested,nEligible)));
assert(all(~active(relativeFrames(:)+onsets),'all'),'Control window overlaps activation.');
end

function [L,T,S]=recordTables(records,frames)
fields={'assemblyDeltaFoF','activeMemberDeltaFoF','inactiveMemberDeltaFoF', ...
 'fractionMembersActive','baselineCorrectedFractionMembersActive'};
X=cell(1,numel(fields));
for k=1:numel(fields)
 X{k}=nan(numel(records),numel(frames));
 for a=1:numel(records),X{k}(a,:)=records(a).(fields{k})(:).';end
end
[L,T,S]=makeOutputTables(records,frames,X{:});
end

function saveFigure(fig,outDir,stem,resolution)
style_assembly_figure(fig);exportgraphics(fig,fullfile(outDir,[stem '.png']),'Resolution',resolution);
% Installed MATLAB R2026a supports SVG with editable vector SEM patches.
style_assembly_figure(fig);exportgraphics(fig,fullfile(outDir,[stem '.svg']),'ContentType','vector');
end

function T=outsideEstimatorComparison(report)
% Decompose kernel minus bar using the all-frame estimator on kernel contributors.
% The window component includes finite random sampling; offsets are averaged here.
A=report.outsideActivity.assemblyTable;C=report.controlWindowSelectionTable;
rows=cell(3,1);order=["surface","molino","pachon"];
for m=1:3
 B=A(A.morph==order(m),:); W=C(C.morph==order(m)&C.obtainedWindowCount>0,:);
 included=false(height(B),1);
 for a=1:height(B)
  included(a)=any(W.fish==B.fish(a)&W.assemblyID==B.assemblyID(a));
 end
 allMean=outsideFishMean(B);matchedMean=outsideFishMean(B(included,:));
 kernelMean=mean(report.controlMorphSummaries(m).fractionMean);
 populationDelta=matchedMean-allMean;windowDelta=kernelMean-matchedMean;
 if isempty(W)
  attribution="no fully inactive windows; kernel unavailable";
 elseif all(included)
  attribution="window selection/frame weighting (including random sampling); identical population";
 else
  attribution="population inclusion and window selection/frame weighting (including random sampling)";
 end
 rows{m}=table(order(m),numel(unique(B.fish)),height(B), ...
  numel(unique(W.fish)),height(W),sum(W.obtainedWindowCount), ...
  nnz(~included),allMean,matchedMean,kernelMean,populationDelta,windowDelta,attribution, ...
  'VariableNames',{'morph','nOutsideFish','nOutsideAssemblies','nWindowFish', ...
  'nWindowAssemblies','nWindows','nOutsideAssembliesWithoutWindows', ...
  'allOutsideFrameMean','allOutsideFrameMeanWindowPopulation','kernelMeanAcrossOffsets', ...
  'populationDelta','windowSelectionDelta','attribution'});
end
T=vertcat(rows{:});
end

function mu=outsideFishMean(T)
if isempty(T),mu=NaN;return,end
G=findgroups(T.fish);
mu=mean(splitapply(@mean,T.fractionActive,G));
end
