function report = compare_assembly_properties_across_morphs(batchResults,varargin)
if nargin<1||isempty(batchResults)
 scriptDir=fileparts(mfilename('fullpath'));
 savedBatch=fullfile(scriptDir,'data_processed','assemblies','session_manifest.mat');
 assert(isfile(savedBatch),'Saved assembly manifest not found: %s',savedBatch);
 saved=load(savedBatch,'manifest','opts');
 assert(isfield(saved,'manifest')&&isfield(saved,'opts'),'Saved manifest must contain manifest and opts.');
 batchResults=struct('manifest',saved.manifest,'options',saved.opts);
end
%COMPARE_ASSEMBLY_PROPERTIES_ACROSS_MORPHS Compare p016 batch outputs by morph.
%
% report = compare_assembly_properties_across_morphs(results)
% report = compare_assembly_properties_across_morphs(results,'nCompactnessPermutations',500)
% report = compare_assembly_properties_across_morphs(results,'dunnCorrection','none')
%
% INPUT
%   batchResults  Output returned by Assem_detection. Only files whose
%                 manifest status is "processed" or "existing" are read.
%
% The fish is the unit of inference. Assembly-level tables are
% exported for inspection, but morph tests use one summary value per fish.
% Spatial compactness is a size-matched permutation z score; negative values
% mean that member ROIs are more compact than random candidate-cell sets.

p=inputParser;
addRequired(p,'batchResults',@(x)isstruct(x)&&isfield(x,'manifest')&&isfield(x,'options'));
addParameter(p,'outputDir','',@(x)ischar(x)||isstring(x));
addParameter(p,'classifiedCandidateRoot','',@(x)ischar(x)||isstring(x));
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
addParameter(p,'nCompactnessPermutations',500,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'ccgMaxLagFrames',5,@(x)isnumeric(x)&&isscalar(x)&&x>=2&&x==round(x));
addParameter(p,'transitionLagFrames',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'randomSeed',1,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'figureVisible','off',@(x)ischar(x)||isstring(x));
addParameter(p,'dunnCorrection','sidak',@(x)(ischar(x)||isstring(x))&& ...
 any(strcmpi(string(x),["sidak","none"])));
parse(p,batchResults,varargin{:}); q=p.Results;
q.dunnCorrection=lower(string(q.dunnCorrection));
if strlength(string(q.classifiedCandidateRoot))==0
 q.classifiedCandidateRoot=fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
  'data_processed','classified_candidates');
end

opts=batchResults.options;
if strlength(string(q.outputDir))==0
 outputName=['morph_comparison_' char(opts.analysisLabel)];
 if q.dunnCorrection=="none",outputName=[outputName '_dunn_uncorrected'];end
 activationTag=sprintf('_K%d_Ron%g_Roff%g_MI%g_a%g_n%d_dur%g_safe%g_seed%d', ...
  q.minActiveMembers,q.onRecruitmentFraction,q.offRecruitmentFraction, ...
  q.minMatchingIndex,q.activationNullAlpha,q.nNullShuffles, ...
  q.minDurationSeconds,q.minNullShiftSeconds,q.activationRandomSeed);
 outputName=[outputName strrep(activationTag,'.','p')];
 outputDir=fullfile(opts.figureRoot,outputName);
else
 outputDir=char(q.outputDir);
end
if ~isfolder(outputDir),mkdir(outputDir);end
rng(q.randomSeed,'twister');

M=batchResults.manifest;
use=M.status=="processed"|M.status=="existing";
M=M(use,:);
if isempty(M),error('assembly:NoBatchOutputs','No processed/existing recordings in the batch manifest.');end

fishRows=struct([]); assemblyRows=struct([]);
episodeRows=struct([]);
transitionResults=struct('morph',{},'fish',{},'counts',{},'probability',{});
failRows=struct('morph',{},'fish',{},'file',{},'message',{});
zeroAssemblyRows=struct('morph',{},'fish',{},'file',{},'reason',{});
behaviorClassPoolTables={};
activationAudits=struct('morph',{},'fish',{},'file',{},'audit',{});

fprintf('\n=== Assembly morph comparison: %d files ===\n',height(M));
for i=1:height(M)
 file=char(M.outputPath(i)); morph=string(M.morph(i)); fish=string(M.recording(i));
 fprintf('[%d/%d] %s/%s\n',i,height(M),morph,fish);
 try
  if ~isfile(file),error('Output file missing: %s',file);end
  X=load(file,'model','activation','S','provenance');
  assert(isfield(X,'model')&&isfield(X,'S'),'Required model/S variables missing.');
  member=logical(X.model.membership);
  assert(isfield(X.model,'candidateRoiIdx')&&isfield(X.S,'candidateRoiIdx'), ...
   'Saved model/session lacks candidate-neuron identity.');
  modelCandidateIDs=double(X.model.candidateRoiIdx(:));
  sessionCandidateIDs=double(X.S.candidateRoiIdx(:));
  assert(isequal(modelCandidateIDs,sessionCandidateIDs)&&size(member,1)==numel(modelCandidateIDs), ...
   'Assembly membership rows do not match the detected candidate-neuron population.');
  if size(member,2)==0
   removeStaleAssemblyPhasePlot(morph,fish,outputDir);
   zeroAssemblyRows(end+1)=struct('morph',morph,'fish',fish,'file',string(file),'reason',"zero assemblies detected"); %#ok<AGROW>
   fprintf('  Excluded from comparison: zero assemblies detected.\n');
   continue
  end
  assert(isfield(X.model,'candidateCellIDs')&&isfield(X.S,'candidateCellIDs')&& ...
   isequal(double(X.model.candidateCellIDs(:)),double(X.S.candidateCellIDs(:)))&& ...
   numel(X.S.candidateCellIDs)==size(member,1),'Candidate IDs are not aligned with membership rows.');
  assert(isfield(X.S,'preferredPhaseRad')&&numel(X.S.preferredPhaseRad)==size(member,1), ...
   'Preferred phases are not aligned with membership rows.');
  suppliedPhase=loadCandidatePhaseFeature(X.S,modelCandidateIDs);
  if ~isempty(suppliedPhase)
   makeAssemblyPhasePlot(member,suppliedPhase,morph,fish,outputDir,q.figureVisible);
  end
  try
   makeAssemblyAnatomyPlot(member,modelCandidateIDs,M(i,:),X,outputDir,q.figureVisible);
  catch anatomyError
   warning('assembly:AnatomyPlotFailed','%s/%s anatomy plot failed: %s',morph,fish,anatomyError.message);
  end
  activity=double(X.S.activity);
  fps=double(X.S.fps);
  if isscalar(fps)&&isfinite(fps)&&fps>0
   fpsForPhysicalTime=fps; fpsUsedFallback=false;
  else
   fpsForPhysicalTime=3.41; fpsUsedFallback=true;
   warning('assembly:InvalidFPSFallback', ...
    '%s/%s has invalid stored fps; physical-time displays use the documented 3.41-Hz fallback.',morph,fish);
  end
  activationOptions=struct('minActiveMembers',q.minActiveMembers, ...
   'onRecruitmentFraction',q.onRecruitmentFraction, ...
   'offRecruitmentFraction',q.offRecruitmentFraction, ...
   'minMatchingIndex',q.minMatchingIndex,'activationNullAlpha',q.activationNullAlpha, ...
   'nNullShuffles',q.nNullShuffles,'minDurationSeconds',q.minDurationSeconds, ...
   'minNullShiftSeconds',q.minNullShiftSeconds,'randomSeed',q.activationRandomSeed);
  currentActivation=compute_assembly_activations(member,activity,fps,activationOptions);
  activationAudits(end+1)=struct('morph',morph,'fish',fish, ... %#ok<AGROW>
   'file',string(file),'audit',currentActivation);
  active=logical(currentActivation.active);
  nCells=size(member,1); nAssemblies=size(member,2); nFrames=size(activity,1);
  assert(size(active,1)==nAssemblies&&size(active,2)==nFrames,'Activation dimensions disagree with model/S.');
  classLabels=loadAlignedBehaviorClasses(M(i,:),X.S,q.classifiedCandidateRoot);
  behaviorResult=analyze_assembly_behavior_classes(member,active,suppliedPhase,classLabels,morph,fish);
  behaviorClassPoolTables{end+1}=behaviorResult.classPoolTable; %#ok<AGROW>
  durationMin=nFrames/fps/60;

  onsets=false(size(active));
  if nFrames>0,onsets(:,1)=active(:,1);end
  if nFrames>1,onsets(:,2:end)=active(:,2:end)&~active(:,1:end-1);end
  episodeCount=sum(onsets,2);

  centroids=loadCentroids(M(i,:),X);
  sizeA=sum(member,1)'; sizeFrac=sizeA/max(nCells,1);
  compactZ=nan(nAssemblies,1);rgRatio=nan(nAssemblies,1);
  occupancy=nan(nAssemblies,1); episodeRate=nan(nAssemblies,1);
  memberSpikeRateDuring=nan(nAssemblies,1);memberSpikeRateOutside=nan(nAssemblies,1);
  withinAssemblyPhaseR=nan(nAssemblies,1); assemblyMeanPhase=nan(nAssemblies,1);
  medianEpisodeDurationByAssembly=nan(nAssemblies,1);
  medianEpisodeDurationSecondsByAssembly=nan(nAssemblies,1);
  neuronOutsideVariance=nan(nAssemblies,1); neuronOutsideSD=nan(nAssemblies,1);
  neuronOutsideCV=nan(nAssemblies,1); neuronOutsideMedian=nan(nAssemblies,1);
  temporalSpikeRateCV=nan(nAssemblies,1);
  medianSpikesPerMemberNeuronPerEpisode=nan(nAssemblies,1);
  withinAssemblyShortLagCCG=nan(nAssemblies,1);
  for a=1:nAssemblies
   ix=find(member(:,a)); occupancy(a)=mean(active(a,:)); episodeRate(a)=episodeCount(a)/durationMin;
   assemblyActiveFrames=active(a,:)';memberEvents=activity(:,ix)>0;
   if any(assemblyActiveFrames)
    memberSpikeRateDuring(a)=sum(memberEvents(assemblyActiveFrames,:),'all')/numel(ix)/(sum(assemblyActiveFrames)/fps/60);
   end
   neuronRates=nan(numel(ix),1);
   if any(~assemblyActiveFrames)&&~isempty(ix)
    memberSpikeRateOutside(a)=sum(memberEvents(~assemblyActiveFrames,:),'all')/numel(ix)/(sum(~assemblyActiveFrames)/fps/60);
    outsideDurationMin=sum(~assemblyActiveFrames)/fps/60;
    neuronRates=sum(memberEvents(~assemblyActiveFrames,:),1)'/outsideDurationMin;
   end
   neuronOutsideMedian(a)=safeMedian(neuronRates);
   if numel(neuronRates)>=2&&all(isfinite(neuronRates))
    neuronOutsideVariance(a)=var(neuronRates,0); neuronOutsideSD(a)=std(neuronRates,0);
   end
   neuronRateMean=safeMean(neuronRates);
   if isfinite(neuronRateMean)&&neuronRateMean>0
    neuronOutsideCV(a)=neuronOutsideSD(a)/neuronRateMean;
   end
   temporalSpikeRateCV(a)=oneMinuteSpikeVariability(memberEvents,fps);
   memberPhases=suppliedPhase(ix); memberPhases=memberPhases(isfinite(memberPhases));
   if ~isempty(memberPhases)
    meanSin=mean(sin(memberPhases)); meanCos=mean(cos(memberPhases));
    withinAssemblyPhaseR(a)=min(1,hypot(meanCos,meanSin));
    assemblyMeanPhase(a)=atan2(meanSin,meanCos);
   end
   episodeBounds=bounds(active(a,:));
   episodeDurations=episodeBounds(:,2)-episodeBounds(:,1)+1;
   medianEpisodeDurationByAssembly(a)=safeMedian(episodeDurations);
   medianEpisodeDurationSecondsByAssembly(a)=safeMedian(episodeDurations/fpsForPhysicalTime);
   episodeTotals=nan(size(episodeDurations));
   for e=1:size(episodeBounds,1)
    firstFrame=episodeBounds(e,1); lastFrame=episodeBounds(e,2);
    episodeTotals(e)=sum(memberEvents(firstFrame:lastFrame,:),'all');
    er=struct('morph',morph,'fish',fish,'assemblyID',a,'episodeID',e, ...
     'totalBinarySpikesPerEpisode',episodeTotals(e), ...
     'spikesPerMemberNeuronPerEpisode',episodeTotals(e)/numel(ix), ...
     'episodeDurationFrames',episodeDurations(e), ...
     'episodeDurationSeconds',episodeDurations(e)/fpsForPhysicalTime, ...
     'assemblySize',numel(ix),'fpsUsedForPhysicalTime',fpsForPhysicalTime, ...
     'fpsFallbackUsed',fpsUsedFallback);
    episodeRows=appendStruct(episodeRows,er);
   end
   medianSpikesPerMemberNeuronPerEpisode(a)=safeMedian(episodeTotals/numel(ix));
   withinAssemblyShortLagCCG(a)=memberCCGMetrics(activity(:,ix),q.ccgMaxLagFrames);
   [rgRatio(a),compactZ(a)]=spatialCompactness(ix,centroids,q.nCompactnessPermutations);
   ar=struct('morph',morph,'fish',fish,'assemblyID',a,'assemblySize',sizeA(a), ...
    'assemblySizeFraction',sizeFrac(a),'activationCoverage',occupancy(a), ...
    'activationEpisodes',episodeCount(a),'activationEpisodesPerMinute',episodeRate(a), ...
    'assemblySpikeRatePerMinute',episodeRate(a), ...
    'memberNeuronSpikeRateDuringActivationPerMinute',memberSpikeRateDuring(a), ...
    'memberNeuronSpikeRateOutsideActivationPerMinute',memberSpikeRateOutside(a), ...
    'spatialRgRatio',rgRatio(a),'spatialCompactnessZ',compactZ(a), ...
    'withinAssemblyPhaseR',withinAssemblyPhaseR(a), ...
    'withinAssemblyCircularVariance',1-withinAssemblyPhaseR(a));
   ar.medianEpisodeDurationFrames=medianEpisodeDurationByAssembly(a);
   ar.medianEpisodeDurationSeconds=medianEpisodeDurationSecondsByAssembly(a);
   ar.fpsUsedForPhysicalTime=fpsForPhysicalTime;
   ar.fpsFallbackUsed=fpsUsedFallback;
   ar.memberNeuronOutsideRateCV=neuronOutsideCV(a);
   ar.memberNeuronOutsideRateMedianSpikesPerMinute=neuronOutsideMedian(a);
   ar.assemblySpikeRateCVPerMinute=temporalSpikeRateCV(a);
   ar.medianSpikesPerMemberNeuronPerActivationEpisode=medianSpikesPerMemberNeuronPerEpisode(a);
   ar.withinAssemblyShortLagCCG=withinAssemblyShortLagCCG(a);
   ar.assemblyMeanPhase=behaviorResult.assemblyTable.assemblyMeanPhase(a);
   assemblyRows=appendStruct(assemblyRows,ar);
  end

  [transitionCounts,transitionProbability]=transitions(active,onsets,q.transitionLagFrames);
  transitionResults(end+1)=struct('morph',morph,'fish',fish, ... %#ok<AGROW>
   'counts',transitionCounts,'probability',transitionProbability);

  meanAssemblySpikeRate=mean(episodeRate);
  assemblySpikeRateSD=safeSampleStd(episodeRate);
  if isfinite(meanAssemblySpikeRate)&&meanAssemblySpikeRate>0
   assemblySpikeRateCV=assemblySpikeRateSD/meanAssemblySpikeRate;
  else
   assemblySpikeRateCV=NaN;
  end
  fr=struct('morph',morph,'fish',fish,'nCandidateCells',nCells,'durationMin',durationMin, ...
   'nAssemblies',nAssemblies, ...
   'activationEpisodesPerMinute',sum(episodeCount)/durationMin, ...
   'meanAssemblySpikeRatePerMinute',meanAssemblySpikeRate, ...
   'assemblySpikeRateSDPerMinute',assemblySpikeRateSD, ...
   'meanAssemblySpikeRatePerMinuteDuringActivationEpisodes',safeMean(memberSpikeRateDuring), ...
   'meanAssemblySpikeRatePerMinuteOutsideActivationEpisodes',safeMean(memberSpikeRateOutside), ...
   'assemblySpikeRateCV',assemblySpikeRateCV, ...
   'medianWithinAssemblyPhaseR',safeMedian(withinAssemblyPhaseR), ...
   'medianWithinAssemblyCircularVariance',safeMedian(1-withinAssemblyPhaseR), ...
   'medianEpisodeDurationFrames',safeMedian(medianEpisodeDurationByAssembly), ...
   'medianEpisodeDurationSecondsForPlot',safeMedian(medianEpisodeDurationSecondsByAssembly), ...
   'medianWithinAssemblyNeuronOutsideRateCV',safeMedian(neuronOutsideCV), ...
   'medianWithinAssemblyNeuronOutsideRateSpikesPerMinute',safeMedian(neuronOutsideMedian), ...
   'medianAssemblySpikeRateCVPerMinute',safeMedian(temporalSpikeRateCV), ...
   'medianSpikesPerMemberNeuronPerActivationEpisode',safeMedian(medianSpikesPerMemberNeuronPerEpisode), ...
   'cwCcwAntiphaseError',behaviorResult.fishTable.cwCcwAntiphaseError);
  fishRows=appendStruct(fishRows,fr);
 catch ME
  warning('assembly:ComparisonFailed','%s/%s failed: %s',morph,fish,ME.message);
  failRows(end+1)=struct('morph',morph,'fish',fish,'file',string(file),'message',string(ME.message)); %#ok<AGROW>
 end
end

if isempty(fishRows),error('assembly:NoValidOutputs','No batch output could be analysed.');end
fishTable=struct2table(fishRows); assemblyTable=toTable(assemblyRows);
episodeTable=toTable(episodeRows);
writetable(fishTable,fullfile(outputDir,'fish_level_metrics.csv'));
writetable(assemblyTable,fullfile(outputDir,'assembly_level_metrics.csv'));
writetable(episodeTable,fullfile(outputDir,'episode_level_total_binary_spikes.csv'));
newAssemblyVariables={'morph','fish','assemblyID','assemblySize', ...
 'memberNeuronOutsideRateCV', ...
 'memberNeuronOutsideRateMedianSpikesPerMinute','assemblySpikeRateCVPerMinute', ...
 'medianEpisodeDurationFrames','medianEpisodeDurationSeconds', ...
 'medianSpikesPerMemberNeuronPerActivationEpisode'};
writetable(assemblyTable(:,newAssemblyVariables),fullfile(outputDir,'assembly_level_new_metrics.csv'));
newFishVariables={'morph','fish','medianEpisodeDurationFrames','medianEpisodeDurationSecondsForPlot', ...
 'medianWithinAssemblyNeuronOutsideRateCV','medianWithinAssemblyNeuronOutsideRateSpikesPerMinute', ...
 'medianAssemblySpikeRateCVPerMinute','medianSpikesPerMemberNeuronPerActivationEpisode'};
writetable(fishTable(:,newFishVariables),fullfile(outputDir,'fish_level_new_metrics.csv'));
makeAssemblyPhaseCirclePlot(assemblyTable,outputDir,q.figureVisible);
if isempty(behaviorClassPoolTables),behaviorClassPoolTable=table;else,behaviorClassPoolTable=vertcat(behaviorClassPoolTables{:});end
writetable(behaviorClassPoolTable,fullfile(outputDir,'fish_behavior_class_pool_phases.csv'));
if isempty(failRows)
 failureTable=table(strings(0,1),strings(0,1),strings(0,1),strings(0,1), ...
  'VariableNames',{'morph','fish','file','message'});
else
 failureTable=struct2table(failRows);
end
writetable(failureTable,fullfile(outputDir,'processing_failures.csv'));
if isempty(zeroAssemblyRows)
 zeroAssemblyTable=table(strings(0,1),strings(0,1),strings(0,1),strings(0,1),'VariableNames',{'morph','fish','file','reason'});
else
 zeroAssemblyTable=struct2table(zeroAssemblyRows);
end
writetable(zeroAssemblyTable,fullfile(outputDir,'zero_assembly_recordings_excluded.csv'));
save(fullfile(outputDir,'transition_matrices.mat'),'transitionResults','-v7.3');
save(fullfile(outputDir,'assembly_activation_audit.mat'),'activationAudits','activationOptions','-v7.3');

metricNames={ ...
 'nAssemblies','activationEpisodesPerMinute', ...
 'meanAssemblySpikeRatePerMinute', ...
 'assemblySpikeRateSDPerMinute', ...
 'meanAssemblySpikeRatePerMinuteDuringActivationEpisodes', ...
 'meanAssemblySpikeRatePerMinuteOutsideActivationEpisodes', ...
 'assemblySpikeRateCV', ...
 'medianWithinAssemblyNeuronOutsideRateCV', ...
 'medianWithinAssemblyNeuronOutsideRateSpikesPerMinute', ...
 'medianAssemblySpikeRateCVPerMinute', ...
 'medianSpikesPerMemberNeuronPerActivationEpisode', ...
 'medianWithinAssemblyPhaseR', ...
 'medianWithinAssemblyCircularVariance', ...
 'medianEpisodeDurationFrames', ...
 'cwCcwAntiphaseError'};
[statistics,omnibusTests,dunnTests]=compareMetrics(fishTable,metricNames,q.dunnCorrection);
writetable(statistics,fullfile(outputDir,'morph_descriptive_summary.csv'));
writetable(omnibusTests,fullfile(outputDir,'morph_omnibus_kw_tests.csv'));
if q.dunnCorrection=="sidak"
 dunnFile='surface_referenced_dunn_sidak_tests.csv';
else
 dunnFile='surface_referenced_dunn_uncorrected_tests.csv';
end
writetable(dunnTests,fullfile(outputDir,dunnFile));
makePlots(fishTable,metricNames,omnibusTests,outputDir,q.figureVisible);
makeBarPlots(fishTable,metricNames,omnibusTests,outputDir,q.figureVisible);
makeAssemblySizePlot(assemblyTable,outputDir,q.figureVisible);

definitions={
 'nAssemblies','Number of detected assemblies.';
 'activationEpisodesPerMinute','Total assembly activation-episode onsets across all assemblies, divided by recording minutes.';
 'meanAssemblySpikeRatePerMinute','Within-fish mean assembly spike rate; an assembly spike is an activation-episode onset under the shared empirical-null recruitment/matching-index hysteresis definition.';
 'meanAssemblySpikeRatePerMinuteDuringActivationEpisodes','Within-fish mean across assemblies of the member-neuron binary event rate during frames when the assembly is active, normalized per neuron per minute.';
 'meanAssemblySpikeRatePerMinuteOutsideActivationEpisodes','Within-fish mean across assemblies of the member-neuron binary event rate during frames when the assembly is inactive, normalized per neuron per minute.';
 'assemblySpikeRateSDPerMinute','Within-fish sample standard deviation of assembly spike rates in activation-episode onsets per minute; undefined for fish with fewer than two assemblies.';
 'assemblySpikeRateCV','Within-fish coefficient of variation of assembly spike rates (sample SD divided by mean); undefined when fewer than two assemblies or mean rate is zero.';
 'medianWithinAssemblyNeuronOutsideRateCV','Within-fish median member-neuron outside-rate SD divided by mean; undefined for non-positive means.';
 'medianWithinAssemblyNeuronOutsideRateSpikesPerMinute','Within-fish median of assembly-level median member-neuron outside-episode rates.';
 'medianAssemblySpikeRateCVPerMinute','Within-fish median temporal CV of total member binary-spike rate across complete non-overlapping one-minute bins.';
 'medianSpikesPerMemberNeuronPerActivationEpisode','Within-fish median of assembly-level median episode spike totals divided by assembly size.';
 'medianWithinAssemblyPhaseR','Within-fish median mean resultant length of member-neuron preferred phases: R = abs(mean(exp(1i*phase))). One means identical phases; zero means phase vectors cancel (not necessarily a uniform distribution). Only finite phases contribute; assemblies without finite phases are excluded.';
 'medianWithinAssemblyCircularVariance','Within-fish median circular variance of member-neuron preferred phases: V = 1 - abs(mean(exp(1i*phase))). Zero means identical phases; larger values mean less alignment, reaching one when phase vectors cancel (not necessarily a uniform distribution). Only finite phases contribute; assemblies without finite phases are excluded.';
 'medianEpisodeDurationFrames','Within-fish median of the per-assembly median contiguous activation-episode duration in frames.';
 'cwCcwAntiphaseError','Absolute deviation in radians of the CW-versus-CCW assembly-pool phase separation from pi; zero is perfect antiphase.'};
definitionTable=cell2table(definitions,'VariableNames',{'metric','definition'});
writetable(definitionTable,fullfile(outputDir,'metric_definitions.csv'));

report=struct('fishLevel',fishTable,'assemblyLevel',assemblyTable, ...
 'activationEpisodes',episodeTable,'behaviorClassPoolPhases',behaviorClassPoolTable, ...
 'morphDescriptiveSummary',statistics,'morphOmnibusKWTests',omnibusTests, ...
 'surfaceReferencedDunnTests',dunnTests,'failures',failureTable, ...
 'zeroAssemblyRecordingsExcluded',zeroAssemblyTable, ...
 'transitionResults',transitionResults,'activationAudits',{activationAudits}, ...
 'outputDir',string(outputDir),'settings',q);
report.settings.activationOptions=activationOptions;
report.settings.activationDefinition="empirical circular-shift coactivation null with recruitment/matching-index onset criteria and hysteretic offsets";
report.settings.neuronActiveDefinition="S.activity > 0 without temporal extension";
report.settings.physicalTimeFPSFallback=3.41;
report.settings.oneMinuteBinDefinition="Complete non-overlapping round(60*fps)-frame bins; trailing partial bin excluded";
report.settings.durationInference="All duration morph p-values use medianEpisodeDurationFrames; seconds are display-only";
if q.dunnCorrection=="sidak"
 report.surfaceReferencedDunnSidakTests=dunnTests;
end
save(fullfile(outputDir,'assembly_morph_comparison_report.mat'),'report','-v7.3');
fprintf('Comparison complete. Outputs: %s\n',outputDir);
end

function y=safeMean(x)
x=double(x(:));x=x(isfinite(x));
if isempty(x),y=NaN;else,y=mean(x);end
end

function y=safeSampleStd(x)
x=double(x(:));x=x(isfinite(x));
if numel(x)<2,y=NaN;else,y=std(x,0);end
end

function A=appendStruct(A,x)
if isempty(A),A=x;else,A(end+1)=x;end %#ok<AGROW>
end

function labels=loadAlignedBehaviorClasses(manifestRow,S,classifiedRoot)
[~,candidateStem,candidateExt]=fileparts(char(manifestRow.candidatePath));
classifiedFile=fullfile(char(classifiedRoot),char(lower(string(manifestRow.morph))), ...
 [candidateStem '_classified' candidateExt]);
assert(isfile(classifiedFile),'Classified candidate file missing: %s',classifiedFile);
C=load(classifiedFile,'candidate_cell_ids','neuron_class_label');
assert(isfield(C,'candidate_cell_ids')&&isfield(C,'neuron_class_label'), ...
 'Classified candidate file lacks candidate_cell_ids or neuron_class_label: %s',classifiedFile);
classifiedIDs=double(C.candidate_cell_ids(:));sessionIDs=double(S.candidateCellIDs(:));
assert(numel(unique(classifiedIDs))==numel(classifiedIDs), ...
 'Duplicate candidate IDs in classified candidate file: %s',classifiedFile);
[found,location]=ismember(sessionIDs,classifiedIDs);
assert(all(found),'%d assembly candidate ID(s) are absent from %s',sum(~found),classifiedFile);
allLabels=string(C.neuron_class_label(:));
assert(numel(allLabels)==numel(classifiedIDs),'Class labels and candidate IDs differ in length: %s',classifiedFile);
labels=allLabels(location);
end

function T=toTable(S)
if isempty(S),T=table;else,T=struct2table(S);end
end

function x=safeMedian(x)
x=x(isfinite(x));if isempty(x),x=NaN;else,x=median(x);end
end

function phaseValues=loadCandidatePhaseFeature(S,modelRoiIDs)
phaseValues=[];
if ~isfield(S,'preferredPhaseRad'),return,end
assert(isfield(S,'candidateRoiIdx')&&isequal(double(S.candidateRoiIdx(:)),modelRoiIDs(:)), ...
 'Preferred phases do not share the assembly membership row order.');
phaseValues=double(S.preferredPhaseRad(:));
assert(numel(phaseValues)==numel(modelRoiIDs)&&all(phaseValues(isfinite(phaseValues))>-pi&phaseValues(isfinite(phaseValues))<=pi), ...
 'Harmonic preferred phases are missing, misaligned, or outside (-pi,pi].');
end

function makeAssemblyPhasePlot(member,phaseValues,morph,fish,outDir,visibility)
nAssemblies=size(member,2);
figHeight=max(420,min(1200,90+55*max(nAssemblies,1)));
fig=figure('Visible',char(visibility),'Color','w','Position',[100 100 760 figHeight]);
ax=axes(fig);hold(ax,'on');
finitePopulation=phaseValues(isfinite(phaseValues));
if isempty(finitePopulation)
 text(ax,0.5,0.5,'No finite harmonic preferred phases', ...
  'Units','normalized','HorizontalAlignment','center');
else
 for a=1:nAssemblies
  v=phaseValues(member(:,a)&isfinite(phaseValues));
  if isempty(v),continue,end
  y=a+linspace(-0.12,0.12,numel(v))';
  scatter(ax,v,y,28,[0.25 0.45 0.75],'filled','MarkerFaceAlpha',0.75);
  plot(ax,circularMean(v),a,'kd','MarkerFaceColor',[0.95 0.65 0.15],'MarkerSize',7);
 end
 xline(ax,circularMean(finitePopulation),'--','Candidate-population circular mean', ...
  'Color',[0.35 0.35 0.35],'LabelVerticalAlignment','bottom');
end
if nAssemblies>0
 ylim(ax,[0.5 nAssemblies+0.5]);yticks(ax,1:nAssemblies);ylabel(ax,'Assembly');
else
 ylim(ax,[0 1]);yticks(ax,[]);
end
xlim(ax,[-pi pi]);xticks(ax,[-pi -pi/2 0 pi/2 pi]);xticklabels(ax,{'-pi','-pi/2','0','pi/2','pi'});xlabel(ax,'First-harmonic preferred phase (rad)');
title(ax,sprintf('%s/%s | member-neuron harmonic preferred phase by assembly',morph,fish),'Interpreter','none');
grid(ax,'on');box(ax,'off');
fileStem=char(morph+"__"+fish+"__assembly_harmonic_preferred_phase");
style_assembly_figure(fig);exportgraphics(fig,fullfile(outDir,string(fileStem)+".png"),'Resolution',300);
style_assembly_figure(fig);exportgraphics(fig,fullfile(outDir,string(fileStem)+".svg"),'ContentType','vector');close(fig);
end

function x=circularMean(v)
x=atan2(mean(sin(v)),mean(cos(v)));if x<=-pi,x=pi;end
end

function makeAssemblyPhaseCirclePlot(assemblyTable,outDir,visibility)
% One unit circle per morph; every point is one assembly from one fish.
morphs=["surface","molino","pachon"];
fig=figure('Visible',char(visibility),'Color','w','Position',[100 100 1200 430]);
layout=tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');
phaseMap=hsv(256);theta=linspace(-pi,pi,721);
for m=1:numel(morphs)
 ax=nexttile(layout,m);hold(ax,'on');
 plot(ax,cos(theta),sin(theta),'-','Color',[0.25 0.25 0.25],'LineWidth',1.4);
 plot(ax,[-1.08 1.08],[0 0],'-','Color',[0.82 0.82 0.82],'LineWidth',0.75);
 plot(ax,[0 0],[-1.08 1.08],'-','Color',[0.82 0.82 0.82],'LineWidth',0.75);
 use=string(assemblyTable.morph)==morphs(m)&isfinite(assemblyTable.assemblyMeanPhase);
 phase=double(assemblyTable.assemblyMeanPhase(use));
 scatter(ax,cos(phase),sin(phase),48,phase,'filled', ...
  'MarkerEdgeColor',[0.15 0.15 0.15],'LineWidth',0.45,'MarkerFaceAlpha',0.82);
 text(ax,1.13,0,'0','HorizontalAlignment','left','VerticalAlignment','middle');
 text(ax,-1.13,0,'\pm\pi','HorizontalAlignment','right','VerticalAlignment','middle','Interpreter','tex');
 text(ax,0,1.13,'\pi/2','HorizontalAlignment','center','VerticalAlignment','bottom','Interpreter','tex');
 text(ax,0,-1.13,'-\pi/2','HorizontalAlignment','center','VerticalAlignment','top','Interpreter','tex');
 nFish=numel(unique(string(assemblyTable.fish(use))));
 title(ax,sprintf('%s  |  %d assemblies, %d fish',morphs(m),numel(phase),nFish),'Interpreter','none');
 axis(ax,'equal');xlim(ax,[-1.28 1.28]);ylim(ax,[-1.28 1.28]);
 axis(ax,'off');colormap(ax,phaseMap);clim(ax,[-pi pi]);
end
cb=colorbar(ax,'eastoutside');cb.Ticks=[-pi -pi/2 0 pi/2 pi];
cb.TickLabels={'-pi','-pi/2','0','pi/2','pi'};
cb.Label.String='Assembly mean preferred phase (rad)';
title(layout,'Assembly preferred phases across morphs');
style_assembly_figure(fig);exportgraphics(fig,fullfile(outDir,'assembly_mean_phase_circles_by_morph.png'),'Resolution',300);
style_assembly_figure(fig);exportgraphics(fig,fullfile(outDir,'assembly_mean_phase_circles_by_morph.svg'),'ContentType','vector');close(fig);
end

function removeStaleAssemblyPhasePlot(morph,fish,outDir)
stem=char(morph+"__"+fish+"__assembly_harmonic_preferred_phase");
for extension={'.png','.svg'}
 file=fullfile(outDir,[stem extension{1}]);
 if isfile(file),delete(file);end
end
end

function centroids=loadCentroids(row,X)
centroids=[]; path='';
if ismember('hdMetricsPath',row.Properties.VariableNames),path=char(row.hdMetricsPath);end
if ~isfile(path)&&isfield(X,'provenance')&&isfield(X.provenance,'options')
 root=X.provenance.options.dataRoot; d=dir(fullfile(root,char(row.morph),char(row.recording),'*hdMetrics_*.mat'));
 d=d(~[d.isdir]);if numel(d)==1,path=fullfile(d.folder,d.name);end
end
if ~isfile(path),return,end
available=who('-file',path);
if ~ismember('hdMetrics',available),return,end
H=load(path,'hdMetrics');
if ~isfield(H,'hdMetrics')||~isfield(H.hdMetrics,'centroids'),return,end
centroids=double(H.hdMetrics.centroids);
if size(centroids,2)~=2&&size(centroids,1)==2,centroids=centroids.';end
if size(centroids,2)<2,centroids=[];return,end
centroids=centroids(:,1:2);
roi=double(X.S.candidateRoiIdx(:)); valid=roi>=1&roi<=size(centroids,1);
mapped=nan(numel(roi),2);mapped(valid,:)=centroids(roi(valid),:);centroids=mapped;
end

function [ratio,z]=spatialCompactness(members,centroids,nPerm)
ratio=NaN;z=NaN;
if isempty(centroids)||numel(members)<2,return,end
usable=find(all(isfinite(centroids),2));members=intersect(members,usable,'stable');
if numel(members)<2||numel(usable)<numel(members),return,end
popRg=radiusGyration(centroids(usable,:));obs=radiusGyration(centroids(members,:));
if isfinite(popRg)&&popRg>0,ratio=obs/popRg;end
if nPerm<2||numel(usable)==numel(members),return,end
null=nan(nPerm,1);
for k=1:nPerm,pick=usable(randperm(numel(usable),numel(members)));null(k)=radiusGyration(centroids(pick,:));end
s=std(null,'omitnan');if isfinite(s)&&s>0,z=(obs-mean(null,'omitnan'))/s;end
end

function r=radiusGyration(x)
x=x(all(isfinite(x),2),:);if size(x,1)<2,r=NaN;else,d=x-mean(x,1);r=sqrt(mean(sum(d.^2,2)));end
end

function shortLagPeak=memberCCGMetrics(traces,maxLag)
shortLagPeak=NaN;
if size(traces,2)<2||size(traces,1)<=maxLag,return,end
peaks=[];
for i=1:size(traces,2)-1
 for j=i+1:size(traces,2)
  lagCorrelation=nan(2*maxLag+1,1);lags=(-maxLag:maxLag)';
  for k=1:numel(lags)
   lag=lags(k);
   if lag>=0,x=traces(1:end-lag,i);y=traces(1+lag:end,j);
   else,x=traces(1-lag:end,i);y=traces(1:end+lag,j);end
   lagCorrelation(k)=pearsonFinite(x,y);
  end
  central=lagCorrelation(abs(lags)<=1);central=central(isfinite(central));
  if ~isempty(central)
   pairPeak=max(central);peaks(end+1,1)=pairPeak; %#ok<AGROW>
  end
 end
end
shortLagPeak=safeMedian(peaks);
end

function r=pearsonFinite(x,y)
ok=isfinite(x)&isfinite(y);x=x(ok);y=y(ok);
if numel(x)<2||std(x)==0||std(y)==0,r=NaN;else,C=corrcoef(x,y);r=C(1,2);end
end

function rateCV=oneMinuteSpikeVariability(memberEvents,fps)
rateCV=NaN;
if ~isscalar(fps)||~isfinite(fps)||fps<=0,return,end
binFrames=round(60*fps);nBins=floor(size(memberEvents,1)/binFrames);
if binFrames<1||nBins<2,return,end
counts=zeros(nBins,1);
for bin=1:nBins
 frames=(bin-1)*binFrames+(1:binFrames);
 counts(bin)=sum(memberEvents(frames,:),'all');
end
rates=counts/(binFrames/fps/60);meanRate=mean(rates);
if meanRate>0,rateCV=std(rates,0)/meanRate;end
end

function [counts,P]=transitions(active,onsets,lag)
n=size(active,1);counts=zeros(n);
for t=1:max(0,size(active,2)-lag)
 s=find(active(:,t));d=find(onsets(:,t+lag));
 for i=1:numel(s),for j=1:numel(d),if s(i)~=d(j),counts(s(i),d(j))=counts(s(i),d(j))+1;end,end,end
end
den=sum(counts,2);P=zeros(n);ok=den>0;P(ok,:)=counts(ok,:)./den(ok);
end

function b=bounds(x)
d=diff([false,logical(x),false]);b=[find(d==1)',find(d==-1)'-1];
end

function [summary,omnibus,dunn]=compareMetrics(T,names,dunnCorrection)
O=struct([]);K=struct([]);D=struct([]);groups=["surface","molino","pachon"];
contrasts=["surface","molino";"surface","pachon"];
for k=1:numel(names)
 name=names{k};y=double(T.(name));g=string(T.morph);
 values=cell(1,3);med=nan(1,3);nn=zeros(1,3);
 for j=1:3
  values{j}=y(g==groups(j)&isfinite(y));
  nn(j)=numel(values{j});med(j)=safeMedian(values{j});
 end
 o=struct('metric',string(name),'surfaceN',nn(1),'molinoN',nn(2),'pachonN',nn(3), ...
  'surfaceMedian',med(1),'molinoMedian',med(2),'pachonMedian',med(3));
 O=appendStruct(O,o);

 meanRanks=nan(1,3);dunnVarianceFactor=NaN;
 eligible=all(nn>=3);
 if eligible
  pooledValues=vertcat(values{:});
  pooledGroups=[repmat(groups(1),nn(1),1);repmat(groups(2),nn(2),1);repmat(groups(3),nn(3),1)];
  pooledRanks=tiedrank(pooledValues);
  for j=1:3,meanRanks(j)=mean(pooledRanks(pooledGroups==groups(j)));end
  [~,~,tieIndex]=unique(pooledValues);
  tieCounts=accumarray(tieIndex,1);
  tieTerm=sum(tieCounts.^3-tieCounts);
  nTotal=numel(pooledValues);
  tieCorrection=1-tieTerm/(nTotal^3-nTotal);
  dunnVarianceFactor=nTotal*(nTotal+1)*tieCorrection/12;
 end
 for c=1:2
  group1=contrasts(c,1);group2=contrasts(c,2);
  i=find(groups==group1,1);j=find(groups==group2,1);
  kwStatistic=NaN;kwDf=NaN;kwP=NaN;
  if nn(i)>=3&&nn(j)>=3
   pairValues=[values{i};values{j}];
   pairGroups=[repmat(group1,nn(i),1);repmat(group2,nn(j),1)];
   try
    [kwP,kwTable]=kruskalwallis(pairValues,cellstr(pairGroups),'off');
    kwDf=kwTable{2,3};kwStatistic=kwTable{2,5};
   catch
    kwP=NaN;kwDf=NaN;kwStatistic=NaN;
   end
  end
  kwRow=struct('metric',string(name),'group1',group1,'group2',group2, ...
   'n1',nn(i),'n2',nn(j),'median1',med(i),'median2',med(j), ...
   'kwStatistic',kwStatistic,'degreesOfFreedom',kwDf,'pOmnibus',kwP, ...
   'test',"two-group Kruskal-Wallis");
  K=appendStruct(K,kwRow);
 end

 raw=nan(2,1);temp=struct([]);
 for c=1:2
  group1=contrasts(c,1);group2=contrasts(c,2);
  i=find(groups==group1,1);j=find(groups==group2,1);dunnZ=NaN;pRaw=NaN;
  if eligible&&isfinite(dunnVarianceFactor)&&dunnVarianceFactor>0
   standardError=sqrt(dunnVarianceFactor*(1/nn(i)+1/nn(j)));
   dunnZ=(meanRanks(i)-meanRanks(j))/standardError;
   pRaw=erfc(abs(dunnZ)/sqrt(2));
  end
  r=struct('metric',string(name),'group1',group1,'group2',group2, ...
   'n1',nn(i),'n2',nn(j),'median1',med(i),'median2',med(j), ...
   'medianDifference',med(i)-med(j),'dunnZ',dunnZ,'pRaw',pRaw,'pAdjusted',NaN, ...
   'test',"Dunn test using ranks pooled across all three morphs", ...
   'correction',dunnCorrectionLabel(dunnCorrection));
  temp=appendStruct(temp,r);raw(c)=pRaw;
 end
 if dunnCorrection=="sidak"
  pAdjusted=1-(1-raw).^2;pAdjusted=min(pAdjusted,1);
 else
  pAdjusted=raw;
 end
 for c=1:2,temp(c).pAdjusted=pAdjusted(c);D=appendStruct(D,temp(c));end
end
summary=toTable(O);omnibus=toTable(K);dunn=toTable(D);
end

function makePlots(T,names,omnibus,outDir,visibility)
labels=string(T.morph);order=["surface","molino","pachon"];
colors=assembly_figure_style().morphColors;
for k=1:numel(names)
 name=names{k};fig=figure('Visible',char(visibility),'Color','w','Position',[100 100 620 480]);ax=axes(fig);hold(ax,'on');
 y=double(T.(name));
 pKW=[contrastKwP(omnibus,name,order(1),order(2)),contrastKwP(omnibus,name,order(1),order(3))];
 for j=1:numel(order)
  v=y(labels==order(j)&isfinite(y));
  if isempty(v),continue,end
  x=j+assembly_deterministic_jitter(numel(v),.11);scatter(ax,x,v,25,'MarkerFaceColor',assembly_mix_with_white(colors(j,:),.38),'MarkerEdgeColor','k','LineWidth',.45,'MarkerFaceAlpha',.82);
  med=median(v);plot(ax,j,med,'d','Color','k','MarkerFaceColor',colors(j,:),'MarkerSize',8,'LineWidth',1.1);
 end
 finiteY=y(isfinite(y));
 if isempty(finiteY),lo=0;hi=1;else,lo=min(finiteY);hi=max(finiteY);end
 span=hi-lo;if ~isfinite(span)||span==0,span=max(abs([lo hi]));end
 if ~isfinite(span)||span==0,span=1;end
 comparisons=[1 2;1 3];
 for c=1:size(comparisons,1)
  a=comparisons(c,1);b=comparisons(c,2);
  pOmnibus=pKW(c);
  bracketY=hi+span*(0.10+0.12*(c-1));tick=0.035*span;
  plot(ax,[a a b b],[bracketY-tick bracketY bracketY bracketY-tick],'-k','LineWidth',1);
  text(ax,(a+b)/2,bracketY+0.01*span,significanceLabel(pOmnibus), ...
   'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold');
 end
 ylim(ax,[lo-0.06*span hi+0.30*span]);
 xlim(ax,[.5 3.5]);xticks(ax,1:3);xticklabels(ax,order);ylabel(ax,name,'Interpreter','none');grid(ax,'on');box(ax,'off');
 title(ax,name,'Interpreter','none');
 text(ax,0.98,0.98,sprintf('KW tests\nSurface-Molino: p=%.3g\nSurface-Pachon: p=%.3g',pKW(1),pKW(2)), ...
  'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
  'BackgroundColor','w','EdgeColor','k','Margin',6,'Interpreter','none');
 style_assembly_figure(fig);exportgraphics(fig,fullfile(outDir,string(name)+".png"),'Resolution',300);style_assembly_figure(fig);exportgraphics(fig,fullfile(outDir,string(name)+".svg"),'ContentType','vector');close(fig);
end
end

function label=significanceLabel(p)
if ~isfinite(p),label='n/a';
elseif p<=0.001,label='***';
elseif p<=0.01,label='**';
elseif p<=0.05,label='*';
else,label='ns';
end
end


function makeBarPlots(T,names,omnibus,outDir,visibility)
labels=string(T.morph);order=["surface","molino","pachon"];
colors=assembly_figure_style().morphColors;
for k=1:numel(names)
 name=names{k};y=double(T.(name));means=nan(1,3);sems=nan(1,3);
 pKW=[contrastKwP(omnibus,name,order(1),order(2)),contrastKwP(omnibus,name,order(1),order(3))];
 fig=figure('Visible',char(visibility),'Color','w','Position',[100 100 620 480]);ax=axes(fig);hold(ax,'on');
 for j=1:3
  v=y(labels==order(j)&isfinite(y));
  if isempty(v),continue,end
  means(j)=mean(v);if numel(v)>1,sems(j)=std(v)/sqrt(numel(v));else,sems(j)=0;end
 end
 bars=bar(ax,1:3,means,.46,'FaceColor','flat','EdgeColor','k','LineWidth',.8,'FaceAlpha',.78);bars.CData=colors;
 errorbar(ax,1:3,means,sems,'k','LineStyle','none','LineWidth',1.7,'CapSize',8);
 for j=1:3
  v=y(labels==order(j)&isfinite(y));if isempty(v),continue,end
  x=j+assembly_deterministic_jitter(numel(v),.11);scatter(ax,x,v,25,'MarkerFaceColor',assembly_mix_with_white(colors(j,:),.38),'MarkerEdgeColor','k','LineWidth',.45,'MarkerFaceAlpha',.82);
 end
 finiteY=y(isfinite(y));
 if isempty(finiteY),lo=0;hi=1;else,lo=min([finiteY;means(:)-sems(:)],[],'omitnan');hi=max([finiteY;means(:)+sems(:)],[],'omitnan');end
 span=hi-lo;if ~isfinite(span)||span==0,span=max(abs([lo hi]));end
 if ~isfinite(span)||span==0,span=1;end
 comparisons=[1 2;1 3];
 for c=1:2
  a=comparisons(c,1);b=comparisons(c,2);
  pOmnibus=pKW(c);
  bracketY=hi+span*(0.10+0.12*(c-1));tick=0.035*span;
  plot(ax,[a a b b],[bracketY-tick bracketY bracketY bracketY-tick],'-k','LineWidth',1);
  text(ax,(a+b)/2,bracketY+0.01*span,significanceLabel(pOmnibus),'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold');
 end
 ylim(ax,[min(0,lo-0.06*span) hi+0.30*span]);xlim(ax,[.5 3.5]);xticks(ax,1:3);xticklabels(ax,order);
 ylabel(ax,[name ' (mean +/- SEM)'],'Interpreter','none');grid(ax,'on');box(ax,'off');
 title(ax,[name ' | mean +/- SEM'],'Interpreter','none');
 text(ax,0.98,0.98,sprintf('KW tests\nSurface-Molino: p=%.3g\nSurface-Pachon: p=%.3g',pKW(1),pKW(2)), ...
  'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top', ...
  'BackgroundColor','w','EdgeColor','k','Margin',6,'Interpreter','none');
 style_assembly_figure(fig);exportgraphics(fig,fullfile(outDir,string(name)+"_bar.png"),'Resolution',300);
 style_assembly_figure(fig);exportgraphics(fig,fullfile(outDir,string(name)+"_bar.svg"),'ContentType','vector');close(fig);
end
end
function makeAssemblySizePlot(assemblyTable,outDir,visibility)
order=["surface","molino","pachon"];displayNames=["Surface","Molino","Pachon"];
colors=assembly_figure_style().morphColors;
g=lower(string(assemblyTable.morph));y=double(assemblyTable.assemblySize);
fig=figure('Visible',char(visibility),'Color','w','Position',[100 100 700 500]);ax=axes(fig);hold(ax,'on');
for j=1:3
 pick=g==order(j)&isfinite(y);v=y(pick);
 if isempty(v),continue,end
 jitter=assembly_deterministic_jitter(numel(v),.11);scatter(ax,j+jitter,v,25,'MarkerFaceColor',assembly_mix_with_white(colors(j,:),.38),'MarkerEdgeColor','k','LineWidth',.45,'MarkerFaceAlpha',.82);
 med=median(v);plot(ax,j,med,'d','Color','k','MarkerFaceColor',colors(j,:),'MarkerSize',8,'LineWidth',1.1);
end
xlim(ax,[.5 3.5]);xticks(ax,1:3);xticklabels(ax,displayNames);ylabel(ax,'Neurons per assembly');
title(ax,'Assembly size by morph (individual assembly is the observational unit)');grid(ax,'on');box(ax,'off');
style_assembly_figure(fig);exportgraphics(fig,fullfile(outDir,'assemblySize_individual_assemblies.png'),'Resolution',300);
style_assembly_figure(fig);exportgraphics(fig,fullfile(outDir,'assemblySize_individual_assemblies.svg'),'ContentType','vector');close(fig);

% Bar-chart companion using the same individual assemblies as observations.
means=nan(1,3);sems=nan(1,3);counts=zeros(1,3);values=cell(1,3);
for j=1:3
 values{j}=y(g==order(j)&isfinite(y));counts(j)=numel(values{j});
 if counts(j)>0,means(j)=mean(values{j});end
 if counts(j)>1,sems(j)=std(values{j},0)/sqrt(counts(j));elseif counts(j)==1,sems(j)=0;end
end
contrasts=[1 2;1 3];testRows=struct([]);pValues=nan(1,2);
for c=1:2
 a=contrasts(c,1);b=contrasts(c,2);pValue=NaN;chiSquare=NaN;degreesOfFreedom=NaN;
 if ~isempty(values{a})&&~isempty(values{b})
  try
   groups=[repmat(order(a),counts(a),1);repmat(order(b),counts(b),1)];
   [pValue,kwTable]=kruskalwallis([values{a};values{b}],cellstr(groups),'off');
   degreesOfFreedom=kwTable{2,3};chiSquare=kwTable{2,5};
  catch
  end
 end
 pValues(c)=pValue;
 testRows=appendStruct(testRows,struct('metric',"assemblySize",'test',"two-group one-way Kruskal-Wallis", ...
  'group1',order(a),'group2',order(b),'n1Assemblies',counts(a),'n2Assemblies',counts(b), ...
  'median1',safeMedian(values{a}),'median2',safeMedian(values{b}), ...
  'chiSquare',chiSquare,'degreesOfFreedom',degreesOfFreedom,'pValue',pValue));
end
writetable(toTable(testRows),fullfile(outDir,'assemblySize_individual_assemblies_kw_tests.csv'));

fig=figure('Visible',char(visibility),'Color','w','Position',[100 100 700 500]);ax=axes(fig);hold(ax,'on');
bars=bar(ax,1:3,means,.46,'FaceColor','flat','EdgeColor','k','LineWidth',.8,'FaceAlpha',.78);bars.CData=colors;
errorbar(ax,1:3,means,sems,'k','LineStyle','none','LineWidth',1.7,'CapSize',8);
for j=1:3
 if counts(j)>0,scatter(ax,j+assembly_deterministic_jitter(counts(j),.11),values{j},25,'MarkerFaceColor',assembly_mix_with_white(colors(j,:),.38),'MarkerEdgeColor','k','LineWidth',.45,'MarkerFaceAlpha',.82);end
end
finiteY=y(isfinite(y));if isempty(finiteY),lo=0;hi=1;else,lo=min([0;finiteY]);hi=max([finiteY;means(:)+sems(:)],[],'omitnan');end
span=hi-lo;if ~isfinite(span)||span<=0,span=1;end
for c=1:2
 a=contrasts(c,1);b=contrasts(c,2);bracketY=hi+span*(.10+.12*(c-1));tick=.035*span;
 plot(ax,[a a b b],[bracketY-tick bracketY bracketY bracketY-tick],'-k','LineWidth',1);
 text(ax,(a+b)/2,bracketY+.01*span,significanceLabel(pValues(c)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold');
end
ylim(ax,[lo hi+.30*span]);xlim(ax,[.5 3.5]);xticks(ax,1:3);xticklabels(ax,displayNames);
ylabel(ax,'Neurons per assembly (mean +/- SEM)');title(ax,'Assembly size by morph (individual assembly is the observational unit)');grid(ax,'on');box(ax,'off');
text(ax,.98,.98,sprintf('One-way KW tests\nSurface-Molino: p=%.3g\nSurface-Pachon: p=%.3g',pValues(1),pValues(2)), ...
 'Units','normalized','HorizontalAlignment','right','VerticalAlignment','top','BackgroundColor','w','EdgeColor','k','Margin',6);
style_assembly_figure(fig);exportgraphics(fig,fullfile(outDir,'assemblySize_individual_assemblies_bar.png'),'Resolution',300);
style_assembly_figure(fig);exportgraphics(fig,fullfile(outDir,'assemblySize_individual_assemblies_bar.svg'),'ContentType','vector');close(fig);
end

function label=dunnCorrectionLabel(correction)
if correction=="sidak"
 label="Šidák-corrected across the two Surface-referenced contrasts";
else
 label="none (raw two-sided Dunn p-value)";
end
end

function p=contrastKwP(omnibus,name,group1,group2)
pick=omnibus.metric==string(name)&omnibus.group1==group1&omnibus.group2==group2;
if ~any(pick),pick=omnibus.metric==string(name)&omnibus.group1==group2&omnibus.group2==group1;end
if any(pick),p=omnibus.pOmnibus(find(pick,1));else,p=NaN;end
end

function makeAssemblyAnatomyPlot(member,roi,row,X,outDir,visibility)
% Map membership rows to original ALL_CELLS ROI masks in image pixel space.
folder='';
if ismember('recordingPath',row.Properties.VariableNames),folder=char(row.recordingPath);end
if ~isfolder(folder)&&isfield(X,'provenance')&&isfield(X.provenance,'options')&&isfield(X.provenance.options,'dataRoot')
 folder=fullfile(char(X.provenance.options.dataRoot),char(row.morph),char(row.recording));
end
if ~isfolder(folder),folder=fullfile('/home/blanche/data',char(row.morph),char(row.recording));end
files=dir(fullfile(folder,'*ALL_CELLS.mat'));names=upper(string({files.name}));
files=files(~contains(names,"SELECTED")&~contains(names,"ANATOMICAL"));
assert(numel(files)==1,'Expected one original ALL_CELLS file in %s; found %d.',folder,numel(files));
A=load(fullfile(files.folder,files.name),'avg','bkg','cells','cell_per');
if isfield(A,'avg'),background=double(A.avg);
elseif isfield(A,'bkg'),background=double(A.bkg);
else,error('ALL_CELLS has no avg or bkg anatomical image.');end
assert(ismatrix(background)&&~isempty(background),'Anatomical image must be a nonempty 2-D matrix.');
if isfield(A,'cells')&&iscell(A.cells),nAll=numel(A.cells);
elseif isfield(A,'cell_per')&&iscell(A.cell_per),nAll=numel(A.cell_per);
else,error('ALL_CELLS has no ROI masks or boundaries.');end
roi=double(roi(:));
assert(size(member,1)==numel(roi)&&numel(unique(roi))==numel(roi)&& ...
 all(isfinite(roi)&roi>=1&roi<=nAll&roi==round(roi)),'Invalid original ROI mapping for anatomical overlay.');
if isfield(X.model,'nTotalRois'),assert(nAll==X.model.nTotalRois,'ALL_CELLS ROI count differs from the assembly model.');end
xy=nan(numel(roi),2);
for j=1:numel(roi)
 r=roi(j);
 if isfield(A,'cells')&&iscell(A.cells)
  idx=double(A.cells{r}(:));
  assert(~isempty(idx)&&all(isfinite(idx)&idx>=1&idx<=numel(background)&idx==round(idx)), ...
   'Invalid mask for original ROI %d.',r);
  [yy,xx]=ind2sub(size(background),idx);xy(j,:)=[mean(xx) mean(yy)];
 else
  boundary=double(A.cell_per{r});
  assert(size(boundary,2)>=2,'Invalid boundary for original ROI %d.',r);
  boundary=boundary(all(isfinite(boundary(:,1:2)),2),1:2);
  assert(~isempty(boundary),'Empty boundary for original ROI %d.',r);
  xy(j,:)=mean(boundary,1);
 end
end
nAssemblies=size(member,2);if nAssemblies==0,return,end
plotDir=fullfile(outDir,'assembly_anatomy');if ~isfolder(plotDir),mkdir(plotDir);end
panelWidth=360;panelHeight=max(240,panelWidth*size(background,1)/size(background,2));
fig=figure('Visible',char(visibility),'Color','w','Position',[50 50 panelWidth*nAssemblies panelHeight+90]);
cleanup=onCleanup(@()close(fig)); %#ok<NASGU>
layout=tiledlayout(fig,1,nAssemblies,'TileSpacing','compact','Padding','compact');
finitePixels=background(isfinite(background));assert(~isempty(finitePixels),'Anatomical image has no finite pixels.');
limits=prctile(finitePixels,[1 99]);if limits(2)<=limits(1),limits=[min(finitePixels) max(finitePixels)];end
if limits(2)<=limits(1),limits=limits(1)+[-.5 .5];end
for a=1:nAssemblies
 ax=nexttile(layout);imagesc(ax,background,limits);colormap(ax,gray(256));hold(ax,'on');
 scatter(ax,xy(member(:,a),1),xy(member(:,a),2),22,[1 1 0],'filled');
 axis(ax,'image');axis(ax,'off');
 title(ax,sprintf('Assembly %d | %d neurons',a,nnz(member(:,a))));
end
title(layout,sprintf('%s / %s',string(row.morph),string(row.recording)),'Interpreter','none');
stem=char(string(row.morph)+"__"+string(row.recording)+"__assembly_anatomy");
style_assembly_figure(fig);exportgraphics(fig,fullfile(plotDir,[stem '.png']),'Resolution',300);
style_assembly_figure(fig);exportgraphics(fig,fullfile(plotDir,[stem '.svg']),'ContentType','vector');

% One anatomical image per fish with a distinct color for each assembly.
overlayFig=figure('Visible',char(visibility),'Color','w','Position',[50 50 1000 750]);
overlayCleanup=onCleanup(@()close(overlayFig)); %#ok<NASGU>
ax=axes(overlayFig);imagesc(ax,background,limits);colormap(ax,gray(256));hold(ax,'on');
assemblyColors=hsv(nAssemblies);
legendHandles=gobjects(nAssemblies,1);legendLabels=strings(nAssemblies,1);
membershipCount=sum(member,2);
remainingMemberships=membershipCount;
for a=1:nAssemblies
 selected=member(:,a);
 exclusive=selected&membershipCount==1;
 shared=selected&membershipCount>1;
 legendHandles(a)=scatter(ax,xy(exclusive,1),xy(exclusive,2),30,assemblyColors(a,:),'filled');
 % Concentric colored rings retain every membership of overlapping neurons.
 if any(shared)
  markerArea=(6+3*(remainingMemberships(shared)-1)).^2;
  scatter(ax,xy(shared,1),xy(shared,2),markerArea,assemblyColors(a,:), ...
   'LineWidth',1.5,'HandleVisibility','off');
 end
 remainingMemberships(selected)=remainingMemberships(selected)-1;
 legendLabels(a)=sprintf('Assembly %d | %d neurons',a,nnz(selected));
end
axis(ax,'image');axis(ax,'off');
title(ax,{sprintf('%s / %s | All assemblies',string(row.morph),string(row.recording)), ...
 'Shared neurons: concentric rings in assembly colors'},'Interpreter','none');
legend(ax,legendHandles,cellstr(legendLabels),'Location','eastoutside','Interpreter','none');
style_assembly_figure(overlayFig);exportgraphics(overlayFig,fullfile(plotDir,[stem '_overlay.png']),'Resolution',300);
style_assembly_figure(overlayFig);exportgraphics(overlayFig,fullfile(plotDir,[stem '_overlay.svg']),'ContentType','vector');

% Truecolor anatomy keeps the phase colormap independent of image intensity.
grayBackground=min(max((background-limits(1))/diff(limits),0),1);
grayBackground(~isfinite(grayBackground))=0;
rgbBackground=repmat(grayBackground,1,1,3);
phaseValues=loadCandidatePhaseFeature(X.S,roi);
phaseFig=figure('Visible',char(visibility),'Color','w','Position',[50 50 1000 750]);
phaseCleanup=onCleanup(@()close(phaseFig)); %#ok<NASGU>
ax=axes(phaseFig);image(ax,rgbBackground);hold(ax,'on');
phaseMap=hsv(256);colormap(ax,phaseMap);clim(ax,[-pi pi]);
remainingMemberships=membershipCount;
hasUndefinedPhase=false;
for a=1:nAssemblies
 selected=member(:,a);
 values=phaseValues(selected&isfinite(phaseValues));
 if isempty(values)
  markerColor=[0.6 0.6 0.6];hasUndefinedPhase=true;
 else
  % Match MATLAB scaled color mapping used by the phase-circle scatter.
  phase=circularMean(values);
  colorIndex=min(256,max(1,1+floor(256*(phase+pi)/(2*pi))));
  markerColor=phaseMap(colorIndex,:);
 end
 exclusive=selected&membershipCount==1;
 shared=selected&membershipCount>1;
 scatter(ax,xy(exclusive,1),xy(exclusive,2),30,markerColor,'filled');
 if any(shared)
  markerArea=(6+3*(remainingMemberships(shared)-1)).^2;
  scatter(ax,xy(shared,1),xy(shared,2),markerArea,markerColor,'LineWidth',1.5);
 end
 remainingMemberships(selected)=remainingMemberships(selected)-1;
end
axis(ax,'image');axis(ax,'off');
cb=colorbar(ax,'eastoutside');cb.Ticks=[-pi -pi/2 0 pi/2 pi];
cb.TickLabels={'-pi','-pi/2','0','pi/2','pi'};
cb.Label.String='Assembly mean preferred phase (rad)';
phaseNote='Shared neurons: concentric rings in assembly phase colors';
if hasUndefinedPhase,phaseNote=[phaseNote '; grey: no finite phase'];end
title(ax,{sprintf('%s / %s | Assembly preferred phases',string(row.morph),string(row.recording)), ...
 phaseNote},'Interpreter','none');
style_assembly_figure(phaseFig);exportgraphics(phaseFig,fullfile(plotDir,[stem '_phase_overlay.png']),'Resolution',300);
style_assembly_figure(phaseFig);exportgraphics(phaseFig,fullfile(plotDir,[stem '_phase_overlay.svg']),'ContentType','vector');

% Every candidate is shown; membership in any assembly is highlighted in red.
candidateFig=figure('Visible',char(visibility),'Color','w','Position',[50 50 1000 750]);
candidateCleanup=onCleanup(@()close(candidateFig)); %#ok<NASGU>
ax=axes(candidateFig);image(ax,rgbBackground);hold(ax,'on');
assigned=any(member,2);
unassignedHandle=scatter(ax,xy(~assigned,1),xy(~assigned,2),30,[0.6 0.6 0.6],'filled');
assignedHandle=scatter(ax,xy(assigned,1),xy(assigned,2),30,[1 0 0],'filled');
axis(ax,'image');axis(ax,'off');
title(ax,sprintf('%s / %s | All candidate neurons',string(row.morph),string(row.recording)), ...
 'Interpreter','none');
legend(ax,[assignedHandle unassignedHandle], ...
 {sprintf('Assigned to an assembly | %d neurons',nnz(assigned)), ...
 sprintf('Unassigned | %d neurons',nnz(~assigned))},'Location','eastoutside');
style_assembly_figure(candidateFig);exportgraphics(candidateFig,fullfile(plotDir,[stem '_candidate_membership_overlay.png']),'Resolution',300);
style_assembly_figure(candidateFig);exportgraphics(candidateFig,fullfile(plotDir,[stem '_candidate_membership_overlay.svg']),'ContentType','vector');

end
