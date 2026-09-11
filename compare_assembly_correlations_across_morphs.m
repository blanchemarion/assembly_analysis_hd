function report=compare_assembly_correlations_across_morphs(batchResults,varargin)
%COMPARE_ASSEMBLY_CORRELATIONS_ACROSS_MORPHS Opposite-phase coupling and turns.
% Fish are independent replicates. Activation episodes use exactly the shared
% defaults and compute_assembly_activations used by the properties comparison.

if nargin<1||isempty(batchResults)
 root=fileparts(mfilename('fullpath'));
 savedBatch=fullfile(root,'data_processed','assemblies','session_manifest.mat');
 assert(isfile(savedBatch),'Saved assembly manifest not found: %s',savedBatch);
 X=load(savedBatch,'manifest','opts');
 assert(isfield(X,'manifest')&&isfield(X,'opts'),'Saved manifest must contain manifest and opts.');
 batchResults=struct('manifest',X.manifest,'options',X.opts);
end
p=inputParser;
addRequired(p,'batchResults',@(x)isstruct(x)&&isfield(x,'manifest')&&isfield(x,'options'));
addParameter(p,'outputDir','',@(x)ischar(x)||isstring(x));
addParameter(p,'phaseDistanceThresholdDeg',135,@(x)isscalar(x)&&isfinite(x)&&x>=0&&x<=180);
addParameter(p,'phaseDistanceBinWidthDeg',30,@(x)isscalar(x)&&isfinite(x)&&x>0&&x<=180);
addParameter(p,'maxLagSeconds',10,@(x)isscalar(x)&&isfinite(x)&&x>0);
addParameter(p,'turnWindowSeconds',10,@(x)isscalar(x)&&isfinite(x)&&x>0);
addParameter(p,'turnAmplitudeThresholdRad',0.239,@(x)isscalar(x)&&isfinite(x)&&x>=0);
addParameter(p,'ahvNativeToPaperSign',-1,@(x)isscalar(x)&&ismember(x,[-1 1]));
addParameter(p,'behaviorTimeMode','auto',@(x)(ischar(x)||isstring(x))&& ...
 ismember(lower(string(x)),["auto","local","startframe","behaviortimesec"]));
addParameter(p,'nOnsetShuffles',1000,@(x)isscalar(x)&&x>=1&&x==round(x));
addParameter(p,'shuffleExclusionSeconds',20,@(x)isscalar(x)&&isfinite(x)&&x>=0);
addParameter(p,'nFishBootstraps',1000,@(x)isscalar(x)&&x>=1&&x==round(x));
addParameter(p,'nLabelPermutations',10000,@(x)isscalar(x)&&x>=1&&x==round(x));
addParameter(p,'minFishPerGroup',3,@(x)isscalar(x)&&x>=2&&x==round(x));
addParameter(p,'randomSeed',1,@(x)isscalar(x)&&isfinite(x)&&x==round(x));
activationDefaults=assembly_activation_defaults();
addParameter(p,'minActiveMembers',activationDefaults.minActiveMembers,@(x)isscalar(x)&&x>=1&&x==round(x));
addParameter(p,'onRecruitmentFraction',activationDefaults.onRecruitmentFraction,@(x)isscalar(x)&&x>=0&&x<=1);
addParameter(p,'offRecruitmentFraction',activationDefaults.offRecruitmentFraction,@(x)isscalar(x)&&x>=0&&x<=1);
addParameter(p,'minMatchingIndex',activationDefaults.minMatchingIndex,@(x)isscalar(x)&&x>=0&&x<=1);
addParameter(p,'activationNullAlpha',activationDefaults.activationNullAlpha,@(x)isscalar(x)&&x>0&&x<=1);
addParameter(p,'nNullShuffles',activationDefaults.nNullShuffles,@(x)isscalar(x)&&x>=1&&x==round(x));
addParameter(p,'minDurationSeconds',activationDefaults.minDurationSeconds,@(x)isscalar(x)&&isfinite(x)&&x>=0);
addParameter(p,'minNullShiftSeconds',activationDefaults.minNullShiftSeconds,@(x)isscalar(x)&&isfinite(x)&&x>=0);
addParameter(p,'activationRandomSeed',activationDefaults.randomSeed,@(x)isscalar(x)&&isfinite(x)&&x==round(x));
addParameter(p,'figureVisible','off',@(x)ischar(x)||isstring(x));
addParameter(p,'pngResolution',300,@(x)isscalar(x)&&x>=150&&x==round(x));
parse(p,batchResults,varargin{:});q=p.Results;q.behaviorTimeMode=lower(string(q.behaviorTimeMode));

opts=batchResults.options;
if strlength(string(q.outputDir))==0
 label="assemblies";if isfield(opts,'analysisLabel'),label=string(opts.analysisLabel);end
 outputDir=fullfile(opts.figureRoot,'assembly_correlations_'+label);
else,outputDir=char(q.outputDir);end
if ~isfolder(outputDir),mkdir(outputDir);end
activationOptions=struct('minActiveMembers',q.minActiveMembers, ...
 'onRecruitmentFraction',q.onRecruitmentFraction,'offRecruitmentFraction',q.offRecruitmentFraction, ...
 'minMatchingIndex',q.minMatchingIndex,'activationNullAlpha',q.activationNullAlpha, ...
 'nNullShuffles',q.nNullShuffles,'minDurationSeconds',q.minDurationSeconds, ...
 'minNullShiftSeconds',q.minNullShiftSeconds,'randomSeed',q.activationRandomSeed);

M=batchResults.manifest;
required=["status","outputPath","morph","recording"];
assert(all(ismember(required,string(M.Properties.VariableNames))), ...
 'Manifest must contain status, outputPath, morph, and recording.');
M=M(ismember(string(M.status),["processed","existing"]),:);
if isempty(M),error('assembly:NoBatchOutputs','No processed/existing recordings in the manifest.');end

pairbots=emptyPairTable();pairCurves=emptyPairCurveTable();phaseDistancePairs=emptyPhaseDistancePairTable();fishCorr=struct([]);
turnEvents=emptyTurnEventTable();turnFish=struct([]);sessionRows=struct([]);
activationAudits=struct('morph',{},'fish',{},'file',{},'audit',{});
shiftRows=struct([]);failRows=struct([]);
zeroAssemblyRows=struct('morph',{},'fish',{},'file',{},'reason',{});
stream=RandStream('mt19937ar','Seed',q.randomSeed);
fprintf('\n=== Assembly correlations and turn timing: %d files ===\n',height(M));
for s=1:height(M)
 morph=lower(string(M.morph(s)));fish=string(M.recording(s));file=char(M.outputPath(s));
 session=struct('morph',morph,'fish',fish,'file',string(file),'status',"excluded", ...
  'reason',"",'nAssemblies',0,'nOppositePhasePairs',0,'nIncludedPairs',0, ...
  'nOverlappingPairsExcluded',0,'nAllTurns',0,'nPositiveTurns',0,'nNegativeTurns',0, ...
  'correlationStatus',"not_run",'turnStatus',"not_run",'behaviorPath',"",'imagingTimeSource',"");
 fprintf('[%d/%d] %s/%s\n',s,height(M),morph,fish);
 try
  X=load(file,'model','S');requireFields(X,{'model','S'},'assembly output');
  requireFields(X.model,{'membership','candidateRoiIdx','candidateCellIDs'},'model');
  requireFields(X.S,{'pca','activity','fps','preferredPhaseRad','candidateRoiIdx','candidateCellIDs'},'S');
  member=logical(X.model.membership);F=double(X.S.pca);activity=double(X.S.activity);
  fps=double(X.S.fps);nFrames=size(F,1);nAssemblies=size(member,2);
  assert(isequal(size(F),size(activity))&&size(F,2)==size(member,1),'Trace/raster/membership dimensions disagree.');
  assert(isequal(double(X.model.candidateRoiIdx(:)),double(X.S.candidateRoiIdx(:)))&& ...
   isequal(double(X.model.candidateCellIDs(:)),double(X.S.candidateCellIDs(:))), ...
   'Candidate identities are not aligned with membership rows.');
  assert(isscalar(fps)&&isfinite(fps)&&fps>0,'Invalid imaging fps.');
  session.nAssemblies=nAssemblies;
  if nAssemblies==0
   session.reason="zero assemblies detected";session.correlationStatus="excluded_zero_assemblies";session.turnStatus="excluded_zero_assemblies";
   zeroAssemblyRows(end+1)=struct('morph',morph,'fish',fish,'file',string(file), ...
    'reason',"zero assemblies detected"); %#ok<AGROW>
   fprintf('  Excluded from comparison: zero assemblies detected.\n');
   sessionRows=appendStruct(sessionRows,session);continue
  end
  phases=double(X.S.preferredPhaseRad(:));
  assert(numel(phases)==size(member,1)&&all(phases(isfinite(phases))>-pi&phases(isfinite(phases))<=pi), ...
   'Preferred phases are missing, misaligned, or outside (-pi,pi].');
  A=compute_assembly_activations(member,activity,fps,activationOptions);
  activationAudits(end+1)=struct('morph',morph,'fish',fish,'file',string(file),'audit',A); %#ok<AGROW>
  onsets=false(nAssemblies,nFrames);
  for a=1:nAssemblies,onsets(a,A.onsetIndices{a})=true;end
  [tCa,timeSource]=loadImagingTime(M(s,:),nFrames);
  session.imagingTimeSource=timeSource;
  segments=continuousSegments(tCa);

  [P,C,D,fishRecord]=assemblyPairAnalysis(member,phases,F,tCa,segments,morph,fish,q);
  pairbots=[pairbots;P];pairCurves=[pairCurves;C];phaseDistancePairs=[phaseDistancePairs;D]; %#ok<AGROW>
  fishCorr=appendStruct(fishCorr,fishRecord);
  session.nOppositePhasePairs=sum(P.phaseEligible);
  session.nIncludedPairs=sum(P.inPrimaryAnalysis);
  session.nOverlappingPairsExcluded=sum(P.exclusionReason=="overlapping_membership");
  if fishRecord.eligible,session.correlationStatus="included";else,session.correlationStatus="insufficient_data";end

  try
   behaviorPath=findBehaviorFile(M(s,:));
   session.behaviorPath=behaviorPath;
   [events,turnInfo]=loadTurns(behaviorPath,tCa,morph,fish,q);
   turnEvents=[turnEvents;events]; %#ok<AGROW>
   session.nAllTurns=turnInfo.nAll;session.nPositiveTurns=turnInfo.nPositive;session.nNegativeTurns=turnInfo.nNegative;
   [records,shiftRecord]=turnOnsetAnalysis(onsets,tCa,segments,events,morph,fish,q,stream);
   for k=1:numel(records),turnFish=appendStruct(turnFish,records(k));end
   shiftRows=appendStruct(shiftRows,shiftRecord);
   if any([records.eligible]),session.turnStatus="included";else,session.turnStatus="insufficient_data";end
  catch turnError
   session.turnStatus="excluded";session.reason=joinReason(session.reason,"turn analysis: "+string(turnError.message));
   warning('assembly:TurnAnalysisExcluded','%s/%s: %s',morph,fish,turnError.message);
  end
  session.status="loaded";
 catch ME
  session.reason=joinReason(session.reason,string(ME.message));
  failRows=appendStruct(failRows,struct('morph',morph,'fish',fish,'file',string(file),'message',string(ME.message)));
  warning('assembly:CorrelationAnalysisFailed','%s/%s: %s',morph,fish,ME.message);
 end
 sessionRows=appendStruct(sessionRows,session);
end

sessionTable=struct2table(sessionRows);failureTable=structToTableOrEmpty(failRows, ...
 {'morph','fish','file','message'},{'string','string','string','string'});
zeroAssemblyTable=structToTableOrEmpty(zeroAssemblyRows, ...
 {'morph','fish','file','reason'},{'string','string','string','string'});
fishCorrTable=fishCorrelationTable(fishCorr);
fishTurnTable=fishTurnScalarTable(turnFish);
fishScalars=outerjoin(fishCorrTable(:,{'morph','fish','zeroLagCorrelation','nPairs'}), ...
 pivotTurnScalars(fishTurnTable),'Keys',{'morph','fish'},'MergeKeys',true,'Type','full');
[corrFishCurves,corrMorphCurves]=summarizeCorrelationCurves(fishCorr,q,stream);
[phaseDistanceFishBins,phaseDistanceMorphBins]=summarizePhaseDistance(phaseDistancePairs,q,stream);
[turnFishCurves,turnMorphCurves]=summarizeTurnCurves(turnFish,q,stream);
statistics=scalarStatistics(fishScalars,q,stream);
morphScalarSummary=morphScalars(fishScalars,q,stream);
turnPeakSummary=morphTurnPeaks(turnMorphCurves);
shiftTable=shiftRecordTable(shiftRows);
parameterTable=parametersTable(q,activationOptions);

writetable(sessionTable,fullfile(outputDir,'session_summary.csv'));
writetable(failureTable,fullfile(outputDir,'processing_failures.csv'));
writetable(zeroAssemblyTable,fullfile(outputDir,'zero_assembly_recordings_excluded.csv'));
writetable(pairbots,fullfile(outputDir,'opposite_phase_pair_exclusions_and_scalars.csv'));
writetable(pairCurves,fullfile(outputDir,'opposite_phase_pair_cross_correlations.csv'));
writetable(fishCorrTable,fullfile(outputDir,'fish_cross_correlation_scalars.csv'));
writetable(corrFishCurves,fullfile(outputDir,'fish_cross_correlation_curves.csv'));
writetable(corrMorphCurves,fullfile(outputDir,'morph_cross_correlation_curves.csv'));
writetable(phaseDistancePairs,fullfile(outputDir,'phase_distance_zero_lag_pairs.csv'));
writetable(phaseDistanceFishBins,fullfile(outputDir,'phase_distance_zero_lag_fish_bins.csv'));
writetable(phaseDistanceMorphBins,fullfile(outputDir,'phase_distance_zero_lag_morph_bins.csv'));
writetable(turnEvents,fullfile(outputDir,'turn_event_alignment_audit.csv'));
writetable(fishTurnTable,fullfile(outputDir,'fish_turn_onset_scalars.csv'));
writetable(turnFishCurves,fullfile(outputDir,'fish_turn_onset_curves.csv'));
writetable(turnMorphCurves,fullfile(outputDir,'morph_turn_onset_curves.csv'));
writetable(fishScalars,fullfile(outputDir,'fish_primary_scalars.csv'));
writetable(morphScalarSummary,fullfile(outputDir,'morph_scalar_summary.csv'));
writetable(turnPeakSummary,fullfile(outputDir,'morph_turn_peak_lags.csv'));
writetable(statistics,fullfile(outputDir,'fish_level_statistics.csv'));
writetable(shiftTable,fullfile(outputDir,'onset_shuffle_offsets.csv'));
writetable(parameterTable,fullfile(outputDir,'analysis_parameters.csv'));

plotCurveSummary(corrMorphCurves,'meanCorrelation','Assembly cross-correlation','Lag (s)', ...
 'Opposite-phase assembly cross-correlations (symmetrized pairs)', ...
 fullfile(outputDir,'opposite_phase_assembly_cross_correlations'),q);
plotPhaseDistance(phaseDistanceMorphBins,fullfile(outputDir,'zero_lag_correlation_by_phase_distance'),q);
plotScalar(fishScalars,'zeroLagCorrelation','Zero-lag Pearson correlation', ...
 'Opposite-phase assembly correlation: fish-level zero lag', ...
 fullfile(outputDir,'opposite_phase_zero_lag_fish_scalars'),q,stream,statistics);
plotTurnCurves(turnMorphCurves,outputDir,q);
for category=["all","positive","negative"]
 field=char("turnExcess_"+category);
 plotScalar(fishScalars,field,'Mean excess onset rate (onsets/s/assembly)', ...
  char("Turn-triggered excess onset rate: "+category+" turns, -3 to +3 s"), ...
  fullfile(outputDir,char("turn_excess_"+category+"_fish_scalars")),q,stream,statistics);
end

report=struct('parameters',q,'activationOptions',activationOptions,'manifest',M, ...
 'sessions',sessionTable,'failures',failureTable,'zeroAssemblyRecordingsExcluded',zeroAssemblyTable, ...
 'activationAudits',{activationAudits}, ...
 'pairAudit',pairbots,'pairCurves',pairCurves,'fishCorrelationScalars',fishCorrTable, ...
 'phaseDistancePairs',phaseDistancePairs,'phaseDistanceFishBins',phaseDistanceFishBins, ...
 'phaseDistanceMorphBins',phaseDistanceMorphBins, ...
 'fishCorrelationCurves',corrFishCurves,'morphCorrelationCurves',corrMorphCurves, ...
 'turnEvents',turnEvents,'fishTurnScalars',fishTurnTable,'fishTurnCurves',turnFishCurves, ...
 'morphTurnCurves',turnMorphCurves,'fishPrimaryScalars',fishScalars,'statistics',statistics, ...
 'morphScalarSummary',morphScalarSummary,'morphTurnPeakLags',turnPeakSummary, ...
 'shuffleOffsets',shiftTable,'outputDir',string(outputDir));
save(fullfile(outputDir,'assembly_correlations_across_morphs_results.mat'),'report','-v7.3');
fprintf('Assembly correlation comparison complete. Outputs: %s\n',outputDir);
end

function [pairTable,curveTable,distanceTable,fishRecord]=assemblyPairAnalysis(member,phase,F,t,segments,morph,fish,q)
nA=size(member,2);lags=(-round(q.maxLagSeconds/(median(diff(t)))):round(q.maxLagSeconds/(median(diff(t)))))';
lagSec=lags*median(diff(t));assemblyPhase=nan(nA,1);
for a=1:nA
 v=phase(member(:,a)&isfinite(phase));
 if ~isempty(v),assemblyPhase(a)=atan2(mean(sin(v)),mean(cos(v)));end
end
signals=nan(size(F,1),nA);
for a=1:nA,signals(:,a)=mean(F(:,member(:,a)),2,'omitmissing');end
pairTable=emptyPairTable();curveTable=emptyPairCurveTable();distanceTable=emptyPhaseDistancePairTable();curves=zeros(0,numel(lags));
for a=1:nA-1
 for b=a+1:nA
  d=rad2deg(abs(atan2(sin(assemblyPhase(a)-assemblyPhase(b)),cos(assemblyPhase(a)-assemblyPhase(b)))));
  eligible=isfinite(d)&&d>=q.phaseDistanceThresholdDeg;shared=nnz(member(:,a)&member(:,b));
  distanceReason="";distanceIncluded=false;distanceZ=NaN;distanceN=0;
  if ~isfinite(d),distanceReason="missing_assembly_phase";
  elseif shared>0,distanceReason="overlapping_membership";
  else
   [distanceZ,distanceN]=laggedCorrelation(signals(:,a),signals(:,b),segments,0);
   if isfinite(distanceZ),distanceIncluded=true;else,distanceReason="invalid_zero_lag_correlation";end
  end
  distanceTable=[distanceTable;table(morph,fish,a,b,assemblyPhase(a),assemblyPhase(b),d,shared, ...
   distanceIncluded,distanceZ,distanceN,distanceReason,'VariableNames',distanceTable.Properties.VariableNames)]; %#ok<AGROW>
  reason="";primary=false;z=NaN;validZero=0;c=nan(size(lags));
  if ~isfinite(d),reason="missing_assembly_phase";
  elseif ~eligible,reason="below_phase_distance_threshold";
  elseif shared>0,reason="overlapping_membership";
  else
   [raw,nValid]=laggedCorrelation(signals(:,a),signals(:,b),segments,lags);
   c=mean([raw,flipud(raw)],2,'omitmissing');validZero=nValid(lags==0);z=c(lags==0);
   if ~isfinite(z),reason="invalid_zero_lag_correlation";
   else,primary=true;curves(end+1,:)=c';end %#ok<AGROW>
  end
  row=table(morph,fish,a,b,assemblyPhase(a),assemblyPhase(b),d,shared,eligible,primary, ...
   z,validZero,reason,'VariableNames',pairTable.Properties.VariableNames);
  pairTable=[pairTable;row]; %#ok<AGROW>
  if primary
   n=numel(lags);curveTable=[curveTable;table(repmat(morph,n,1),repmat(fish,n,1), ...
    repmat(a,n,1),repmat(b,n,1),lags,lagSec,c, ...
    'VariableNames',curveTable.Properties.VariableNames)]; %#ok<AGROW>
  end
 end
end
mu=nan(size(lagSec));if ~isempty(curves),mu=mean(curves,1,'omitmissing')';end
fishRecord=struct('morph',morph,'fish',fish,'lagSeconds',lagSec,'correlation',mu, ...
 'zeroLagCorrelation',mu(lags==0),'nPairs',size(curves,1),'eligible',~isempty(curves));
end

function [r,nValid]=laggedCorrelation(x,y,segments,lags)
r=nan(size(lags));nValid=zeros(size(lags));n=numel(x);
for j=1:numel(lags)
 k=lags(j);
 if k>=0,i=(1:n-k)';h=i+k;else,h=(1:n+k)';i=h-k;end
 ok=segments(i)==segments(h)&segments(i)>0&isfinite(x(i))&isfinite(y(h));
 xx=x(i(ok));yy=y(h(ok));nValid(j)=numel(xx);
 if numel(xx)>=3&&std(xx)>0&&std(yy)>0,C=corrcoef(xx,yy);r(j)=C(1,2);end
end
end

function [events,info]=loadTurns(path,tCa,morph,fish,q)
X=load(path,'pass2FileResult');assert(isfield(X,'pass2FileResult'),'Behavior file lacks pass2FileResult.');
B=X.pass2FileResult;assert(isfield(B,'fps')&&isfinite(double(B.fps(1)))&&double(B.fps(1))>0,'Invalid behavior fps.');
nBeh=inferBehaviorLength(B);assert(nBeh>=2,'Could not infer behavior length.');
[tb,timeSource]=chooseBehaviorTime(B,nBeh,double(B.fps(1)),tCa,q.behaviorTimeMode);
if isfield(B,'dtheta')&&~isempty(B.dtheta),angle=cleanAngleUnitsToRad(B.dtheta);source="dtheta";
elseif isfield(B,'LI')&&~isempty(B.LI),angle=cleanAngleUnitsToRad(B.LI);source="LI fallback";
else,error('Behavior file contains neither dtheta nor LI.');end
assert(isfield(B,'swimFrames')&&~isempty(B.swimFrames),'swimFrames is required.');
[startFrame,~]=parseSwimFrames(double(B.swimFrames),numel(angle));
n=min(numel(angle),numel(startFrame));angle=angle(1:n);startFrame=round(startFrame(1:n));
paper=q.ahvNativeToPaperSign*angle;
valid=isfinite(startFrame)&startFrame>=1&startFrame<=nBeh&isfinite(paper)&abs(paper)>q.turnAmplitudeThresholdRad;
times=nan(n,1);times(valid)=tb(startFrame(valid));
idx=round(interp1(tCa(:),(1:numel(tCa))',times,'linear',NaN));
included=valid&isfinite(idx)&idx>=1&idx<=numel(tCa);
reason=repmat("below_threshold_or_invalid",n,1);reason(valid)="outside_imaging_support";reason(included)="included";
direction=repmat("excluded",n,1);direction(included&paper>0)="positive";direction(included&paper<0)="negative";
events=table(repmat(morph,n,1),repmat(fish,n,1),(1:n)',startFrame,angle,paper,times,idx, ...
 included,direction,reason,repmat(source,n,1),repmat(string(timeSource),n,1), ...
 'VariableNames',emptyTurnEventTable().Properties.VariableNames);
info=struct('nAll',nnz(included),'nPositive',nnz(included&paper>0),'nNegative',nnz(included&paper<0));
end

function [records,shiftRecord]=turnOnsetAnalysis(onsets,t,segments,events,morph,fish,q,stream)
dt=median(diff(t));lags=(-round(q.turnWindowSeconds/dt):round(q.turnWindowSeconds/dt))';
lagSec=lags*dt;pooled=mean(onsets,1);included=events.included;idx=events.imagingIndex;
groups={idx(included),idx(included&events.direction=="positive"),idx(included&events.direction=="negative")};
names=["all","positive","negative"];observed=nan(3,numel(lags));exposure=zeros(3,numel(lags));
for c=1:3,[observed(c,:),exposure(c,:)]=triggeredRate(pooled,groups{c},segments,lags,dt);end
safe=max(1,ceil(q.shuffleExclusionSeconds/dt));allowed=safe:(numel(t)-safe);
allowed=allowed(allowed>0&allowed<numel(t));
nullSum=zeros(size(observed));nullN=zeros(size(observed));used=zeros(0,1);
if ~isempty(allowed)
 if numel(allowed)>=q.nOnsetShuffles
  used=allowed(randperm(stream,numel(allowed),q.nOnsetShuffles));
 else,used=allowed(randi(stream,numel(allowed),q.nOnsetShuffles,1));end
 for b=1:q.nOnsetShuffles
  shifted=circshift(pooled,[0 used(b)]);
  for c=1:3
   rate=triggeredRate(shifted,groups{c},segments,lags,dt);
   ok=isfinite(rate);nullSum(c,ok)=nullSum(c,ok)+rate(ok);nullN(c,ok)=nullN(c,ok)+1;
  end
 end
end
nullMean=nullSum./nullN;nullMean(nullN==0)=NaN;corrected=observed-nullMean;
records=struct([]);primary=lagSec>=-3&lagSec<=3;
for c=1:3
 curve=corrected(c,:);finitePrimary=isfinite(curve(primary));
 scalar=NaN;if any(finitePrimary),v=curve(primary);scalar=mean(v(isfinite(v)));end
 peak=NaN;if any(isfinite(curve)),[~,j]=max(curve,[],'omitmissing');peak=lagSec(j);end
 record=struct('morph',morph,'fish',fish,'category',names(c),'lagSeconds',lagSec, ...
  'observedRate',observed(c,:)','shuffleMeanRate',nullMean(c,:)','correctedRate',curve', ...
  'validTurnExposure',exposure(c,:)','nTurns',numel(groups{c}),'meanExcessMinus3ToPlus3',scalar, ...
  'peakLagSeconds',peak,'nShuffles',numel(used),'eligible',numel(groups{c})>0&&isfinite(scalar)&&~isempty(allowed));
 records=appendStruct(records,record);
end
shiftRecord=struct('morph',morph,'fish',fish,'nAllowedOffsets',numel(allowed), ...
 'nShuffles',numel(used),'offsetFramesJson',string(jsonencode(used(:)')), ...
 'offsetSecondsJson',string(jsonencode((used(:)'*dt))));
end

function [rate,exposure]=triggeredRate(signal,turnIdx,segments,lags,dt)
rate=nan(1,numel(lags));exposure=zeros(1,numel(lags));turnIdx=round(turnIdx(isfinite(turnIdx)));
for j=1:numel(lags)
 target=turnIdx+lags(j);ok=target>=1&target<=numel(signal);
 target=target(ok);source=turnIdx(ok);ok=segments(target)==segments(source)&segments(target)>0;
 target=target(ok);exposure(j)=numel(target);
 if exposure(j)>0,rate(j)=sum(signal(target),'omitmissing')/exposure(j)/dt;end
end
end

function [fishTable,morphTable]=summarizePhaseDistance(P,q,stream)
fishTable=table(strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1), ...
 nan(0,1),zeros(0,1),'VariableNames',{'morph','fish','binLowerDeg','binUpperDeg', ...
 'binCenterDeg','meanZeroLagCorrelation','nPairs'});
morphTable=table(strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1),nan(0,1),nan(0,1), ...
 nan(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'morph','binLowerDeg','binUpperDeg', ...
 'binCenterDeg','meanZeroLagCorrelation','ciLower','ciUpper','nFish','nPairs'});
edges=0:q.phaseDistanceBinWidthDeg:180;if edges(end)<180,edges=[edges 180];end
P=P(P.included&isfinite(P.zeroLagCorrelation),:);
for e=1:numel(edges)-1
 if e==numel(edges)-1,pick=P.circularPhaseDistanceDeg>=edges(e)&P.circularPhaseDistanceDeg<=edges(e+1);
 else,pick=P.circularPhaseDistanceDeg>=edges(e)&P.circularPhaseDistanceDeg<edges(e+1);end
 U=P(pick,:);if isempty(U),continue,end
 [G,morph,fish]=findgroups(U.morph,U.fish);
 mu=splitapply(@(x)mean(x,'omitmissing'),U.zeroLagCorrelation,G);
 np=splitapply(@numel,U.zeroLagCorrelation,G);n=numel(mu);center=(edges(e)+edges(e+1))/2;
 fishTable=[fishTable;table(morph,fish,repmat(edges(e),n,1),repmat(edges(e+1),n,1), ...
  repmat(center,n,1),mu,np,'VariableNames',fishTable.Properties.VariableNames)]; %#ok<AGROW>
end
for morph=["surface","molino","pachon"]
 for e=1:numel(edges)-1
  U=fishTable(fishTable.morph==morph&fishTable.binLowerDeg==edges(e),:);if isempty(U),continue,end
  [lo,hi]=bootstrapCI(U.meanZeroLagCorrelation,q.nFishBootstraps,stream);
  morphTable=[morphTable;table(morph,edges(e),edges(e+1),(edges(e)+edges(e+1))/2, ...
   mean(U.meanZeroLagCorrelation,'omitmissing'),lo,hi,height(U),sum(U.nPairs), ...
   'VariableNames',morphTable.Properties.VariableNames)]; %#ok<AGROW>
 end
end
end

function [fishTable,morphTable]=summarizeCorrelationCurves(records,q,stream)
fishTable=table();valid=find(arrayfun(@(x)x.eligible,records));
for j=valid(:)'
 r=records(j);n=numel(r.lagSeconds);fishTable=[fishTable;table(repmat(r.morph,n,1), ...
  repmat(r.fish,n,1),r.lagSeconds,r.correlation,repmat(r.nPairs,n,1), ...
  'VariableNames',{'morph','fish','lagSeconds','correlation','nPairs'})]; %#ok<AGROW>
end
morphTable=summarizeCurveRecords(records,'correlation',q.maxLagSeconds,q.nFishBootstraps,stream,"meanCorrelation");
end

function [fishTable,morphTable]=summarizeTurnCurves(records,q,stream)
fishTable=table();
for j=1:numel(records)
 r=records(j);n=numel(r.lagSeconds);
 fishTable=[fishTable;table(repmat(r.morph,n,1),repmat(r.fish,n,1),repmat(r.category,n,1), ...
  r.lagSeconds,r.observedRate,r.shuffleMeanRate,r.correctedRate,r.validTurnExposure, ...
  repmat(r.nTurns,n,1),repmat(r.eligible,n,1), ...
  'VariableNames',{'morph','fish','category','lagSeconds','observedRate','shuffleMeanRate', ...
  'correctedRate','validTurnExposure','nTurns','eligible'})]; %#ok<AGROW>
end
morphTable=table();
for category=["all","positive","negative"]
 pick=arrayfun(@(x)x.category==category&&x.eligible,records);
 for spec={["observedRate","observed"],["correctedRate","corrected"]}
  T=summarizeCurveRecords(records(pick),char(spec{1}(1)),q.turnWindowSeconds,q.nFishBootstraps,stream,"meanRate");
  if ~isempty(T),T.category=repmat(category,height(T),1);T.signal=repmat(spec{1}(2),height(T),1);morphTable=[morphTable;T];end %#ok<AGROW>
 end
end
if ~isempty(morphTable),morphTable=movevars(morphTable,{'category','signal'},'After','morph');end
end

function T=summarizeCurveRecords(records,field,maxSeconds,nBoot,stream,meanName)
T=table();order=["surface","molino","pachon"];
for m=1:3
 R=records(arrayfun(@(x)x.morph==order(m)&&x.eligible,records));if isempty(R),continue,end
 fps=arrayfun(@(x)1/median(diff(x.lagSeconds)),R);[~,ref]=min(abs(fps-median(fps)));
 grid=R(ref).lagSeconds(:);grid=grid(abs(grid)<=maxSeconds+eps);
 Y=nan(numel(R),numel(grid));
 for f=1:numel(R),Y(f,:)=interp1(R(f).lagSeconds,double(R(f).(field)),grid,'linear',NaN);end
 mu=mean(Y,1,'omitmissing');lo=nan(size(mu));hi=lo;
 for k=1:numel(grid),[lo(k),hi(k)]=bootstrapCI(Y(:,k),nBoot,stream);end
 nFish=sum(isfinite(Y),1)';
 U=table(repmat(order(m),numel(grid),1),grid,mu',lo',hi',nFish, ...
  'VariableNames',{'morph','lagSeconds',char(meanName),'ciLower','ciUpper','nFish'});
 T=[T;U]; %#ok<AGROW>
end
end

function T=fishCorrelationTable(R)
T=table(strings(0,1),strings(0,1),nan(0,1),zeros(0,1),false(0,1), ...
 'VariableNames',{'morph','fish','zeroLagCorrelation','nPairs','eligible'});
for j=1:numel(R),T=[T;table(R(j).morph,R(j).fish,R(j).zeroLagCorrelation,R(j).nPairs,R(j).eligible, ...
 'VariableNames',T.Properties.VariableNames)];end
end

function T=fishTurnScalarTable(R)
T=table(strings(0,1),strings(0,1),strings(0,1),zeros(0,1),nan(0,1),nan(0,1),false(0,1), ...
 'VariableNames',{'morph','fish','category','nTurns','meanExcessMinus3ToPlus3','peakLagSeconds','eligible'});
for j=1:numel(R),T=[T;table(R(j).morph,R(j).fish,R(j).category,R(j).nTurns, ...
 R(j).meanExcessMinus3ToPlus3,R(j).peakLagSeconds,R(j).eligible,'VariableNames',T.Properties.VariableNames)];end
end

function W=pivotTurnScalars(T)
keys=unique(T(:,{'morph','fish'}));W=keys;
for category=["all","positive","negative"]
 U=T(T.category==category,{'morph','fish','meanExcessMinus3ToPlus3'});
 U.Properties.VariableNames{3}=char("turnExcess_"+category);W=outerjoin(W,U,'Keys',{'morph','fish'},'MergeKeys',true,'Type','full');
end
end

function T=scalarStatistics(F,q,stream)
metrics=["zeroLagCorrelation","turnExcess_all","turnExcess_positive","turnExcess_negative"];
contrasts=["surface","molino";"surface","pachon"];rows=struct([]);
for metric=metrics
 temp=struct([]);
 for c=1:2
  a=double(F.(metric)(F.morph==contrasts(c,1)));a=a(isfinite(a));
  b=double(F.(metric)(F.morph==contrasts(c,2)));b=b(isfinite(b));
  status="included";pPerm=NaN;pKW=NaN;stat=NaN;kwStat=NaN;kwDf=NaN;
  if numel(a)<q.minFishPerGroup||numel(b)<q.minFishPerGroup,status="insufficient_fish";
  else
   stat=mean(a)-mean(b);pool=[a;b];nA=numel(a);null=zeros(q.nLabelPermutations,1);
   for k=1:q.nLabelPermutations,ix=randperm(stream,numel(pool));null(k)=mean(pool(ix(1:nA)))-mean(pool(ix(nA+1:end)));end
   pPerm=(1+sum(abs(null)>=abs(stat)))/(q.nLabelPermutations+1);
   try
    pairGroups=[repmat(contrasts(c,1),numel(a),1);repmat(contrasts(c,2),numel(b),1)];
    [pKW,K]=kruskalwallis(pool,cellstr(pairGroups),'off');
    kwDf=K{2,3};kwStat=K{2,5};
   catch,pKW=NaN;kwDf=NaN;end
  end
  temp=appendStruct(temp,struct('metric',metric,'prespecifiedPrimary',metric=="zeroLagCorrelation"||metric=="turnExcess_all", ...
   'group1',contrasts(c,1),'group2',contrasts(c,2),'n1Fish',numel(a),'n2Fish',numel(b), ...
   'mean1',safeMean(a),'mean2',safeMean(b),'meanDifference',stat, ...
   'median1',safeMedian(a),'median2',safeMedian(b),'permutationP',pPerm, ...
   'permutationPAdjusted',NaN,'kwTest',"two-group Kruskal-Wallis", ...
   'kwStatistic',kwStat,'kwDegreesFreedom',kwDf,'kwP',pKW,'kwPAdjusted',NaN, ...
   'degreesOfFreedom',kwDf,'pValue',pKW, ...
   'multipleTestingCorrection',"none",'status',status));
 end
 for c=1:2
  temp(c).permutationPAdjusted=temp(c).permutationP;
  temp(c).kwPAdjusted=temp(c).kwP;
  rows=appendStruct(rows,temp(c));
 end
end
T=struct2table(rows);
end

function T=morphScalars(F,q,stream)
metrics=["zeroLagCorrelation","turnExcess_all","turnExcess_positive","turnExcess_negative"];
order=["surface","molino","pachon"];rows=struct([]);
for metric=metrics
 for morph=order
  v=double(F.(metric)(F.morph==morph));v=v(isfinite(v));[lo,hi]=bootstrapCI(v,q.nFishBootstraps,stream);
  rows=appendStruct(rows,struct('metric',metric,'morph',morph,'nFish',numel(v), ...
   'mean',safeMean(v),'median',safeMedian(v),'bootstrapCILower',lo,'bootstrapCIUpper',hi, ...
   'status',string(ternary(~isempty(v),"included","insufficient_data"))));
 end
end
T=struct2table(rows);
end

function T=morphTurnPeaks(C)
T=table(strings(0,1),strings(0,1),strings(0,1),zeros(0,1),nan(0,1),nan(0,1), ...
 'VariableNames',{'morph','category','signal','nFishAtPeak','peakLagSeconds','peakMeanRate'});
if isempty(C),return,end
for morph=["surface","molino","pachon"]
 for category=["all","positive","negative"]
  for signal=["observed","corrected"]
   U=C(C.morph==morph&C.category==category&C.signal==signal,:);if isempty(U),continue,end
   [peak,j]=max(U.meanRate,[],'omitmissing');T=[T;table(morph,category,signal,U.nFish(j),U.lagSeconds(j),peak, ...
    'VariableNames',T.Properties.VariableNames)]; %#ok<AGROW>
  end
 end
end
end

function plotPhaseDistance(T,stem,q)
fig=figure('Visible',char(q.figureVisible),'Color','w','Position',[100 100 820 530]);ax=axes(fig);hold(ax,'on');
colors=morphColors();order=["surface","molino","pachon"];
if isempty(T),text(ax,.5,.5,'Insufficient data','Units','normalized','HorizontalAlignment','center');
else
 for m=1:3
  U=T(T.morph==order(m),:);if isempty(U),continue,end
  edges=0:q.phaseDistanceBinWidthDeg:180;if edges(end)<180,edges=[edges 180];end
  centers=(edges(1:end-1)+edges(2:end))/2;
  y=nan(size(centers));lo=y;hi=y;[tf,loc]=ismember(U.binCenterDeg,centers);
  y(loc(tf))=U.meanZeroLagCorrelation(tf);lo(loc(tf))=U.ciLower(tf);hi(loc(tf))=U.ciUpper(tf);
  fillBandWithGaps(ax,centers,lo,hi,colors(m,:));
  labels=arrayfun(@(nf,np)sprintf('%d fish, %d pairs',nf,np),U.nFish,U.nPairs,'UniformOutput',false);
  plot(ax,centers,y,'-o','LineWidth',2.5,'Color',colors(m,:), ...
   'MarkerFaceColor',colors(m,:),'DisplayName',char(order(m)),'ButtonDownFcn',[]);
  for k=1:height(U),text(ax,U.binCenterDeg(k),U.meanZeroLagCorrelation(k),['  ' labels{k}], ...
    'Color',colors(m,:),'FontSize',7,'Rotation',25,'Clipping','on');end
 end
end
xlim(ax,[0 180]);ticks=0:q.phaseDistanceBinWidthDeg:180;if ticks(end)<180,ticks=[ticks 180];end
xticks(ax,ticks);yline(ax,0,':k','HandleVisibility','off');
xlabel(ax,'Assembly preferred-phase distance (deg)');ylabel(ax,'Zero-lag Pearson correlation');
title(ax,'Assembly correlation versus preferred-phase distance | fish-weighted mean and 95% CI');
if ~isempty(T),legend(ax,'Location','best');end
grid(ax,'on');box(ax,'off');saveFigure(fig,stem,q.pngResolution);close(fig);
end

function fillBandWithGaps(ax,x,lo,hi,color)
ok=isfinite(lo)&isfinite(hi);starts=find(ok&[true ~ok(1:end-1)]);stops=find(ok&[~ok(2:end) true]);
for j=1:numel(starts)
 ix=starts(j):stops(j);
 if numel(ix)==1
  dx=.35*median(diff(x));if ~isfinite(dx),dx=1;end
  xx=[x(ix)-dx x(ix)+dx];ll=[lo(ix) lo(ix)];hh=[hi(ix) hi(ix)];
 else,xx=x(ix);ll=lo(ix);hh=hi(ix);end
 fill(ax,[xx fliplr(xx)],[ll fliplr(hh)],color,'FaceAlpha',.13,'EdgeColor','none','HandleVisibility','off');
end
end

function plotCurveSummary(T,field,yLabel,xLabel,titleText,stem,q)
fig=figure('Visible',char(q.figureVisible),'Color','w','Position',[100 100 820 530]);ax=axes(fig);hold(ax,'on');
colors=morphColors();order=["surface","molino","pachon"];
if isempty(T)||~ismember('morph',T.Properties.VariableNames)
 text(ax,.5,.5,'Insufficient data','Units','normalized','HorizontalAlignment','center');
else
 for m=1:3
  U=T(T.morph==order(m),:);if isempty(U),continue,end
  fill(ax,[U.lagSeconds;flipud(U.lagSeconds)],[U.ciLower;flipud(U.ciUpper)],colors(m,:), ...
   'FaceAlpha',.13,'EdgeColor','none','HandleVisibility','off');
  plot(ax,U.lagSeconds,U.(field),'LineWidth',2.5,'Color',colors(m,:),'DisplayName',sprintf('%s (max n=%d fish)',order(m),max(U.nFish)));
 end
end
xline(ax,0,'--k','HandleVisibility','off');yline(ax,0,':k','HandleVisibility','off');xlabel(ax,xLabel);ylabel(ax,yLabel);
title(ax,titleText);if ~isempty(T),legend(ax,'Location','best');end
grid(ax,'on');box(ax,'off');saveFigure(fig,stem,q.pngResolution);close(fig);
end

function plotTurnCurves(T,outDir,q)
fig=figure('Visible',char(q.figureVisible),'Color','w','Position',[100 100 1100 1050]);
layout=tiledlayout(fig,3,2,'TileSpacing','compact','Padding','compact');
colors=morphColors();order=["surface","molino","pachon"];categories=["all","positive","negative"];signals=["observed","corrected"];
for c=1:3
 for z=1:2
  ax=nexttile(layout);hold(ax,'on');
  if isempty(T)||~all(ismember({'category','signal'},T.Properties.VariableNames))
   U=table();text(ax,.5,.5,'Insufficient data','Units','normalized','HorizontalAlignment','center');
  else,U=T(T.category==categories(c)&T.signal==signals(z),:);
  end
  if ~isempty(U)
   for m=1:3
    V=U(U.morph==order(m),:);if isempty(V),continue,end
    fill(ax,[V.lagSeconds;flipud(V.lagSeconds)],[V.ciLower;flipud(V.ciUpper)],colors(m,:), ...
     'FaceAlpha',.13,'EdgeColor','none','HandleVisibility','off');
    plot(ax,V.lagSeconds,V.meanRate,'Color',colors(m,:),'LineWidth',2.5,'DisplayName',char(order(m)));
   end
  end
  xline(ax,0,'--k','Turn','HandleVisibility','off');if z==2,yline(ax,0,':k','HandleVisibility','off');end
  xlabel(ax,'Lag from turn onset (s)');ylabel(ax,'Onsets/s/assembly');title(ax,categories(c)+" turns | "+signals(z));
  grid(ax,'on');box(ax,'off');if c==1&&z==1,legend(ax,'Location','best');end
 end
end
title(layout,'Turn-triggered assembly activation onsets | fish mean and bootstrap 95% CI');
saveFigure(fig,fullfile(outDir,'turn_triggered_assembly_onset_rates'),q.pngResolution);close(fig);
end

function plotScalar(T,field,yLabel,titleText,stem,q,stream,statistics)
fig=figure('Visible',char(q.figureVisible),'Color','w','Position',[100 100 650 500]);ax=axes(fig);hold(ax,'on');
colors=morphColors();order=["surface","molino","pachon"];
for m=1:3
 v=double(T.(field)(T.morph==order(m)));v=v(isfinite(v));if isempty(v),continue,end
 [lo,hi]=bootstrapCI(v,q.nFishBootstraps,stream);mu=mean(v);
 bar(ax,m,mu,.46,'FaceColor',colors(m,:),'FaceAlpha',.78,'EdgeColor','k','LineWidth',.8);
 errorbar(ax,m,mu,mu-lo,hi-mu,'k','LineStyle','none','LineWidth',1.7,'CapSize',8);
 scatter(ax,m+assembly_deterministic_jitter(numel(v),.11),v,25,'MarkerFaceColor',assembly_mix_with_white(colors(m,:),.38),'MarkerEdgeColor','k','LineWidth',.45,'MarkerFaceAlpha',.82);
end
xticks(ax,1:3);xticklabels(ax,{'Surface','Molino','Pachon'});ylabel(ax,yLabel);title(ax,titleText);
U=statistics(statistics.metric==string(field),:);addKWBrackets(ax,U);
grid(ax,'on');box(ax,'off');saveFigure(fig,stem,q.pngResolution);close(fig);
end

function addKWBrackets(ax,S)
if isempty(S),return,end
yl=ylim(ax);span=diff(yl);if ~isfinite(span)||span<=0,span=1;end
base=yl(2)+.08*span;step=.13*span;
for k=1:height(S)
 x1=1;x2=find(["surface","molino","pachon"]==S.group2(k),1);if isempty(x2),continue,end
 y=base+(k-1)*step;p=S.pValue(k);stars=pStars(p);
 plot(ax,[x1 x1 x2 x2],[y y+.025*span y+.025*span y],'k-','LineWidth',1.1,'Clipping','off');
 text(ax,mean([x1 x2]),y+.035*span,sprintf('%s  KW p=%s',stars,formatP(p)), ...
  'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8,'Clipping','off');
end
ylim(ax,[yl(1) base+(height(S)-1)*step+.16*span]);
end

function s=pStars(p)
if ~isfinite(p),s='n/a';elseif p<.001,s='***';elseif p<.01,s='**';elseif p<.05,s='*';else,s='ns';end
end

function s=formatP(p)
if ~isfinite(p),s='NA';elseif p<.001,s=sprintf('%.2g',p);else,s=sprintf('%.3f',p);end
end

function [lo,hi]=bootstrapCI(v,nBoot,stream)
v=v(isfinite(v));lo=NaN;hi=NaN;if isempty(v),return,end
if numel(v)==1,lo=v;hi=v;return,end
B=zeros(nBoot,1);for b=1:nBoot,B(b)=mean(v(randi(stream,numel(v),numel(v),1)));end
z=prctile(B,[2.5 97.5]);lo=z(1);hi=z(2);
end

function [t,source]=loadImagingTime(row,n)
assert(ismember('candidatePath',row.Properties.VariableNames)&&isfile(char(row.candidatePath)), ...
 'A valid manifest candidatePath is required for the behavior-imaging alignment.');
C=load(char(row.candidatePath),'time_s');
assert(isfield(C,'time_s')&&numel(C.time_s)>=n,'Candidate time_s missing or shorter than assembly traces.');
t=double(C.time_s(1:n));source="candidate time_s";
assert(all(isfinite(t))&&all(diff(t)>0),'Imaging timestamps must be finite and strictly increasing.');
end

function segment=continuousSegments(t)
dt=median(diff(t));breaks=[true;diff(t)>1.5*dt|diff(t)<=0];segment=cumsum(breaks);
end

function path=findBehaviorFile(row)
if ismember('behaviorPath',row.Properties.VariableNames)&&isfile(char(row.behaviorPath)),path=string(row.behaviorPath);return,end
assert(ismember('recordingPath',row.Properties.VariableNames)&&isfolder(char(row.recordingPath)), ...
 'Manifest lacks a valid recordingPath for behavior discovery.');
patterns={'*swimResults_pass2*.mat','*swimResults*.mat','*BEHAVIOR*.mat','*behavior*.mat'};D=struct([]);
for k=1:numel(patterns),x=dir(fullfile(char(row.recordingPath),patterns{k}));x=x(~[x.isdir]);D=[D;x(:)];end %#ok<AGROW>
assert(~isempty(D),'No behavior file matching the behavior-model patterns.');
[~,u]=unique(lower(string(fullfile({D.folder},{D.name}))),'stable');D=D(u);
score=10*contains(lower(string({D.name})),'swimresults_pass2')+3*contains(lower(string({D.name})),'pre_first_event');
[~,ord]=sortrows([score(:),[D.datenum]'],[-1 -2]);D=D(ord);valid=false(numel(D),1);
for k=1:numel(D),w=whos('-file',fullfile(D(k).folder,D(k).name));valid(k)=any(strcmp({w.name},'pass2FileResult'));end
D=D(valid);assert(~isempty(D),'No discovered behavior file contains pass2FileResult.');path=string(fullfile(D(1).folder,D(1).name));
end

function n=inferBehaviorLength(B)
n=NaN;for name={'behaviorTimeSec','heading_est','tail_angle'},if isfield(B,name{1})&&~isempty(B.(name{1})),n=numel(B.(name{1}));return,end,end
if isfield(B,'nFrames')&&isfinite(double(B.nFrames(1))),n=round(double(B.nFrames(1)));end
end

function [t,source]=chooseBehaviorTime(B,n,fps,tCa,mode)
local=(0:n-1)'/fps;startFrame=1;
if isfield(B,'startFrame')&&isfinite(double(B.startFrame(1))),startFrame=round(double(B.startFrame(1)));end
fromStart=((startFrame:startFrame+n-1)'-1)/fps;explicit=[];
if isfield(B,'behaviorTimeSec')&&numel(B.behaviorTimeSec)==n,explicit=double(B.behaviorTimeSec(:));end
switch mode
 case "local",t=local;source="local: (frame-1)/fps";
 case "startframe",t=fromStart;source="startFrame-adjusted";
 case "behaviortimesec",assert(~isempty(explicit),'behaviorTimeSec requested but unavailable.');t=explicit;source="pass2FileResult.behaviorTimeSec";
 otherwise
  C={local,fromStart};N=["local: (frame-1)/fps","startFrame-adjusted"];if ~isempty(explicit),C{end+1}=explicit;N(end+1)="behaviorTimeSec";end
  best=1;coverage=-Inf;distance=Inf;
  for k=1:numel(C),f=C{k}(isfinite(C{k}));if numel(f)<2,continue,end
   cv=sum(isfinite(tCa)&tCa>=min(f)&tCa<=max(f));ds=abs(min(tCa)-min(f));
   if cv>coverage||(cv==coverage&&ds<distance),best=k;coverage=cv;distance=ds;end
  end
  t=C{best};source="auto -> "+N(best);
end
assert(numel(t)==n&&all(diff(t(isfinite(t)))>0),'Invalid behavior time vector.');
end

function a=cleanAngleUnitsToRad(a)
a=double(a(:));f=a(isfinite(a));if ~isempty(f)&&max(abs(f))>2*pi+.5,a=deg2rad(a);end;a=mod(a+pi,2*pi)-pi;
end

function [startFrame,endFrame]=parseSwimFrames(S,n)
if size(S,1)==2,startFrame=S(1,:)';endFrame=S(2,:)';elseif size(S,2)==2,startFrame=S(:,1);endFrame=S(:,2);else,error('swimFrames must be 2 x nBouts or nBouts x 2.');end
k=min([n,numel(startFrame),numel(endFrame)]);startFrame=startFrame(1:k);endFrame=endFrame(1:k);
end

function T=parametersTable(q,a)
S=q;if isfield(S,'batchResults'),S=rmfield(S,'batchResults');end
S.activationOptions=a;names=fieldnames(S);value=strings(numel(names),1);
for k=1:numel(names),try,value(k)=string(jsonencode(S.(names{k})));catch,value(k)=string(S.(names{k}));end,end
T=table(string(names),value,'VariableNames',{'parameter','value'});
extra=table(["activationEpisodeDefinition";"assemblySignal";"pairSymmetry";"pairAggregation"; ...
 "turnRateNormalization";"turnAggregation";"curveStatus";"multipleTestingCorrection";"turnThresholdSource"], ...
 ["compute_assembly_activations with the recorded shared activationOptions; episode onset is the first frame of final duration-filtered active run"; ...
 "unmasked S.pca deltaF/F arithmetic mean across member neurons"; ...
 "unordered-pair curve is mean of forward curve and its time reversal"; ...
 "pairs within fish, then fish equally within morph"; ...
 "assembly onsets divided by valid turn exposure, native imaging-frame duration, and assembly count"; ...
 "assemblies within fish, then fish equally within morph"; ...
 "full lag curves and positive/negative turn analyses are exploratory"; ...
 "No multiple-testing correction; independent two-group one-way Kruskal-Wallis tests match compare_assembly_properties_across_morphs"; ...
 "0.239 rad from compare_class_features_across_morphs C.TurnThreshold; native angle multiplied by -1 so positive is CCW"], ...
 'VariableNames',{'parameter','value'});
T=[T;extra];
end

function T=shiftRecordTable(S)
T=structToTableOrEmpty(S,{'morph','fish','nAllowedOffsets','nShuffles','offsetFramesJson','offsetSecondsJson'}, ...
 {'string','string','double','double','string','string'});
end

function T=structToTableOrEmpty(S,names,types)
if ~isempty(S),T=struct2table(S);return,end
vars=cell(size(names));for k=1:numel(names),if types{k}=="string",vars{k}=strings(0,1);else,vars{k}=zeros(0,1);end,end
T=table(vars{:},'VariableNames',names);
end

function T=emptyPairTable()
T=table(strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),nan(0,1),nan(0,1), ...
 nan(0,1),zeros(0,1),false(0,1),false(0,1),nan(0,1),zeros(0,1),strings(0,1), ...
 'VariableNames',{'morph','fish','assembly1','assembly2','phase1Rad','phase2Rad', ...
 'circularPhaseDistanceDeg','nSharedMembers','phaseEligible','inPrimaryAnalysis', ...
 'zeroLagCorrelation','nValidZeroLagSamples','exclusionReason'});
end

function T=emptyPhaseDistancePairTable()
T=table(strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),nan(0,1),nan(0,1), ...
 nan(0,1),zeros(0,1),false(0,1),nan(0,1),zeros(0,1),strings(0,1), ...
 'VariableNames',{'morph','fish','assembly1','assembly2','phase1Rad','phase2Rad', ...
 'circularPhaseDistanceDeg','nSharedMembers','included','zeroLagCorrelation', ...
 'nValidZeroLagSamples','exclusionReason'});
end

function T=emptyPairCurveTable()
T=table(strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),zeros(0,1),nan(0,1),nan(0,1), ...
 'VariableNames',{'morph','fish','assembly1','assembly2','lagFrames','lagSeconds','symmetrizedCorrelation'});
end

function T=emptyTurnEventTable()
T=table(strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),nan(0,1),nan(0,1),nan(0,1), ...
 nan(0,1),false(0,1),strings(0,1),strings(0,1),strings(0,1),strings(0,1), ...
 'VariableNames',{'morph','fish','boutID','behaviorOnsetFrame','nativeAngleRad','paperAngleRad', ...
 'behaviorOnsetTimeSec','imagingIndex','included','direction','exclusionReason','angleSource','behaviorTimeSource'});
end

function C=morphColors(),C=assembly_figure_style().morphColors;end
function x=safeMean(x),x=x(isfinite(x));if isempty(x),x=NaN;else,x=mean(x);end,end
function x=safeMedian(x),x=x(isfinite(x));if isempty(x),x=NaN;else,x=median(x);end,end
function y=ternary(condition,a,b),if condition,y=a;else,y=b;end,end
function s=joinReason(a,b),if strlength(a)==0,s=b;else,s=a+"; "+b;end,end
function S=appendStruct(S,r),if isempty(S),S=r;else,S(end+1)=r;end,end
function requireFields(S,names,label),missing=names(~isfield(S,names));assert(isempty(missing),'%s missing: %s',label,strjoin(missing,', '));end
function saveFigure(fig,stem,resolution)
pngFile=[char(stem) '.png'];svgFile=[char(stem) '.svg'];
try
 style_assembly_figure(fig);exportgraphics(fig,pngFile,'Resolution',resolution);
catch ME
 if ~exportOptionUnsupported(ME),rethrow(ME);end
 print(fig,pngFile,'-dpng',sprintf('-r%d',resolution));
end
try
 style_assembly_figure(fig);exportgraphics(fig,svgFile,'ContentType','vector');
catch ME
 if ~exportOptionUnsupported(ME),rethrow(ME);end
 print(fig,svgFile,'-dsvg');
end
end

function tf=exportOptionUnsupported(ME)
message=lower(string(ME.message));
tf=contains(message,"additional, unrecognized inputs")|| ...
 contains(message,"unrecognized parameter")||contains(message,"too many input");
end
