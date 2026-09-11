function A=compute_assembly_activations(membership,activityRaster,fps,options)
%COMPUTE_ASSEMBLY_ACTIVATIONS Empirical-null assembly activation with hysteresis.
% S.activity is used as stored (binary threshold >0); events are not extended.
% Each null realization independently circularly shifts every complete neuron
% trace, preserving per-neuron rates, run durations, and autocorrelation.
if nargin<4||isempty(options),options=struct;end
D=assembly_activation_defaults();names=fieldnames(D);
for i=1:numel(names),if ~isfield(options,names{i}),options.(names{i})=D.(names{i});end,end
validateInputs(membership,activityRaster,fps,options);
membership=logical(membership);B=activityRaster>0;
[nFrames,nNeurons]=size(B);nAssemblies=size(membership,2);
assemblySize=sum(membership,1)';
K=(double(B)*double(membership))';
recruitment=K./assemblySize;
populationActive=sum(B,2)';
denominator=K+populationActive;
matchingIndex=zeros(size(K));good=denominator>0;
matchingIndex(good)=2*K(good)./denominator(good);

safeShiftFrames=max(1,ceil(options.minNullShiftSeconds*fps));
allowedShifts=safeShiftFrames:(nFrames-safeShiftFrames);
allowedShifts=allowedShifts(allowedShifts>0&allowedShifts<nFrames);
assert(~isempty(allowedShifts),['No admissible circular shift: recording must exceed twice ' ...
 'minNullShiftSeconds after conversion to frames.']);
maxAssemblySize=max([assemblySize;0]);
nullCountHistogram=zeros(nAssemblies,maxAssemblySize+1);
nullShifts=zeros(options.nNullShuffles,nNeurons,'uint32');
oldRng=rng;cleanup=onCleanup(@()rng(oldRng));rng(options.randomSeed,'twister'); %#ok<NASGU>
for s=1:options.nNullShuffles
 shifted=false(nFrames,nNeurons);
 for neuron=1:nNeurons
  shift=allowedShifts(randi(numel(allowedShifts)));
  nullShifts(s,neuron)=shift;
  shifted(:,neuron)=circshift(B(:,neuron),shift);
 end
 nullK=double(shifted)*double(membership);
 for a=1:nAssemblies
  nullCountHistogram(a,:)=nullCountHistogram(a,:)+ ...
   accumarray(nullK(:,a)+1,1,[maxAssemblySize+1 1])';
 end
end
nNullSamples=options.nNullShuffles*nFrames;
nullUpperTailProbability=zeros(size(nullCountHistogram));
nullCountThreshold=zeros(nAssemblies,1);
empiricalP=ones(nAssemblies,nFrames);
for a=1:nAssemblies
 tail=fliplr(cumsum(fliplr(nullCountHistogram(a,:))));
 tail=(tail+1)/(nNullSamples+1);
 nullUpperTailProbability(a,:)=tail;
 threshold=find(tail<=options.activationNullAlpha,1,'first');
 if isempty(threshold),nullCountThreshold(a)=assemblySize(a)+1;
 else,nullCountThreshold(a)=threshold-1;end
 empiricalP(a,:)=tail(K(a,:)+1);
end

startEligible=K>=options.minActiveMembers & ...
 recruitment>=options.onRecruitmentFraction & ...
 matchingIndex>=options.minMatchingIndex & ...
 empiricalP<=options.activationNullAlpha;
active=false(nAssemblies,nFrames);
for a=1:nAssemblies
 isActive=false;
 for t=1:nFrames
  if ~isActive
   if startEligible(a,t),isActive=true;active(a,t)=true;end
  elseif recruitment(a,t)<options.offRecruitmentFraction || ...
    matchingIndex(a,t)<options.minMatchingIndex
   isActive=false;
  else
   active(a,t)=true;
  end
 end
end
minDurationFrames=max(1,round(options.minDurationSeconds*fps));
for a=1:nAssemblies
 episodeBounds=logicalRuns(active(a,:));
 for e=1:size(episodeBounds,1)
  if episodeBounds(e,2)-episodeBounds(e,1)+1<minDurationFrames
   active(a,episodeBounds(e,1):episodeBounds(e,2))=false;
  end
 end
end
onsetIndices=cell(nAssemblies,1);offsetIndices=cell(nAssemblies,1);
for a=1:nAssemblies
 episodeBounds=logicalRuns(active(a,:));
 onsetIndices{a}=episodeBounds(:,1);offsetIndices{a}=episodeBounds(:,2);
end
A=struct('K',K,'recruitment',recruitment,'populationActive',populationActive, ...
 'matchingIndex',matchingIndex,'empiricalP',empiricalP, ...
 'nullCountThreshold',nullCountThreshold,'nullCountHistogram',nullCountHistogram, ...
 'nullUpperTailProbability',nullUpperTailProbability,'active',active, ...
 'onsetIndices',{onsetIndices},'offsetIndices',{offsetIndices}, ...
 'startEligible',startEligible,'assemblySize',assemblySize, ...
 'minDurationFrames',minDurationFrames,'safeShiftFrames',safeShiftFrames, ...
 'allowedShifts',allowedShifts,'nullShifts',nullShifts,'options',options, ...
 'activityDefinition','S.activity > 0 without temporal extension');
end

function validateInputs(membership,activityRaster,fps,o)
assert(ismatrix(membership)&&ismatrix(activityRaster),'Inputs must be matrices.');
assert(size(activityRaster,2)==size(membership,1),'Raster columns must match membership rows.');
assert(all(sum(logical(membership),1)>0),'Assemblies must contain at least one member.');
assert(isscalar(fps)&&isfinite(fps)&&fps>0,'fps must be positive and finite.');
assert(isscalar(o.minActiveMembers)&&o.minActiveMembers>=1&&o.minActiveMembers==round(o.minActiveMembers));
assert(isscalar(o.onRecruitmentFraction)&&o.onRecruitmentFraction>=0&&o.onRecruitmentFraction<=1);
assert(isscalar(o.offRecruitmentFraction)&&o.offRecruitmentFraction>=0&&o.offRecruitmentFraction<=o.onRecruitmentFraction);
assert(isscalar(o.minMatchingIndex)&&o.minMatchingIndex>=0&&o.minMatchingIndex<=1);
assert(isscalar(o.activationNullAlpha)&&o.activationNullAlpha>0&&o.activationNullAlpha<=1);
assert(isscalar(o.nNullShuffles)&&o.nNullShuffles>=1&&o.nNullShuffles==round(o.nNullShuffles));
assert(isscalar(o.minDurationSeconds)&&isfinite(o.minDurationSeconds)&&o.minDurationSeconds>=0);
assert(isscalar(o.minNullShiftSeconds)&&isfinite(o.minNullShiftSeconds)&&o.minNullShiftSeconds>=0);
assert(isscalar(o.randomSeed)&&isfinite(o.randomSeed)&&o.randomSeed==round(o.randomSeed));
end

function bounds=logicalRuns(x)
d=diff([false,logical(x),false]);bounds=[find(d==1)',find(d==-1)'-1];
end
