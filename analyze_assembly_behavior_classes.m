function result=analyze_assembly_behavior_classes(member,active,preferredPhaseRad,classLabels,morph,fish)
%ANALYZE_ASSEMBLY_BEHAVIOR_CLASSES Map aligned AHV cell classes to assemblies.
% Inputs are for one fish and must share candidate-ID order.
member=logical(member);active=logical(active);
preferredPhaseRad=double(preferredPhaseRad(:));classLabels=string(classLabels(:));
morph=string(morph);fish=string(fish);nCells=size(member,1);nAssemblies=size(member,2);
assert(numel(preferredPhaseRad)==nCells&&numel(classLabels)==nCells,'Inputs must be aligned.');
assert(size(active,1)==nAssemblies,'active must have one row per assembly.');
canonical=repmat("",nCells,1);
canonical(strcmpi(classLabels,"CW"))="CW";canonical(strcmpi(classLabels,"CCW"))="CCW";
canonical(strcmpi(classLabels,"Symmetric"))="Symmetric";
canonical(strcmpi(classLabels,"NotFit")|ismissing(classLabels)|strlength(classLabels)==0)="NotFit";
unexpected=canonical=="";
assert(~any(unexpected),'Unexpected class label(s): %s',strjoin(unique(classLabels(unexpected)),', '));
classOrder=["CW","CCW","Symmetric"];
assemblyID=(1:nAssemblies)';assemblySize=sum(member,1)';nValidClassMembers=zeros(nAssemblies,1);
dominantClass=repmat("NotFit",nAssemblies,1);purityScore=nan(nAssemblies,1);assemblyMeanPhase=nan(nAssemblies,1);
for a=1:nAssemblies
 ix=member(:,a);labels=canonical(ix);labels=labels(labels~="NotFit");nValidClassMembers(a)=numel(labels);
 if ~isempty(labels)
  counts=arrayfun(@(c)sum(labels==c),classOrder);[largest,k]=max(counts);
  dominantClass(a)=classOrder(k);purityScore(a)=largest/numel(labels);
 end
 assemblyMeanPhase(a)=circularMeanStandard(preferredPhaseRad(ix&isfinite(preferredPhaseRad)));
end
assemblyTable=table(repmat(morph,nAssemblies,1),repmat(fish,nAssemblies,1),assemblyID,assemblySize, ...
 nValidClassMembers,dominantClass,purityScore,assemblyMeanPhase,'VariableNames',{'morph','fish', ...
 'assemblyID','assemblySize','nValidClassMembers','dominantClass','purityScore','assemblyMeanPhase'});
if nAssemblies>=2
 C=corr(double(active'),'Rows','pairwise');[assemblyA,assemblyB]=find(triu(true(nAssemblies),1));
 temporalCorrelation=C(sub2ind([nAssemblies,nAssemblies],assemblyA,assemblyB));
 classA=dominantClass(assemblyA);classB=dominantClass(assemblyB);relationship=repmat("Other",numel(assemblyA),1);
 relationship((classA=="CW"&classB=="CCW")|(classA=="CCW"&classB=="CW"))="Opponent";
 relationship((classA=="CW"&classB=="CW")|(classA=="CCW"&classB=="CCW"))="Same";
 pairTable=table(repmat(morph,numel(assemblyA),1),repmat(fish,numel(assemblyA),1),assemblyA,assemblyB, ...
  classA,classB,relationship,temporalCorrelation,'VariableNames',{'morph','fish','assemblyA','assemblyB', ...
  'classA','classB','relationship','temporalCorrelation'});
else
 pairTable=table(strings(0,1),strings(0,1),zeros(0,1),zeros(0,1),strings(0,1),strings(0,1), ...
  strings(0,1),zeros(0,1),'VariableNames',{'morph','fish','assemblyA','assemblyB','classA','classB', ...
  'relationship','temporalCorrelation'});
end
poolMeanPhase=nan(3,1);poolN=zeros(3,1);
for k=1:3
 values=assemblyMeanPhase(dominantClass==classOrder(k)&isfinite(assemblyMeanPhase));
 poolN(k)=numel(values);poolMeanPhase(k)=circularMeanStandard(values);
end
classPoolTable=table(repmat(morph,3,1),repmat(fish,3,1),classOrder',poolN,poolMeanPhase, ...
 'VariableNames',{'morph','fish','dominantClass','nAssemblies','circularMeanPhase'});
cwPhase=poolMeanPhase(1);ccwPhase=poolMeanPhase(2);
if isfinite(cwPhase)&&isfinite(ccwPhase)
 cwCcwPhaseSeparation=abs(circularDifference(cwPhase,ccwPhase));cwCcwAntiphaseError=abs(pi-cwCcwPhaseSeparation);
else
 cwCcwPhaseSeparation=NaN;cwCcwAntiphaseError=NaN;
end
validPurity=purityScore(isfinite(purityScore));
if isempty(validPurity),medianAssemblyPurity=NaN;else,medianAssemblyPurity=median(validPurity);end
fishTable=table(morph,fish,nAssemblies,medianAssemblyPurity,cwPhase,ccwPhase,cwCcwPhaseSeparation, ...
 cwCcwAntiphaseError,'VariableNames',{'morph','fish','nAssemblies','medianAssemblyPurity', ...
 'cwPoolMeanPhase','ccwPoolMeanPhase','cwCcwPhaseSeparation','cwCcwAntiphaseError'});
result=struct('assemblyTable',assemblyTable,'pairTable',pairTable,'fishTable',fishTable, ...
 'classPoolTable',classPoolTable);
end
function value=circularMeanStandard(phases)
phases=phases(isfinite(phases));if isempty(phases),value=NaN;return,end
value=atan2(mean(sin(phases)),mean(cos(phases)));if value<=-pi,value=pi;end
end
function value=circularDifference(a,b)
value=atan2(sin(a-b),cos(a-b));
end
