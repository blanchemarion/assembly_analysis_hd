function report=compare_assembly_spatial_compactness_across_morphs(batchResults,varargin)
%COMPARE_ASSEMBLY_SPATIAL_COMPACTNESS_ACROSS_MORPHS Hansen-style analysis.
% One detected assembly is one observation; assemblies remain nested in fish.
% Membership and activation are never redefined. Distances are unweighted.
% Hansen's H<1.5 cutoff was selected by visual inspection of a bimodal
% distribution and is bin-dependent; fixed edges are therefore saved/audited.
if nargin<1||isempty(batchResults)
 root=fileparts(mfilename('fullpath'));
 saved=load(fullfile(root,'data_processed','assemblies','session_manifest.mat'),'manifest','opts');
 batchResults=struct('manifest',saved.manifest,'options',saved.opts);
end
p=inputParser;addRequired(p,'batchResults',@(x)isstruct(x)&&isfield(x,'manifest')&&isfield(x,'options'));
addParameter(p,'outputDir','',@(x)ischar(x)||isstring(x));
addParameter(p,'distanceBinWidthUm',10,@positiveScalar);
addParameter(p,'maxDistanceUm',500,@positiveScalar);
addParameter(p,'sensitivityBinWidthsUm',[5 10 20],@(x)isnumeric(x)&&isvector(x)&&all(isfinite(x)&x>0));
addParameter(p,'compactEntropyCutoff',1.5,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
addParameter(p,'nBootstrap',2000,@positiveInteger);
addParameter(p,'nClusterBootstrap',2000,@positiveInteger);
addParameter(p,'randomSeed',1,@(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x==round(x));
addParameter(p,'figureVisible','off',@(x)ischar(x)||isstring(x));
addParameter(p,'pngResolution',300,@(x)isnumeric(x)&&isscalar(x)&&x>=150);
parse(p,batchResults,varargin{:});q=p.Results;
q.sensitivityBinWidthsUm=unique(double(q.sensitivityBinWidthsUm(:)'),'stable');
if ~ismember(q.distanceBinWidthUm,q.sensitivityBinWidthsUm)
 q.sensitivityBinWidthsUm=sort([q.sensitivityBinWidthsUm q.distanceBinWidthUm]);
end
edges=fixedEdges(q.distanceBinWidthUm,q.maxDistanceUm);
binSpec=binSpecification(edges,'um');
if strlength(string(q.outputDir))==0
 outputName=sprintf('assembly_spatial_compactness_%s_bin%gum_max%gum', ...
  char(batchResults.options.analysisLabel),q.distanceBinWidthUm,q.maxDistanceUm);
 outputDir=fullfile(batchResults.options.figureRoot,strrep(outputName,'.','p'));
else,outputDir=char(q.outputDir);end
if ~isfolder(outputDir),mkdir(outputDir);end
oldRng=rng;cleanup=onCleanup(@()rng(oldRng));rng(q.randomSeed,'twister'); %#ok<NASGU>
M=batchResults.manifest;M=M(M.status=="processed"|M.status=="existing",:);
assert(all(ismember({'morph','recording','recordingPath','outputPath'},M.Properties.VariableNames)), ...
 'Manifest lacks required recording/output paths.');
rows=struct([]);exclusions=struct([]);distanceStore={};mappingAudit=struct([]);
fprintf('\n=== Hansen-style assembly spatial compactness: %d recordings ===\n',height(M));
for i=1:height(M)
 morph=lower(string(M.morph(i)));fish=string(M.recording(i));resultFile=char(M.outputPath(i));
 fprintf('[%d/%d] %s/%s\n',i,height(M),morph,fish);
 try
  assert(isfile(resultFile),'Assembly result is missing: %s',resultFile);
  X=load(resultFile,'model','S');assert(isfield(X,'model')&&isfield(X,'S'),'Saved model/S missing.');
  requireFields(X.model,{'membership','candidateRoiIdx','candidateCellIDs','nTotalRois'},'model');
  requireFields(X.S,{'candidateRoiIdx','candidateCellIDs'},'S');
  membership=logical(X.model.membership);roi=double(X.model.candidateRoiIdx(:));
  cellID=double(X.model.candidateCellIDs(:));
  assert(size(membership,1)==numel(roi)&&numel(roi)==numel(cellID), ...
   'Membership rows, candidateRoiIdx and candidateCellIDs are not aligned.');
  assert(isequal(roi,double(X.S.candidateRoiIdx(:)))&& ...
   isequal(cellID,double(X.S.candidateCellIDs(:))), ...
   'Model and session candidate identifiers disagree.');
  assert(numel(unique(roi))==numel(roi)&&numel(unique(cellID))==numel(cellID), ...
   'Candidate original ROI identifiers or cell IDs are not unique.');
  if size(membership,2)==0
   exclusions=appendStruct(exclusions,exclusionRow(morph,fish,NaN,NaN,NaN, ...
    "zero assemblies detected","",resultFile));
   fprintf('  Excluded from comparison: zero assemblies detected.\n');
   continue
  end
  allCellsFile=findAllCellsFile(char(M.recordingPath(i)));
  [coordinates,nAll,pixelLengthX,pixelLengthY,coordinateSource,distanceCheckError]=loadAllCoordinates(allCellsFile,roi);
  assert(nAll==double(X.model.nTotalRois),'ALL_CELLS ROI count differs from model.nTotalRois.');
  audit=struct('morph',morph,'fish',fish,'nMembershipRows',size(membership,1), ...
   'nAllRois',nAll,'minOriginalRoi',min(roi),'maxOriginalRoi',max(roi), ...
   'uniqueOriginalRoiCount',numel(unique(roi)),'pixelLengthXUm',pixelLengthX, ...
   'pixelLengthYUm',pixelLengthY,'coordinateSource',string(coordinateSource), ...
   'maxStoredDistanceErrorUm',distanceCheckError,'allCellsFile',string(allCellsFile), ...
   'assemblyResultFile',string(resultFile));
  mappingAudit=appendStruct(mappingAudit,audit);
  for a=1:size(membership,2)
   evaluation=evaluate_hansen_assembly(membership(:,a),roi,coordinates,edges,q.compactEntropyCutoff);
   memberRoi=evaluation.memberOriginalRoi;assemblySize=numel(memberRoi);nValid=sum(evaluation.validCoordinate);
   if ~evaluation.evaluated
    exclusions=appendStruct(exclusions,exclusionRow(morph,fish,a,assemblySize,nValid,evaluation.reason,allCellsFile,resultFile));
    continue
   end
   metrics=evaluation.metrics;
   row=struct('morph',morph,'fish',fish,'recording',fish,'assemblyID',a, ...
    'assemblySize',assemblySize,'nValidMemberCoordinates',nValid, ...
    'nPairwiseDistances',numel(metrics.distances), ...
    'meanPairwiseDistanceUm',metrics.meanPairwiseDistance, ...
    'medianPairwiseDistanceUm',metrics.medianPairwiseDistance, ...
    'distanceDistributionEntropy',metrics.distanceEntropy, ...
    'isCompactHansen',metrics.isCompactHansen,'coordinateUnits',"um", ...
    'histogramBinSpecification',binSpec, ...
    'histogramCounts',join(string(metrics.histogramCounts),'|'), ...
    'histogramProbability',join(string(metrics.histogramProbability),'|'), ...
    'memberOriginalRoiIndices',join(string(memberRoi'),'|'), ...
    'allCellsFile',string(allCellsFile),'assemblyResultFile',string(resultFile));
   rows=appendStruct(rows,row);distanceStore{end+1,1}=metrics.distances; %#ok<AGROW>
  end
 catch ME
  warning('assembly:HansenSpatialFailed','%s/%s: %s',morph,fish,ME.message);
  exclusions=appendStruct(exclusions,exclusionRow(morph,fish,NaN,NaN,NaN, ...
   "recording failure: "+string(ME.message),"",resultFile));
 end
end
assert(~isempty(rows),'No assembly had at least two valid mapped coordinates.');
assemblyTable=struct2table(rows);exclusionTable=structTable(exclusions,exclusionTemplate());
mappingAuditTable=struct2table(mappingAudit);
[descriptiveTable,compactTable]=descriptiveStatistics(assemblyTable,q.nBootstrap,q.compactEntropyCutoff);
[kwTable,dunnTable]=rankTests(assemblyTable);
pairwiseKWTable=pairwiseKWTests(assemblyTable);
clusterBootstrapTable=fishClusterBootstrap(assemblyTable,q.nClusterBootstrap);
[sensitivityTable,sensitivityClassTable]=binSensitivity(assemblyTable,distanceStore,q.sensitivityBinWidthsUm,q.maxDistanceUm,q.compactEntropyCutoff,q.distanceBinWidthUm);
makeScatter(assemblyTable,q.compactEntropyCutoff,outputDir,q);
makeDistributions(assemblyTable,pairwiseKWTable,outputDir,q);
makeCompactPlot(compactTable,outputDir,q);
makeSensitivityPlot(sensitivityTable,outputDir,q);
writetable(assemblyTable,fullfile(outputDir,'hansen_assembly_level_metrics.csv'));
writetable(exclusionTable,fullfile(outputDir,'hansen_assembly_exclusions.csv'));
writetable(mappingAuditTable,fullfile(outputDir,'hansen_roi_mapping_audit.csv'));
writetable(descriptiveTable,fullfile(outputDir,'hansen_morph_descriptive_statistics.csv'));
writetable(compactTable,fullfile(outputDir,'hansen_compact_classification_summary.csv'));
writetable(kwTable,fullfile(outputDir,'hansen_kruskal_wallis_tests.csv'));
writetable(pairwiseKWTable,fullfile(outputDir,'hansen_pairwise_kruskal_wallis_tests.csv'));
writetable(dunnTable,fullfile(outputDir,'hansen_planned_dunn_sidak_tests.csv'));
writetable(clusterBootstrapTable,fullfile(outputDir,'hansen_fish_cluster_bootstrap.csv'));
writetable(sensitivityTable,fullfile(outputDir,'hansen_bin_width_sensitivity.csv'));
writetable(sensitivityClassTable,fullfile(outputDir,'hansen_assembly_bin_sensitivity.csv'));
settings=q;settings.distanceBinEdges=edges;settings.histogramBinSpecification=binSpec;
settings.coordinateUnits="um";settings.coordinateMapping="membership row -> model.candidateRoiIdx -> full ALL_CELLS ROI";
settings.coordinateSource="centroid of ALL_CELLS.cells ROI pixels, scaled by pixelLengthX/Y";
settings.observationalUnit="assembly (assemblies are nested within fish)";
settings.methodQualification="Hansen-style, not exact numerical reproduction because published bin edges are not fully specified";
settings.compactCutoffQualification="H < 1.5; Hansen selected cutoff by visual inspection of a bimodal distribution; classification depends on binning";
settingsTable=settingsToTable(settings);
writetable(settingsTable,fullfile(outputDir,'hansen_analysis_settings.csv'));
report=struct('assemblyLevel',assemblyTable,'exclusions',exclusionTable,'mappingAudit',mappingAuditTable, ...
 'descriptiveStatistics',descriptiveTable,'compactSummary',compactTable, ...
 'kruskalWallis',kwTable,'pairwiseKruskalWallis',pairwiseKWTable,'plannedDunnSidak',dunnTable, ...
 'fishClusterBootstrap',clusterBootstrapTable,'binWidthSensitivity',sensitivityTable, ...
 'assemblyBinSensitivity',sensitivityClassTable,'pairwiseDistances',{distanceStore}, ...
 'settings',settings,'settingsTable',settingsTable,'outputDir',string(outputDir));
save(fullfile(outputDir,'hansen_spatial_compactness_report.mat'),'report','-v7.3');
fprintf('Hansen-style spatial comparison complete: %s\n',outputDir);
end

function tf=positiveScalar(x),tf=isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0;end
function tf=positiveInteger(x),tf=positiveScalar(x)&&x==round(x);end
function requireFields(S,names,label),missing=names(~isfield(S,names));assert(isempty(missing),'%s missing: %s',label,strjoin(missing,', '));end
function S=appendStruct(S,x),if isempty(S),S=x;else,S(end+1)=x;end,end %#ok<AGROW>
function edges=fixedEdges(width,maximum),edges=[0:width:maximum Inf];if edges(end-1)<maximum,edges=[edges(1:end-1) maximum Inf];end,end
function s=binSpecification(edges,units),finite=edges(isfinite(edges));s=sprintf('fixed edges %g:%g:%g %s plus +Inf overflow',finite(1),finite(2)-finite(1),finite(end),units);end

function file=findAllCellsFile(folder)
f=dir(fullfile(folder,'*ALL_CELLS.mat'));names=upper(string({f.name}));
f=f(~contains(names,"SELECTED")&~contains(names,"NON_SELECTED")&~contains(names,"ANATOMICAL"));
assert(numel(f)==1,'Expected exactly one canonical ALL_CELLS file in %s; found %d.',folder,numel(f));file=fullfile(f.folder,f.name);
end

function [coordinates,nAll,px,py,source,maxError]=loadAllCoordinates(file,candidateRoiIdx)
A=load(file,'cells','cell_per','avg','pixelLengthX','pixelLengthY','distances','cell_number');
assert(isfield(A,'pixelLengthX')&&isfield(A,'pixelLengthY')&&positiveScalar(A.pixelLengthX)&&positiveScalar(A.pixelLengthY), ...
 'Reliable pixel calibration is unavailable in %s.',file);px=double(A.pixelLengthX);py=double(A.pixelLengthY);
assert(isfield(A,'avg')&&ismatrix(A.avg),'ALL_CELLS avg image is required to decode ROI masks.');
if isfield(A,'cells'),nAll=numel(A.cells);coordinates=nan(nAll,2);
 for r=1:nAll
  idx=double(A.cells{r}(:));idx=idx(isfinite(idx)&idx>=1&idx<=numel(A.avg)&idx==round(idx));
  if ~isempty(idx),[y,x]=ind2sub(size(A.avg),idx);coordinates(r,:)=[mean(x)*px mean(y)*py];end
 end
 source='centroid of full-ROI ALL_CELLS.cells masks scaled by pixelLengthX/Y';
elseif isfield(A,'cell_per'),nAll=numel(A.cell_per);coordinates=nan(nAll,2);
 for r=1:nAll
  xy=double(A.cell_per{r});if size(xy,2)>=2&&any(all(isfinite(xy(:,1:2)),2)),coordinates(r,:)=mean(xy(:,1:2),1,'omitmissing').*[px py];end
 end
 source='centroid of full-ROI ALL_CELLS.cell_per boundaries scaled by pixelLengthX/Y';
else,error('ALL_CELLS has neither cells nor cell_per ROI geometry.');end
if isfield(A,'cell_number'),assert(double(A.cell_number)==nAll,'cell_number disagrees with ROI geometry count.');end
assert(numel(unique(candidateRoiIdx))==numel(candidateRoiIdx)&&all(candidateRoiIdx>=1&candidateRoiIdx<=nAll&candidateRoiIdx==round(candidateRoiIdx)), ...
 'Saved original ROI identifiers are nonunique or outside ALL_CELLS range.');
maxError=NaN;
if isfield(A,'distances')
 assert(isequal(size(A.distances),[nAll nAll]),'ALL_CELLS distances is not full-ROI ordered.');
 ids=candidateRoiIdx(all(isfinite(coordinates(candidateRoiIdx,:)),2));ids=ids(1:min(50,numel(ids)));
 if numel(ids)>=2
  C=coordinates(ids,:);computed=sqrt((C(:,1)-C(:,1)').^2+(C(:,2)-C(:,2)').^2);
  maxError=max(abs(computed-double(A.distances(ids,ids))),[],'all');
  assert(maxError<1e-6,'Mask-derived calibrated centroids disagree with stored full-ROI distances.');
 end
end
end

function x=exclusionRow(morph,fish,id,n,nValid,reason,allFile,resultFile)
x=struct('morph',string(morph),'fish',string(fish),'recording',string(fish),'assemblyID',id, ...
 'assemblySize',n,'nValidMemberCoordinates',nValid,'reason',string(reason), ...
 'allCellsFile',string(allFile),'assemblyResultFile',string(resultFile));
end
function x=exclusionTemplate(),x=exclusionRow("","",NaN,NaN,NaN,"","","");end
function T=structTable(S,template),if isempty(S),T=struct2table(template);T(1,:)=[];else,T=struct2table(S);end,end

function [D,C]=descriptiveStatistics(T,nBoot,cutoff)
morphs=["surface","molino","pachon"];metrics=["meanPairwiseDistanceUm","distanceDistributionEntropy"];
rows=cell(6,1);r=0;compactRows=cell(3,1);
for m=1:3
 pick=T.morph==morphs(m);fish=unique(T.fish(pick));n=sum(pick);
 compact=sum(T.isCompactHansen(pick));compactRows{m}=table(morphs(m),n,numel(fish),compact,100*compact/max(n,1),cutoff, ...
  'VariableNames',{'morph','nAssemblies','nFish','nCompact','percentCompact','entropyCutoff'});
 for k=1:2
  r=r+1;x=double(T.(metrics(k))(pick));boots=nan(nBoot,1);
  for b=1:nBoot,boots(b)=mean(x(randi(numel(x),numel(x),1)));end
  ci=prctile(boots,[2.5 97.5]);rows{r}=table(morphs(m),metrics(k),n,numel(fish),mean(x),median(x),std(x),std(x)/sqrt(n),ci(1),ci(2), ...
   'VariableNames',{'morph','metric','nAssemblies','nFish','mean','median','SD','SEM','bootstrapCI95Low','bootstrapCI95High'});
 end
end
D=vertcat(rows{:});C=vertcat(compactRows{:});
end

function [K,D]=rankTests(T)
metrics=["meanPairwiseDistanceUm","distanceDistributionEntropy"];groups=["surface","molino","pachon"];
kr=cell(2,1);dr=cell(4,1);zrow=0;
for k=1:2
 y=double(T.(metrics(k)));g=string(T.morph);[p,tbl]=kruskalwallis(y,cellstr(g),'off');kr{k}=table(metrics(k),tbl{2,5},tbl{2,3},p,'VariableNames',{'metric','chiSquare','df','pValue'});
 [ranks,tieFactor]=tiedRanksWithFactor(y);N=numel(y);means=nan(3,1);ns=zeros(3,1);
 for j=1:3,sel=g==groups(j);means(j)=mean(ranks(sel));ns(j)=sum(sel);end
 for j=2:3
  zrow=zrow+1;se=sqrt(N*(N+1)/12*tieFactor*(1/ns(1)+1/ns(j)));z=(means(1)-means(j))/se;pRaw=erfc(abs(z)/sqrt(2));pSidak=1-(1-pRaw)^2;
  dr{zrow}=table(metrics(k),groups(1),groups(j),ns(1),ns(j),means(1)-means(j),z,pRaw,min(pSidak,1), ...
   'VariableNames',{'metric','group1','group2','n1','n2','meanRankDifference','z','pRaw','pSidak'});
 end
end
K=vertcat(kr{:});D=vertcat(dr{:});
end
function [r,tieFactor]=tiedRanksWithFactor(y)
r=tiedrank(y);[~,~,id]=unique(y);c=accumarray(id,1);N=numel(y);tieFactor=1-sum(c.^3-c)/(N^3-N);end

function T=fishClusterBootstrap(A,nBoot)
metrics=["meanPairwiseDistanceUm","distanceDistributionEntropy"];contrasts=["molino","pachon"];rows=cell(4,1);r=0;
for k=1:2
 for c=1:2
  r=r+1;x=clusterContrast(A,metrics(k),contrasts(c),nBoot);rows{r}=x;
 end
end
T=vertcat(rows{:});
end
function row=clusterContrast(A,metric,other,nBoot)
base="surface";fishA=unique(A.fish(A.morph==base));fishB=unique(A.fish(A.morph==other));boot=nan(nBoot,1);
for b=1:nBoot
 a=sampleFishValues(A,metric,base,fishA);d=sampleFishValues(A,metric,other,fishB);boot(b)=mean(a)-mean(d);
end
ci=prctile(boot,[2.5 97.5]);p=2*min(mean(boot<=0),mean(boot>=0));
row=table(metric,base,other,numel(fishA),numel(fishB),mean(double(A.(metric)(A.morph==base)))-mean(double(A.(metric)(A.morph==other))),ci(1),ci(2),min(p,1), ...
 'VariableNames',{'metric','group1','group2','nFish1','nFish2','observedPooledMeanDifference','clusterBootstrapCI95Low','clusterBootstrapCI95High','clusterBootstrapTwoSidedP'});
end
function values=sampleFishValues(A,metric,morph,fish)
pick=fish(randi(numel(fish),numel(fish),1));values=[];
for i=1:numel(pick),values=[values;double(A.(metric)(A.morph==morph&A.fish==pick(i)))];end %#ok<AGROW>
end

function [S,C]=binSensitivity(A,distanceStore,widths,maxDistance,cutoff,primaryWidth)
n=height(A);entropy=nan(n,numel(widths));compact=false(n,numel(widths));
for w=1:numel(widths)
 edges=fixedEdges(widths(w),maxDistance);
 for i=1:n,counts=histcounts(distanceStore{i},edges);P=counts/sum(counts);entropy(i,w)=-sum(P(P>0).*log(P(P>0)));compact(i,w)=entropy(i,w)<cutoff;
 end
end
morphs=["surface","molino","pachon"];rows=cell(numel(widths)*3,1);r=0;
for w=1:numel(widths),for m=1:3,r=r+1;pick=A.morph==morphs(m);primary=compact(:,widths==primaryWidth);if isempty(primary),primary=compact(:,w);end
 rows{r}=table(morphs(m),widths(w),sum(pick),sum(compact(pick,w)),100*mean(compact(pick,w)),sum(compact(pick,w)~=primary(pick)), ...
  'VariableNames',{'morph','binWidthUm','nAssemblies','nCompact','percentCompact','nClassificationsChangedFromPrimary'});end,end
S=vertcat(rows{:});C=table(A.morph,A.fish,A.assemblyID);
for w=1:numel(widths),C.(sprintf('entropy_bin%gum',widths(w)))=entropy(:,w);C.(sprintf('compact_bin%gum',widths(w)))=compact(:,w);end
end

function T=settingsToTable(S)
names=fieldnames(S);values=strings(numel(names),1);
for i=1:numel(names)
 v=S.(names{i});
 if isnumeric(v)||islogical(v),values(i)=join(string(v(:)'),',');
 elseif ischar(v)||isstring(v),values(i)=join(string(v(:)'),',');
 else,values(i)=string(jsonencode(v));end
end
T=table(string(names),values,'VariableNames',{'setting','value'});
end

function makeScatter(T,cutoff,out,q)
[c,n]=plotInfo();f=figure('Visible',char(q.figureVisible),'Color','w');ax=axes(f);hold(ax,'on');
for m=1:3,p=T.morph==lower(n(m));scatter(ax,T.meanPairwiseDistanceUm(p),T.distanceDistributionEntropy(p),25,'MarkerFaceColor',assembly_mix_with_white(c(m,:),.38),'MarkerEdgeColor','k','LineWidth',.45,'MarkerFaceAlpha',.82,'DisplayName',sprintf('%s (n=%d)',n(m),sum(p)));end
yline(ax,cutoff,'--k','H = 1.5','HandleVisibility','off');xlabel(ax,'Mean pairwise distance (\mum)');ylabel(ax,'Distance-distribution entropy H');legend(ax,'Location','best');grid(ax,'on');title(ax,'Hansen-style spatial compactness');saveFigure(f,out,'hansen_entropy_vs_mean_distance',q.pngResolution);
end
function makeDistributions(T,K,out,q)
[c,n]=plotInfo();metrics={'meanPairwiseDistanceUm','distanceDistributionEntropy'};labels={'Mean pairwise distance (\mum)','Distance-distribution entropy H'};
f=figure('Visible',char(q.figureVisible),'Color','w','Position',[100 100 1050 450]);
for k=1:2
 ax=subplot(1,2,k);hold(ax,'on');
 for m=1:3
  p=T.morph==lower(n(m));x=double(T.(metrics{k})(p));x=x(isfinite(x));
  if isempty(x),continue,end
  j=assembly_deterministic_jitter(numel(x),.11);scatter(ax,m+j,x,25,'MarkerFaceColor',assembly_mix_with_white(c(m,:),.38),'MarkerEdgeColor','k','LineWidth',.45,'MarkerFaceAlpha',.82);
  errorbar(ax,m,mean(x),std(x)/sqrt(numel(x)),'kd','MarkerFaceColor',c(m,:),'LineWidth',1.7,'CapSize',8);
 end
 limits=ylim(ax);lo=limits(1);hi=limits(2);span=hi-lo;
 for m=2:3
  pValue=K.pValue(K.metric==string(metrics{k})&K.group2==lower(n(m)));
  bracketY=hi+span*(.10+.12*(m-2));tick=.035*span;
  plot(ax,[1 1 m m],[bracketY-tick bracketY bracketY bracketY-tick],'-k','LineWidth',1);
  text(ax,(1+m)/2,bracketY+.01*span,significanceLabel(pValue),'HorizontalAlignment','center', ...
   'VerticalAlignment','bottom','FontWeight','bold');
 end
 ylim(ax,[lo hi+.30*span]);xlim(ax,[.5 3.5]);xticks(ax,1:3);xticklabels(ax,n);ylabel(ax,labels{k});grid(ax,'on');box(ax,'off');
 title(ax,'Two-group KW tests (uncorrected)');
end
saveFigure(f,out,'hansen_metric_distributions',q.pngResolution);
end
function K=pairwiseKWTests(T)
% Planned two-group tests on individual assemblies, matching the properties plots.
metrics=["meanPairwiseDistanceUm","distanceDistributionEntropy"];groups=["surface","molino","pachon"];
rows=cell(4,1);r=0;g=lower(string(T.morph));
for k=1:numel(metrics)
 y=double(T.(metrics(k)));
 for m=2:3
  a=y(g==groups(1)&isfinite(y));b=y(g==groups(m)&isfinite(y));
  p=NaN;chiSquare=NaN;df=NaN;
  if ~isempty(a)&&~isempty(b)
   if all([a;b]==a(1))
    p=1;chiSquare=0;df=1;
   else
    labels=[repmat(groups(1),numel(a),1);repmat(groups(m),numel(b),1)];
    [p,tbl]=kruskalwallis([a;b],cellstr(labels),'off');chiSquare=tbl{2,5};df=tbl{2,3};
   end
  end
  r=r+1;rows{r}=table(metrics(k),groups(1),groups(m),numel(a),numel(b),chiSquare,df,p, ...
   'VariableNames',{'metric','group1','group2','n1Assemblies','n2Assemblies','chiSquare','df','pValue'});
 end
end
K=vertcat(rows{:});
end
function label=significanceLabel(p)
if ~isfinite(p),label='n/a';
elseif p<=0.001,label='***';
elseif p<=0.01,label='**';
elseif p<=0.05,label='*';
else,label='ns';
end
end
function makeCompactPlot(T,out,q)
[c,n]=plotInfo();f=figure('Visible',char(q.figureVisible),'Color','w');ax=axes(f);bar(ax,1:3,T.percentCompact,.46,'FaceColor','flat','FaceAlpha',.78,'EdgeColor','k','LineWidth',.8);ax.Children.CData=c;xticks(ax,1:3);xticklabels(ax,n);ylabel(ax,'Compact assemblies (%)');title(ax,'Hansen-style compact classification (H < 1.5)');grid(ax,'on');saveFigure(f,out,'hansen_percent_compact_by_morph',q.pngResolution);
end
function makeSensitivityPlot(T,out,q)
[c,n]=plotInfo();f=figure('Visible',char(q.figureVisible),'Color','w');ax=axes(f);hold(ax,'on');for m=1:3,p=T.morph==lower(n(m));plot(ax,T.binWidthUm(p),T.percentCompact(p),'-o','Color',c(m,:),'LineWidth',2.5,'DisplayName',n(m));end;xlabel(ax,'Fixed histogram bin width (\mum)');ylabel(ax,'Compact assemblies (%)');title(ax,'Sensitivity of H < 1.5 classification to bin width');legend(ax,'Location','best');grid(ax,'on');saveFigure(f,out,'hansen_bin_width_sensitivity',q.pngResolution);
end
function [colors,names]=plotInfo(),colors=assembly_figure_style().morphColors;names=["Surface","Molino","Pachon"];end
function saveFigure(f,out,name,res),style_assembly_figure(f);exportgraphics(f,fullfile(out,[name '.png']),'Resolution',res);exportgraphics(f,fullfile(out,[name '.svg']),'ContentType','vector');close(f);end
