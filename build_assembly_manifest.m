function manifest=build_assembly_manifest(opts)
% Match disk recordings to the independent QC keep-list by morph and rec ID.
assert(isfile(opts.qcFile),'QC keep-list not found: %s',opts.qcFile);
Q=readtable(opts.qcFile,'TextType','string'); required=["group","file","keep"];
assert(all(ismember(required,string(Q.Properties.VariableNames))),'QC file must contain group, file and keep columns.');
Q=Q(Q.keep==1,:); qMorph=arrayfun(@normMorph,Q.group); qRec=arrayfun(@recId,Q.file);
if any(qMorph==""|qRec==""), error('assembly:InvalidQCEntry','A kept QC row lacks a parseable morph or recording ID.'); end
qKey=qMorph+"__"+qRec; [u,~,g]=unique(qKey); counts=accumarray(g,1);
if any(counts>1), error('assembly:AmbiguousQCMatch','Duplicate kept QC key(s): %s',strjoin(u(counts>1),', ')); end
rows={}; found=false(height(Q),1);
for m=1:numel(opts.morphs)
 morph=string(opts.morphs{m}); d=dir(fullfile(opts.dataRoot,morph,'rec*')); d=d([d.isdir]); [~,z]=sort({d.name}); d=d(z);
 for j=1:numel(d)
  rec=string(d(j).name); p=string(fullfile(d(j).folder,d(j).name)); key=morph+"__"+recId(rec); qi=find(qKey==key);
  status="excluded"; reason="not_listed_in_QC_ENTRY_KEPT_files.csv";
  if numel(qi)>1, error('assembly:AmbiguousQCMatch','Ambiguous QC match for %s',key); end
  if numel(qi)==1, status="pending"; reason="included_qc_keep"; found(qi)=true; end
  rf=files(p,'*RASTER*.mat',{'SELECTED_ROIS','NON_SELECTED_ROIS','ANATOMICAL','WITH_ASSEMBLIES','PCA_PROMAX','ENSEMBLE','CLASSIFIED','DIRECT_MODEL'}); hf=files(p,'*candidate_neurons_clean.mat',{'CLASSIFIED','DIRECT_MODEL'});
  if status=="pending"&&numel(rf)~=1, status="excluded"; reason="expected_one_raster_found_"+numel(rf); end
  if status=="pending"&&numel(hf)~=1, status="excluded"; reason="expected_one_hdmetrics_found_"+numel(hf); end
  r=""; h=""; if numel(rf)==1,r=string(fullfile(rf.folder,rf.name));end, if numel(hf)==1,h=string(fullfile(hf.folder,hf.name));end
  out=string(fullfile(opts.processedRoot,opts.analysisLabel,morph+'__'+rec+'.mat'));
  rows(end+1,:)={morph,rec,p,r,h,h,out,NaN,NaN,NaN,NaN,status,reason,"disk"}; %#ok<AGROW>
 end
end
for i=find(~found)', warning('assembly:MissingQCRecording','Kept QC entry not found on disk: %s',qKey(i)); rows(end+1,:)={qMorph(i),qRec(i),"",string(Q.file(i)),string(Q.file(i)),string(Q.file(i)),"",NaN,NaN,NaN,NaN,"missing","qc_entry_not_found_on_disk","qc"}; end %#ok<AGROW>
manifest=cell2table(rows,'VariableNames',{'morph','recording','recordingPath','rasterPath','hdMetricsPath','candidatePath','outputPath','nCells','nFrames','durationSec','fps','status','exclusionReason','matchSource'});
if isfinite(opts.maxSessionsPerMorph)
 for m=1:numel(opts.morphs), ix=find(manifest.status=="pending"&manifest.morph==opts.morphs{m}); if numel(ix)>opts.maxSessionsPerMorph, drop=ix(opts.maxSessionsPerMorph+1:end); manifest.status(drop)="excluded"; manifest.exclusionReason(drop)="smoke_subset_limit"; end, end
end
end
function f=files(p,pat,banned), f=dir(fullfile(p,pat)); f=f(~[f.isdir]); keep=true(size(f)); for i=1:numel(f),u=upper(f(i).name);for b=1:numel(banned),keep(i)=keep(i)&&~contains(u,banned{b});end,end,f=f(keep);end
function m=normMorph(x),x=lower(string(x));if contains(x,'surface'),m="surface";elseif contains(x,'molino'),m="molino";elseif contains(x,'pachon'),m="pachon";else,m="";end,end
function r=recId(x),t=regexp(lower(char(x)),'rec0*([0-9]+)','tokens','once');if isempty(t),r="";else,r="rec"+string(str2double(t{1}));end,end
