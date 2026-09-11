function result=run_assembly_batch(opts)
ensureDir(opts.processedRoot); ensureDir(opts.figureRoot); outDir=fullfile(opts.processedRoot,opts.analysisLabel); ensureDir(outDir);
manifest=build_assembly_manifest(opts); logPath=fullfile(opts.processedRoot,sprintf('processing_log_%s.txt',opts.runTag));
fid=fopen(logPath,'a'); assert(fid>0,'Cannot open processing log.'); cleaner=onCleanup(@()fclose(fid)); fprintf(fid,'\n=== %s ===\n',datestr(now,30)); %#ok<NASGU>
batchTimer=tic; nSessions=height(manifest);
fprintf("Manifest contains %d sessions. Log: %s\n",nSessions,logPath); drawnow;
for i=1:height(manifest)
 if manifest.status(i)~="pending"
  fprintf(fid,'EXCLUDED %s/%s: %s\n',manifest.morph(i),manifest.recording(i),manifest.exclusionReason(i)); fprintf("[%d/%d] EXCLUDED %s/%s: %s\n",i,nSessions,manifest.morph(i),manifest.recording(i),manifest.exclusionReason(i)); drawnow; continue
 end
 sessionTimer=tic; fprintf("[%d/%d] Processing %s/%s...\n",i,nSessions,manifest.morph(i),manifest.recording(i)); drawnow;
 try
  S=load_assembly_session(manifest(i,:),opts); manifest.nCells(i)=numel(S.candidateRoiIdx); manifest.nFrames(i)=S.nFrames; manifest.durationSec(i)=S.durationSec; manifest.fps(i)=S.fps;
  out=char(manifest.outputPath(i)); if isfile(out)&&~opts.overwriteExisting, manifest.status(i)="existing"; fprintf("  Existing output skipped (%.1f s).\n",toc(sessionTimer)); drawnow; continue, end
  seed=opts.randomSeed+i-1; model=discover_assemblies(S.pca,S.activity,opts,seed); model.candidateRoiIdx=S.candidateRoiIdx; model.candidateCellIDs=S.candidateCellIDs; model.preferredPhaseRad=S.preferredPhaseRad; model.nTotalRois=S.nTotalRois;
  activation=call_assembly_activations(model,S.activity,struct('mode','bh-fdr','pThreshold',.05,'minEngagement',.20,'minDurationFrames',1,'mergeGapFrames',0));
  provenance=struct('trainFrames',(1:S.nFrames)','evaluationFrames',(1:S.nFrames)','purpose','batch_descriptive_not_validation','seed',seed,'options',opts);
  save(out,'model','activation','provenance','S','-v7.3'); manifest.status(i)="processed"; fprintf(fid,'OK %s/%s assemblies=%d\n',manifest.morph(i),manifest.recording(i),size(model.weights,2));
  fprintf("  Complete: %d assemblies (%.1f s).\n",size(model.weights,2),toc(sessionTimer)); drawnow;
 catch ME
  manifest.status(i)="failed"; manifest.exclusionReason(i)=string(ME.message); fprintf(fid,'FAILED %s/%s %s\n%s\n',manifest.morph(i),manifest.recording(i),ME.message,getReport(ME,'extended','hyperlinks','off'));
  fprintf("  FAILED after %.1f s: %s\n",toc(sessionTimer),ME.message); drawnow;
 end
end
writetable(manifest,fullfile(opts.processedRoot,'session_manifest.csv')); save(fullfile(opts.processedRoot,'session_manifest.mat'),'manifest','opts'); result=struct('manifest',manifest,'logPath',logPath,'options',opts);
fprintf("Batch complete in %.1f min: %d processed, %d existing, %d failed.\n",toc(batchTimer)/60,sum(manifest.status=="processed"),sum(manifest.status=="existing"),sum(manifest.status=="failed")); drawnow;
end
function ensureDir(p), if ~isfolder(p), mkdir(p); end, end
