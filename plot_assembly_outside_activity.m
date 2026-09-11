function result=plot_assembly_outside_activity(report,outDir,visibility,resolution)
% All inactive frames, including recordings without complete onset windows.
% Binary-positive frames/minute matches the properties-comparison estimator.
rows=cell(0,6);
for i=1:numel(report.activationAudits)
 r=report.activationAudits(i); A=r.audit;
 records=report.perAssembly;
 ix=find(string({records.file})==string(r.file),1);
 if isempty(ix)
  X=load(r.file,'S');fps=double(X.S.fps);
 else
  fps=records(ix).fps;
 end
 for a=1:size(A.active,1)
  outside=~A.active(a,:);n=sum(outside);
  if n==0,continue,end
  fraction=mean(A.recruitment(a,outside));
  rows(end+1,:)={string(r.morph),string(r.fish),a,n,fraction,fraction*fps*60}; %#ok<AGROW>
 end
end
assemblyTable=cell2table(rows,'VariableNames',{'morph','fish','assemblyID','nOutsideFrames','fractionActive','binaryEventsPerNeuronPerMinute'});
[G,morph,fish]=findgroups(assemblyTable.morph,assemblyTable.fish);
fractionActive=splitapply(@mean,assemblyTable.fractionActive,G);
binaryEventsPerNeuronPerMinute=splitapply(@mean,assemblyTable.binaryEventsPerNeuronPerMinute,G);
fishTable=table(morph,fish,fractionActive,binaryEventsPerNeuronPerMinute);
writetable(assemblyTable,fullfile(outDir,'outside_activity_all_frames_per_assembly.csv'));
writetable(fishTable,fullfile(outDir,'outside_activity_all_frames_per_fish.csv'));
colors=assembly_figure_style().morphColors;
order=["surface","molino","pachon"];
fig=figure('Visible',char(visibility),'Color','w','Position',[100 100 1000 450]);
layout=tiledlayout(fig,1,2,'TileSpacing','compact');
fields={'fractionActive','binaryEventsPerNeuronPerMinute'};
labels={'Fraction of assembly members active','Binary events per neuron per minute'};
summaryRows=cell(0,6);
for k=1:2
 ax=nexttile(layout);hold(ax,'on');
 for m=1:3
  v=fishTable.(fields{k})(fishTable.morph==order(m));n=numel(v);
  if n==0,continue,end
  mu=mean(v);sem=std(v)/sqrt(n);
  bar(ax,m,mu,.46,'FaceColor',colors(m,:),'FaceAlpha',.78,'EdgeColor','k','LineWidth',.8);
  scatter(ax,m+assembly_deterministic_jitter(n,.11),v,25,'MarkerFaceColor',assembly_mix_with_white(colors(m,:),.38),'MarkerEdgeColor','k','LineWidth',.45,'MarkerFaceAlpha',.82);
  errorbar(ax,m,mu,sem,'k','LineStyle','none','LineWidth',1.7,'CapSize',8);
  summaryRows(end+1,:)={order(m),string(fields{k}),n,mu,std(v),sem}; %#ok<AGROW>
 end
 xticks(ax,1:3);xticklabels(ax,{'Surface','Molino','Pachon'});
 ylabel(ax,labels{k});grid(ax,'on');ylim(ax,[0 Inf]);
end
title(layout,'All outside-episode frames | equal fish weights | mean +/- SEM');
stem=fullfile(outDir,'outside_activity_all_frames_across_morphs');
style_assembly_figure(fig);exportgraphics(fig,[stem '.png'],'Resolution',resolution);savefig(fig,[stem '.fig']);close(fig);
summaryTable=cell2table(summaryRows,'VariableNames',{'morph','metric','nFish','mean','SD','SEM'});
writetable(summaryTable,fullfile(outDir,'outside_activity_all_frames_morph_summary.csv'));
result=struct('assemblyTable',assemblyTable,'fishTable',fishTable,'summaryTable',summaryTable);
end
