function summary=plot_assembly_onset_diagnostics(report,outDir,visibility,resolution)
% Explain onset gating using the same windows and fish weights as the main plot.
frames=report.settings.relativeFrames(:);E=report.episodeSelection;
rows=cell(0,9);
for i=1:numel(E)
 e=E(i);t=e.retainedOnsetIndices;if isempty(t),continue,end
 j=find(string({report.activationAudits.file})==string(e.file));
 assert(isscalar(j));A=report.activationAudits(j).audit;a=e.assemblyID;
 ix=frames+t;
 fields={'recruitment','matchingIndex','active','startEligible'};
 V=zeros(numel(frames),4);
 for k=1:4
  Z=A.(fields{k});V(:,k)=mean(reshape(double(Z(a,ix)),size(ix)),2);
 end
 for k=1:numel(frames)
  rows(end+1,:)={string(e.morph),string(e.fish),a,frames(k),numel(t),V(k,1),V(k,2),V(k,3),V(k,4)}; %#ok<AGROW>
 end
end
T=cell2table(rows,'VariableNames',{'morph','fish','assemblyID','relativeFrame','nWindows','recruitment','matchingIndex','activationProbability','startEligibleProbability'});
writetable(T,fullfile(outDir,'onset_activation_diagnostics_per_assembly.csv'));
order=["surface","molino","pachon"];colors=assembly_figure_style().morphColors;
fields={'recruitment','matchingIndex','activationProbability','startEligibleProbability'};
labels={'Active member fraction','Matching index','Fraction of windows classified active','Fraction of windows meeting onset criteria'};
fig=figure('Visible',char(visibility),'Color','w','Position',[100 100 1100 750]);layout=tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
summaryRows=cell(0,7);
for k=1:4
 ax=nexttile(layout);hold(ax,'on');
 for m=1:3
  mu=nan(size(frames));se=mu;
  for f=1:numel(frames)
   Q=T(T.morph==order(m)&T.relativeFrame==frames(f),:);
   if isempty(Q),continue,end
   G=findgroups(Q.fish);v=splitapply(@mean,Q.(fields{k}),G);
   mu(f)=mean(v);se(f)=std(v)/sqrt(numel(v));
   summaryRows(end+1,:)={order(m),string(fields{k}),frames(f),numel(v),mu(f),std(v),se(f)}; %#ok<AGROW>
  end
  fill(ax,[frames;flipud(frames)],[mu-se;flipud(mu+se)],colors(m,:),'FaceAlpha',.18,'EdgeColor','none','HandleVisibility','off');
  plot(ax,frames,mu,'Color',colors(m,:),'LineWidth',2.5,'DisplayName',order(m));
 end
 if k==1,yline(ax,report.settings.onRecruitmentFraction,':k','Onset threshold','HandleVisibility','off');end
 if k==2,yline(ax,report.settings.minMatchingIndex,':k','MI threshold','HandleVisibility','off');end
 xline(ax,0,'--k','HandleVisibility','off');ylim(ax,[0 1.05]);xlim(ax,[frames(1) frames(end)]);xticks(ax,frames);
 xlabel(ax,'Frames from detected onset');ylabel(ax,labels{k});grid(ax,'on');legend(ax,'Location','best');
end
title(layout,'Onset criteria and detected state | equal fish weights | mean +/- SEM');
style_assembly_figure(fig);exportgraphics(fig,fullfile(outDir,'onset_activation_diagnostics.png'),'Resolution',resolution);close(fig);
summary=cell2table(summaryRows,'VariableNames',{'morph','metric','relativeFrame','nFish','mean','SD','SEM'});
writetable(summary,fullfile(outDir,'onset_activation_diagnostics_summary.csv'));
end
