function out=call_assembly_activations(modelOrEvaluation,activityOrOpts,opts)
% Calling criterion is separate from the membership/loading threshold.
if isstruct(modelOrEvaluation)&&isfield(modelOrEvaluation,'pValue')
 E=modelOrEvaluation; opts=activityOrOpts;
else
 model=modelOrEvaluation; activity=activityOrOpts;
 if ~isfield(model,'candidateRoiIdx'), model.candidateRoiIdx=(1:size(activity,2))'; end
 E=apply_frozen_assembly_template(model,activity,model.candidateRoiIdx);
end
p=E.pValue; sig=false(size(p));
switch lower(opts.mode)
 case {'bh-fdr','bh_fdr','fdr'}
  for a=1:size(p,1), sig(a,:)=bh(p(a,:),opts.pThreshold); end
 otherwise, error('assembly:BHFDROnly','Activation inference uses BH-FDR only.');
end
active=sig&E.engagement>=opts.minEngagement; episodes=cell(size(active,1),1);
for a=1:size(active,1), active(a,:)=mergeGaps(active(a,:),opts.mergeGapFrames); active(a,:)=removeShort(active(a,:),opts.minDurationFrames); episodes{a}=runs(active(a,:)); end
out=struct('active',active,'engagement',E.engagement,'pValue',p,'episodes',{episodes}, ...
 'criteria',opts,'templateFingerprint',E.templateFingerprint,'membershipFingerprint',[size(E.membership),sum(E.membership(:))]);
end
function h=bh(p,q), [s,o]=sort(p); k=find(s<=q*(1:numel(p))/numel(p),1,'last'); h=false(size(p)); if ~isempty(k), h(o(1:k))=true; end, end
function x=mergeGaps(x,g), if g<1, return, end, z=runs(~x); for i=1:size(z,1), if z(i,1)>1&&z(i,2)<numel(x)&&z(i,2)-z(i,1)+1<=g, x(z(i,1):z(i,2))=true; end, end, end
function x=removeShort(x,n), z=runs(x); for i=1:size(z,1), if z(i,2)-z(i,1)+1<n, x(z(i,1):z(i,2))=false; end, end, end
function z=runs(x), d=diff([false,x(:)',false]); z=[find(d==1)',find(d==-1)'-1]; end
