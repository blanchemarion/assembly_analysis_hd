function [W,C,info]=prune_weak_assemblies(W,C,A,opts,seed)
% Size-matched synchrony/coherence null computed on training activity only.
rng(seed,'twister'); B=A>0; na=size(W,2); keep=false(1,na); sync=zeros(1,na); coh=zeros(1,na); syncNull=cell(1,na); cohNull=cell(1,na);
for a=1:na
 members=find(W(:,a)~=0); sync(a)=syncMetric(B,members,opts.pruneSyncBinFrames); coh(a)=cohMetric(B,members); sn=zeros(opts.pruneNullIterations,1); cn=sn;
 for q=1:opts.pruneNullIterations, ids=randperm(size(B,2),numel(members)); sn(q)=syncMetric(B,ids,opts.pruneSyncBinFrames); cn(q)=cohMetric(B,ids); end
 syncNull{a}=sn; cohNull{a}=cn; keep(a)=sync(a)>=max(opts.pruneMinSyncFraction,prctile(sn,95))&&coh(a)>=prctile(cn,95);
end
W=W(:,keep); C=C(:,keep); info=struct('kept',keep,'sync',sync,'coherence',coh,'syncNull',{syncNull},'coherenceNull',{cohNull});
end
function v=syncMetric(B,ids,w), if isempty(ids), v=0; return, end, x=conv2(double(B(:,ids)),ones(w,1),'same')>0; v=max(sum(x,2))/numel(ids); end
function v=cohMetric(B,ids), if numel(ids)<2, v=0; return, end, R=corr(double(B(:,ids))); x=R(triu(true(size(R)),1)); x=x(isfinite(x)); if isempty(x), v=0; else, v=mean(x); end, end
