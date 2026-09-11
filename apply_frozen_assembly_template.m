function E=apply_frozen_assembly_template(model,activity,candidateRoiIdx)
% Apply, never refit, a training template to held-out candidate-ordered data.
if size(activity,2)~=size(model.weights,1), error('assembly:CandidateCountChanged','Candidate count changed.'); end
if ~isequal(candidateRoiIdx(:),model.candidateRoiIdx(:)), error('assembly:CandidateOrderChanged','Candidate identity/order changed.'); end
B=activity>0; pop=sum(B,2); E=zeros(size(model.weights,2),size(B,1)); P=ones(size(E)); coherence=NaN(size(model.weights,2),1);
for a=1:size(E,1)
 ids=find(model.membership(:,a)); k=sum(B(:,ids),2); den=numel(ids)+pop; good=den>0; E(a,good)=(2*k(good)./den(good))';
 if numel(ids)>1, R=corr(double(B(:,ids))); x=R(triu(true(size(R)),1)); coherence(a)=mean(x(isfinite(x)),'omitnan'); end
 for t=find(k>0)', P(a,t)=hyperUpper(size(B,2),numel(ids),pop(t),k(t)); end
end
E=struct('engagement',E,'pValue',P,'coherence',coherence,'templateFingerprint',model.trainFingerprint,'weights',model.weights,'membership',model.membership);
end
function p=hyperUpper(N,K,n,k), hi=min(K,n); x=k:hi; lp=gammaln(K+1)-gammaln(x+1)-gammaln(K-x+1)+gammaln(N-K+1)-gammaln(n-x+1)-gammaln(N-K-n+x+1)-gammaln(N+1)+gammaln(n+1)+gammaln(N-n+1); p=min(1,sum(exp(lp))); end
