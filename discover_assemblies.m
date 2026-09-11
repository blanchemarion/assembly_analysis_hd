function model=discover_assemblies(F,activity,opts,seed)
% PCA-Promax discovery using training data only.
assert(size(F,1)==size(activity,1)&&size(F,2)==size(activity,2)); rng(seed,'twister');
[Z,mu,sd]=zscore0(F); C=(Z.'*Z)/max(1,size(Z,1)-1); C=(C+C.')/2; [V,D]=eig(C,'vector'); [eigval,o]=sort(real(D),'descend'); V=real(V(:,o));
[threshold,nullEig]=shiftNull(Z,opts); pass=eigval(:)>threshold(:); firstFail=find(~pass,1,'first'); if isempty(firstFail), nPC=numel(pass); else, nPC=firstFail-1; end
% Consecutive leading ranks only: an isolated later rank can never survive.
if nPC==0, W=zeros(size(F,2),0); components=W; else, components=promax(V(:,1:nPC),opts.promaxPower); components=normalizeSigns(components); W=thresholdComponents(components,opts.loadingThreshold,opts.minAssemblySize); components=W; end
if opts.mergeSimilarAssemblies&&size(W,2)>1, [W,components]=mergeSimilar(W,components,opts); end
pruneInfo=struct('kept',true(1,size(W,2)),'sync',[],'coherence',[]);
if opts.pruneWeakAssemblies&&~isempty(W), [W,components,pruneInfo]=prune_weak_assemblies(W,components,activity,opts,seed+100000); end
model=struct('weights',W,'membership',W~=0,'components',components,'eigenvalues',eigval, ...
 'nullThresholdByRank',threshold,'nullEigenvalues',nullEig,'nLeadingPCs',nPC, ...
 'loadingThreshold',opts.loadingThreshold,'preprocessMean',mu,'preprocessStd',sd, ...
 'pcaSelectionRule','consecutive_leading_significant_ranks','pruneInfo',pruneInfo,'seed',seed,'trainFingerprint',fingerprint(F));
end
function [th,N]=shiftNull(Z,opts)
[T,Nc]=size(Z); N=zeros(opts.nNullIter,Nc); lo=max(1,round(opts.minShiftFraction*T)); hi=max(lo,T-lo); X=zeros(size(Z));
for q=1:opts.nNullIter, for c=1:Nc, if hi>lo, sh=randi([lo hi]); else, sh=lo; end, X(:,c)=circshift(Z(:,c),sh); end
 C=(X.'*X)/max(1,T-1); N(q,:)=sort(max(real(eig((C+C.')/2)),0),'descend');
end
th=prctile(N,opts.nullPercentile,1)';
end
function R=promax(X,p)
if size(X,2)<2, R=X; return, end
try, R=rotatefactors(X,'Method','promax','Power',p,'Normalize','off'); return, catch, end
[nr,nc]=size(X); Q=eye(nc); old=0;
for k=1:500, L=X*Q; target=L.^3-(1/nr)*L*diag(sum(L.^2,1)); [U,S,V]=svd(X.'*target,'econ'); Q=U*V.'; val=sum(diag(S)); if k>1&&val<old*(1+1e-6), break, end, old=val; end
L=X*Q; target=sign(L).*abs(L).^p; B=L\target; B=B./max(sqrt(sum(B.^2,1)),eps); R=X*Q*B;
end
function X=normalizeSigns(X)
for c=1:size(X,2), X(:,c)=X(:,c)/max(norm(X(:,c)),eps); if max(X(:,c))<=abs(min(X(:,c))), X(:,c)=-X(:,c); end, end
end
function W=thresholdComponents(C,th,minN)
Z=zscore0(C); keep=sum(Z>=th,1)>=minN; C=C(:,keep); Z=Z(:,keep); W=C; W(Z<th)=0; W=W./max(vecnorm(W),eps);
end
function [W,C]=mergeSimilar(W,C,o)
while size(W,2)>1
 S=abs(W.'*W); S(triu(true(size(S))))=0; [v,ix]=max(S(:)); if v<=o.similarityThreshold, break, end, [a,b]=ind2sub(size(S),ix); z=C(:,a)+sign(C(:,a)'*C(:,b))*C(:,b); z=z/max(norm(z),eps); kill=sort([a b],'descend'); W(:,kill)=[]; C(:,kill)=[]; nw=thresholdComponents(z,o.loadingThreshold,o.minAssemblySize); if ~isempty(nw), W(:,end+1)=nw; C(:,end+1)=z; end
end
end
function [Z,mu,sd]=zscore0(X), mu=mean(X,1,'omitnan'); sd=std(X,0,1,'omitnan'); sd(sd==0|~isfinite(sd))=1; Z=(X-mu)./sd; Z(~isfinite(Z))=0; end
function h=fingerprint(X), h=[size(X),sum(X(:),'omitnan'),sum(X(:).^2,'omitnan')]; end
