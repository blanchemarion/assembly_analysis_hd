function M=hansen_spatial_metrics(coordinates,distanceBinEdges)
%HANSEN_SPATIAL_METRICS Unweighted pair-distance metrics with fixed bins.
coordinates=double(coordinates);distanceBinEdges=double(distanceBinEdges(:)');
assert(ismatrix(coordinates)&&size(coordinates,2)==2,'Coordinates must be N-by-2.');
assert(size(coordinates,1)>=2&&all(isfinite(coordinates),'all'),'At least two finite coordinates are required.');
assert(numel(distanceBinEdges)>=2&&distanceBinEdges(1)==0&& ...
 all(diff(distanceBinEdges)>0),'Distance-bin edges must increase strictly from zero.');
n=size(coordinates,1);distances=zeros(n*(n-1)/2,1);k=0;
for i=1:n-1
 d=sqrt(sum((coordinates(i+1:end,:)-coordinates(i,:)).^2,2));
 distances(k+(1:numel(d)))=d;k=k+numel(d);
end
counts=histcounts(distances,distanceBinEdges);
P=counts/sum(counts);positive=P>0;
entropy=-sum(P(positive).*log(P(positive)));
M=struct('distances',distances,'histogramCounts',counts, ...
 'histogramProbability',P,'meanPairwiseDistance',mean(distances), ...
 'medianPairwiseDistance',median(distances),'distanceEntropy',entropy, ...
 'isCompactHansen',entropy<1.5,'distanceBinEdges',distanceBinEdges);
end
