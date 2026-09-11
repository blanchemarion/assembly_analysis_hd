function E=evaluate_hansen_assembly(membershipColumn,candidateRoiIdx,allCoordinates,distanceBinEdges,cutoff)
%EVALUATE_HANSEN_ASSEMBLY Map members and either calculate or cleanly exclude.
[roi,coordinates,valid]=map_assembly_members_to_all_cells(membershipColumn,candidateRoiIdx,allCoordinates);
E=struct('evaluated',false,'reason',"",'memberOriginalRoi',roi, ...
 'memberCoordinates',coordinates,'validCoordinate',valid,'metrics',struct);
if numel(roi)<2||sum(valid)<2
 E.reason=sprintf('requires >=2 members with valid coordinates; assemblySize=%d, nValid=%d',numel(roi),sum(valid));
 return
end
E.metrics=hansen_spatial_metrics(coordinates(valid,:),distanceBinEdges);
E.metrics.isCompactHansen=E.metrics.distanceEntropy<cutoff;E.evaluated=true;
end
