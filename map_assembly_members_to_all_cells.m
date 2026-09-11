function [memberOriginalRoi,memberCoordinates,validCoordinate]=map_assembly_members_to_all_cells(membershipColumn,candidateRoiIdx,allCoordinates)
%MAP_ASSEMBLY_MEMBERS_TO_ALL_CELLS Map membership rows through original ROI IDs.
membershipColumn=logical(membershipColumn(:));candidateRoiIdx=double(candidateRoiIdx(:));
allCoordinates=double(allCoordinates);
assert(numel(membershipColumn)==numel(candidateRoiIdx), ...
 'Membership rows and candidateRoiIdx must align one-to-one.');
assert(numel(unique(candidateRoiIdx))==numel(candidateRoiIdx), ...
 'candidateRoiIdx must be unique.');
assert(all(isfinite(candidateRoiIdx)&candidateRoiIdx==round(candidateRoiIdx)&candidateRoiIdx>=1), ...
 'candidateRoiIdx must contain positive finite integers.');
assert(size(allCoordinates,2)==2,'Complete ROI coordinates must be N-by-2.');
assert(all(candidateRoiIdx<=size(allCoordinates,1)), ...
 'candidateRoiIdx exceeds the complete ALL_CELLS ROI coordinate range.');
memberOriginalRoi=candidateRoiIdx(membershipColumn);
memberCoordinates=allCoordinates(memberOriginalRoi,:);
validCoordinate=all(isfinite(memberCoordinates),2);
end
