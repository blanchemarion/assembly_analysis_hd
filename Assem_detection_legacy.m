%% ========================================================================
% BATCH PCA-PROMAX ASSEMBLY DETECTION
%
% For each fish folder:
%
%   1. Find the original file ending in RASTER.mat
%   2. Find the hdMetrics_*.mat file
%   3. Read candidate cells from:
%
%          hdMetrics.r1pi.candidate_idx
%
%   4. Keep only the first 30 minutes
%   5. Keep only candidate cells
%   6. Run PCA + circular-shift null model + promax
%   7. Save assembly results inside each fish folder
%
% Destination structure:
%
%   /mnt/data3/rec_resultsV2/1surface/rec*/
%   /mnt/data3/rec_resultsV2/2molino/rec*/
%   /mnt/data3/rec_resultsV2/3pachon/rec*/
%
% Important:
%   candidate_idx must contain MATLAB-style ROI indices:
%   1, 2, 3, ..., nROIs
% ========================================================================

clear;
close all;
clc;

%% ========================= USER SETTINGS ================================

groupFolders = {
    '/mnt/disk1/rec_results_v3/1surface'
    '/mnt/disk1/rec_results_v3/2molino'
    '/mnt/disk1/rec_results_v3/3pachon'
};

% Field used to identify assemblies
dataFieldForPCA = 'deltaFoF';

% Binary activity used for assembly engagement
dataFieldForEngagement = 'raster';

% Recording segment used for assembly detection
analysisDurationMin = 30;

% Which cells enter assembly detection:
%
%   'candidates'  only hdMetrics.r1pi.candidate_idx. Assemblies then
%                 describe the internal organization of a population that
%                 is already selected for head-direction tuning.
%
%   'allCells'    every ROI in the mask. Assemblies then describe the
%                 imaged region as a whole, and head-direction content
%                 becomes something to test afterwards rather than
%                 something imposed beforehand.
%
%   'ask'         choose at runtime.
%
% These are different scientific questions, not two settings of one
% analysis, and the results are not directly comparable.
cellSelectionMode = 'ask';

switch lower(cellSelectionMode)

    case 'ask'

        choice = menu( ...
            'Which cells should enter assembly detection?', ...
            'Only head-direction candidate cells', ...
            'All cells in the mask');

        if choice == 0
            error('Cell selection was cancelled by the user.');
        end

        if choice == 1
            cellSelectionMode = 'candidates';
        else
            cellSelectionMode = 'allCells';
        end

    case {'candidates','allcells'}
        % Nothing to resolve.

    otherwise
        error('Unknown cellSelectionMode: %s',cellSelectionMode);
end

% Token written into every output name, so the two modes never overwrite
% each other.
switch lower(cellSelectionMode)
    case 'candidates'
        cellSelectionToken = 'CANDIDATES';
    case 'allcells'
        cellSelectionToken = 'ALLCELLS';
end

fprintf('\nCell selection: %s (output token %s)\n', ...
    cellSelectionMode,cellSelectionToken);

if strcmpi(cellSelectionMode,'allCells')
    fprintf(['\n[COST] The null model runs one eigendecomposition of an ' ...
        'N x N matrix\nper iteration. With the full mask N is far ' ...
        'larger than the candidate set and\ncost grows as N^3, so a run ' ...
        'that takes minutes on candidates can take hours\non all cells. ' ...
        'Consider lowering opts.nNullIter and smoke-testing one fish ' ...
        'first.\n']);
end

opts = struct();

opts.nNullIter          = 200;
opts.nullPercentile     = 95;

% How lambda_k observed is judged against the circular-shift null:
%
%   'maxPerIteration'  one threshold for all ranks, taken from the
%                      distribution of the largest null eigenvalue. This is
%                      a family-wise test calibrated for lambda_1, so it is
%                      very conservative for ranks 2 and above.
%
%   'parallelAnalysis' Horn parallel analysis: lambda_k observed is compared
%                      with the percentile of lambda_k under the null, rank
%                      by rank. Recommended.
opts.nullEigenvalueMode = 'parallelAnalysis';

% Minimum circular shift, as a fraction of the recording. Very small shifts
% barely decorrelate slow calcium transients and inflate the null.
opts.minShiftFraction = 0.10;

% Promax rotation
opts.promaxPower = 4;

% Assembly membership
opts.loadingThreshold       = 1.5;
opts.minAssemblySize        = 2;
opts.mergeSimilarAssemblies = true;
opts.similarityThreshold    = 0.60;

% Outputs
% Weak-assembly pruning, following Romano et al. 2017.
%
% After promax and merging, each assembly is compared with random cell sets
% on two statistics:
%
%   synchrony    the largest fraction of member cells simultaneously active
%                within a short bin. Answers: do these cells ever fire
%                together as a group?
%
%   coherence    the mean pairwise correlation between member cells in the
%                binary raster. Answers: are they correlated at all?
%
% An assembly is discarded when it falls below the 95th percentile of its
% own size-matched null on either statistic. Without this step a component
% that survives the eigenvalue test can still describe cells that never
% actually coactivate.
opts.pruneWeakAssemblies  = true;
opts.pruneNullIterations  = 1000;
opts.pruneSyncBinFrames   = 3;

% Romano additionally requires at least two thirds of the assembly to be
% co-active in some bin, regardless of what the null says.
opts.pruneMinSyncFraction = 2/3;

% Romano's third criterion removes assemblies whose combined z-score across
% the two statistics falls below -1, computed ACROSS the assemblies of that
% fish. It is relative rather than absolute, so it always removes something
% when several assemblies exist and it is unstable with few of them. Off by
% default; turn on only to reproduce Romano exactly.
opts.pruneRelativeZ       = false;

% Outputs
opts.saveOutput       = true;
opts.saveDiagnostics  = true;
opts.closeFigures     = true;
opts.overwriteExisting = true;

% Reproducibility
opts.randomSeed = 1;

%% ========================= FIND SESSIONS =================================

sessionList = struct( ...
    'groupName', {}, ...
    'fishName', {}, ...
    'fishPath', {}, ...
    'rasterFile', {}, ...
    'hdMetricsFile', {});

fprintf('\n============================================================\n');
fprintf('SEARCHING FISH FOLDERS\n');
fprintf('============================================================\n');

for iGroup = 1:numel(groupFolders)

    groupPath = groupFolders{iGroup};

    if ~isfolder(groupPath)
        warning('Group folder does not exist: %s', groupPath);
        continue
    end

    [~, groupName] = fileparts(groupPath);

    fishFolders = dir(fullfile(groupPath, 'rec*'));
    fishFolders = fishFolders([fishFolders.isdir]);

    fprintf('\nGroup: %s\n', groupName);
    fprintf('Fish folders found: %d\n', numel(fishFolders));

    for iFish = 1:numel(fishFolders)

        fishName = fishFolders(iFish).name;
        fishPath = fullfile(fishFolders(iFish).folder, fishName);

        fprintf('\n  Fish: %s\n', fishName);

        %% ---------------- Find original RASTER file ----------------

        rasterMatches = dir(fullfile(fishPath, '*RASTER.mat'));
        rasterMatches = rasterMatches(~[rasterMatches.isdir]);

        % Exclude outputs generated by assembly analyses
        excludeRaster = false(size(rasterMatches));

        for k = 1:numel(rasterMatches)

            upperName = upper(rasterMatches(k).name);

            excludeRaster(k) = ...
                contains(upperName, 'WITH_ASSEMBLIES') || ...
                contains(upperName, 'PCA_PROMAX')      || ...
                contains(upperName, 'ENSEMBLE');
        end

        rasterMatches = rasterMatches(~excludeRaster);

        if isempty(rasterMatches)
            fprintf('    [SKIP] No original *RASTER.mat file found.\n');
            continue
        end

        if numel(rasterMatches) > 1

            fprintf('    [SKIP] Multiple original *RASTER.mat files found:\n');

            for k = 1:numel(rasterMatches)
                fprintf('      %s\n', rasterMatches(k).name);
            end

            continue
        end

        rasterFile = fullfile( ...
            rasterMatches(1).folder, ...
            rasterMatches(1).name);

        %% ---------------- Find hdMetrics file -----------------------

        hdMatches = dir(fullfile(fishPath, 'hdMetrics_*.mat'));
        hdMatches = hdMatches(~[hdMatches.isdir]);

        % Exclude derived files if any
        excludeHd = false(size(hdMatches));

        for k = 1:numel(hdMatches)

            upperName = upper(hdMatches(k).name);

            excludeHd(k) = ...
                contains(upperName, 'PCA_PROMAX') || ...
                contains(upperName, 'WITH_ASSEMBLIES');
        end

        hdMatches = hdMatches(~excludeHd);

        if isempty(hdMatches)
            fprintf('    [SKIP] No hdMetrics_*.mat file found.\n');
            continue
        end

        if numel(hdMatches) > 1

            fprintf('    [SKIP] Multiple hdMetrics_*.mat files found:\n');

            for k = 1:numel(hdMatches)
                fprintf('      %s\n', hdMatches(k).name);
            end

            continue
        end

        hdMetricsFile = fullfile( ...
            hdMatches(1).folder, ...
            hdMatches(1).name);

        %% ---------------- Register session --------------------------

        newSession.groupName    = groupName;
        newSession.fishName     = fishName;
        newSession.fishPath     = fishPath;
        newSession.rasterFile   = rasterFile;
        newSession.hdMetricsFile = hdMetricsFile;

        sessionList(end+1) = newSession; %#ok<SAGROW>

        fprintf('    Raster    : %s\n', rasterMatches(1).name);
        fprintf('    hdMetrics : %s\n', hdMatches(1).name);
    end
end

nSessions = numel(sessionList);

fprintf('\n============================================================\n');
fprintf('VALID SESSIONS FOUND: %d\n', nSessions);
fprintf('============================================================\n');

if nSessions == 0
    error('No valid fish sessions were found.');
end

%% ========================= BATCH SUMMARY ================================

batchSummary = table();

batchSummary.group             = strings(nSessions,1);
batchSummary.fish              = strings(nSessions,1);
batchSummary.rasterFile        = strings(nSessions,1);
batchSummary.hdMetricsFile     = strings(nSessions,1);
batchSummary.outputFile        = strings(nSessions,1);

batchSummary.success           = false(nSessions,1);
batchSummary.errorMessage      = strings(nSessions,1);

batchSummary.nTotalFrames      = NaN(nSessions,1);
batchSummary.nFramesUsed       = NaN(nSessions,1);
batchSummary.durationUsedMin   = NaN(nSessions,1);

batchSummary.nTotalROIs        = NaN(nSessions,1);
batchSummary.nCandidateCells   = NaN(nSessions,1);
batchSummary.nAssemblies       = NaN(nSessions,1);

batchSummary.lambdaThreshold   = NaN(nSessions,1);
batchSummary.loadingThreshold  = NaN(nSessions,1);

batchSummary.assemblySize      = cell(nSessions,1);

%% ========================= PROCESS EACH FISH =============================

for iSession = 1:nSessions

    groupName     = sessionList(iSession).groupName;
    fishName      = sessionList(iSession).fishName;
    fishPath      = sessionList(iSession).fishPath;
    rasterFile    = sessionList(iSession).rasterFile;
    hdMetricsFile = sessionList(iSession).hdMetricsFile;

    batchSummary.group(iSession)         = string(groupName);
    batchSummary.fish(iSession)          = string(fishName);
    batchSummary.rasterFile(iSession)    = string(rasterFile);
    batchSummary.hdMetricsFile(iSession) = string(hdMetricsFile);

    fprintf('\n\n============================================================\n');
    fprintf('SESSION %d / %d\n', iSession, nSessions);
    fprintf('Group : %s\n', groupName);
    fprintf('Fish  : %s\n', fishName);
    fprintf('============================================================\n');

    try

        %% ===================== LOAD RASTER ===========================

        R = load(rasterFile);

        if ~isfield(R, dataFieldForPCA)
            error('Raster file does not contain "%s".', ...
                dataFieldForPCA);
        end

        if ~isfield(R, dataFieldForEngagement)
            error('Raster file does not contain "%s".', ...
                dataFieldForEngagement);
        end

        %% ===================== LOAD CANDIDATES =======================

H = load(hdMetricsFile);

        if ~isfield(H, 'hdMetrics') || ...
                ~isstruct(H.hdMetrics)

            error('File does not contain a valid hdMetrics structure.');
        end

        switch lower(cellSelectionMode)

            case 'candidates'

                if ~isfield(H.hdMetrics, 'r1pi') || ...
                        ~isstruct(H.hdMetrics.r1pi)

                    error('Missing hdMetrics.r1pi.');
                end

                if ~isfield(H.hdMetrics.r1pi, 'candidate_idx')

                    error('Missing hdMetrics.r1pi.candidate_idx.');
                end

                candidateCellIds = double( ...
                    H.hdMetrics.r1pi.candidate_idx(:));

                candidateCellIds = candidateCellIds( ...
                    isfinite(candidateCellIds));

                candidateCellIds = round(candidateCellIds);

                candidateCellIds = unique( ...
                    candidateCellIds, ...
                    'stable');

                if isempty(candidateCellIds)
                    error('hdMetrics.r1pi.candidate_idx is empty.');
                end

                fprintf('Candidate cells listed in hdMetrics: %d\n', ...
                    numel(candidateCellIds));

            case 'allcells'

                % The ROI count is only known once the raster has been
                % oriented, so the list is filled in below.
                candidateCellIds = [];

                fprintf('Using every ROI in the mask.\n');
        end

        %% ===================== FIND FPS ==============================

        fps = findImagingFps(R);

        if ~isfinite(fps) || fps <= 0

            error(['Could not determine imaging fps. ', ...
                'The first 30 minutes cannot be selected safely.']);
        end

        fprintf('Imaging fps: %.6f Hz\n', fps);

        %% ===================== ORIENT MATRICES =======================

        F_pca_all = makeTimeByNeuronMatrix( ...
            R.(dataFieldForPCA), ...
            dataFieldForPCA);

        F_eng_all = makeTimeByNeuronMatrix( ...
            R.(dataFieldForEngagement), ...
            dataFieldForEngagement);

        if ~isequal(size(F_pca_all), size(F_eng_all))

            error(['deltaFoF and raster have different oriented sizes.\n' ...
                'PCA matrix       : %d x %d\n' ...
                'Engagement matrix: %d x %d'], ...
                size(F_pca_all,1), size(F_pca_all,2), ...
                size(F_eng_all,1), size(F_eng_all,2));
        end

        nTotalFrames = size(F_pca_all,1);
        nTotalROIs   = size(F_pca_all,2);

        fprintf('Complete recording: %d frames x %d ROIs\n', ...
            nTotalFrames, nTotalROIs);


        %% ===================== VALIDATE CANDIDATES ===================

        if strcmpi(cellSelectionMode,'allCells')

            candidateCellIds = (1:nTotalROIs).';

            fprintf('Cells entering the analysis: %d (full mask)\n', ...
                nTotalROIs);
        end

        validCandidate = ...
            candidateCellIds >= 1 & ...
            candidateCellIds <= nTotalROIs;

        invalidCandidateIds = ...
            candidateCellIds(~validCandidate);

        candidateCellIds = ...
            candidateCellIds(validCandidate);

        if ~isempty(invalidCandidateIds)

            warning(['Removed %d candidate IDs outside ', ...
                'the valid ROI range 1:%d.'], ...
                numel(invalidCandidateIds), ...
                nTotalROIs);

            fprintf('Invalid candidate IDs:\n');
            disp(invalidCandidateIds(:).');
        end

        nCandidateCells = numel(candidateCellIds);

        if nCandidateCells < 2
            error(['Only %d valid candidate cells were found. ', ...
                'At least two are required.'], ...
                nCandidateCells);
        end

        fprintf('Valid candidate cells: %d / %d total ROIs\n', ...
            nCandidateCells, nTotalROIs);

        %% ===================== FIRST 30 MINUTES ======================

        analysisDurationSec = analysisDurationMin * 60;

        nRequestedFrames = round( ...
            analysisDurationSec * fps);

        nFramesUsed = min( ...
            nRequestedFrames, ...
            nTotalFrames);

        if nTotalFrames < nRequestedFrames

            warning(['Recording is shorter than %.1f minutes. ', ...
                'Using all available data: %.3f minutes.'], ...
                analysisDurationMin, ...
                nTotalFrames / fps / 60);
        end

        frameIdx = 1:nFramesUsed;

        durationUsedMin = ...
            nFramesUsed / fps / 60;

        fprintf('Frames used: %d / %d\n', ...
            nFramesUsed, nTotalFrames);

        fprintf('Duration used: %.3f minutes\n', ...
            durationUsedMin);

        %% ===================== SELECT CANDIDATES =====================

        % Final matrices are:
        %
        %   time x candidate cell
        %
        % Non-candidate cells are completely excluded from PCA,
        % promax rotation, assembly definition and engagement.

        F_pca = F_pca_all( ...
            frameIdx, ...
            candidateCellIds);

        F_eng = F_eng_all( ...
            frameIdx, ...
            candidateCellIds);

        fprintf('Final PCA matrix: %d frames x %d candidate cells\n', ...
            size(F_pca,1), ...
            size(F_pca,2));

        %% ===================== REMOVE INVALID CELLS ==================

        % Candidate cells with no variability cannot contribute to PCA.
        % They are removed only from the PCA analysis, while their original
        % ROI identities are preserved.

        candidateStd = std(F_pca,0,1,'omitnan');

        validVariance = ...
            isfinite(candidateStd) & ...
            candidateStd > 0;

        excludedZeroVarianceIds = ...
            candidateCellIds(~validVariance);

        if any(~validVariance)

            warning(['Removed %d candidate cells with zero or invalid ', ...
                'variance during the first 30 minutes.'], ...
                sum(~validVariance));
        end

        candidateCellIdsUsed = ...
            candidateCellIds(validVariance);

        F_pca = F_pca(:,validVariance);
        F_eng = F_eng(:,validVariance);

        nCandidateCellsUsed = numel(candidateCellIdsUsed);

        if nCandidateCellsUsed < 2

            error(['Only %d candidate cells with valid variance ', ...
                'remain after filtering.'], ...
                nCandidateCellsUsed);
        end

        fprintf('Candidate cells entering PCA: %d\n', ...
            nCandidateCellsUsed);

        %% ===================== RUN PCA-PROMAX ========================

        rng(opts.randomSeed + iSession - 1);

        [assembliesCandidate, ...
            engagement, ...
            engSignificance, ...
            info] = ...
            pcaPromaxEnsembles( ...
                F_pca, ...
                F_eng, ...
                opts);

        nAssemblies = ...
            size(assembliesCandidate,2);

        assemblySize = ...
            sum(assembliesCandidate ~= 0,1);

        fprintf('\nDetected assemblies: %d\n', ...
            nAssemblies);

        if nAssemblies > 0

            fprintf('Assembly sizes:\n');
            disp(assemblySize);
        end

        %% ===================== MAP TO ALL ROIs =======================

        % assembliesCandidate:
        %
        %   nCandidateCellsUsed x nAssemblies
        %
        % assembliesAllROIs:
        %
        %   nTotalROIs x nAssemblies
        %
        % Non-candidate cells are always zero.

        assembliesAllROIs = zeros( ...
            nTotalROIs, ...
            nAssemblies);

        assembliesAllROIs( ...
            candidateCellIdsUsed,:) = ...
            assembliesCandidate;

        assemblyMembershipCandidate = ...
            assembliesCandidate ~= 0;

        assemblyMembershipAllROIs = ...
            assembliesAllROIs ~= 0;

        %% ===================== ASSEMBLY MEMBERS ======================

        assemblyMembersCandidatePosition = ...
            cell(1,nAssemblies);

        assemblyMembersOriginalRoiIdx = ...
            cell(1,nAssemblies);

        for a = 1:nAssemblies

            % Position inside candidateCellIdsUsed
            assemblyMembersCandidatePosition{a} = ...
                find(assemblyMembershipCandidate(:,a));

            % Original ROI index in the complete raster
            assemblyMembersOriginalRoiIdx{a} = ...
                candidateCellIdsUsed( ...
                    assemblyMembersCandidatePosition{a});
        end

%% ========================================================================
% ANATOMICAL ASSEMBLY MAPS — 3 x 3 GRID
%
% Requires:
%   H.hdMetrics.roiImage
%   H.hdMetrics.centroids
%   assemblyMembersOriginalRoiIdx
%   nAssemblies
%   fishName
%   fishPath
%
% Each subplot shows:
%   - Anatomical field of view
%   - All candidate cells as small gray points
%   - Assembly members highlighted in yellow
%
% If there are more than 9 assemblies, multiple figures are generated.
% ========================================================================


if ~isfield(H.hdMetrics, 'roiImage') || ...
        isempty(H.hdMetrics.roiImage)

    warning('[ANATOMY] Missing hdMetrics.roiImage for %s.', ...
        fishName);

elseif ~isfield(H.hdMetrics, 'centroids') || ...
        isempty(H.hdMetrics.centroids)

    warning('[ANATOMY] Missing hdMetrics.centroids for %s.', ...
        fishName);

elseif nAssemblies == 0

    fprintf('[ANATOMY] No assemblies to plot for %s.\n', ...
        fishName);

else

    %% ---------------- Load anatomical information ------------------

    roiImage = H.hdMetrics.roiImage;
    centroids = double(H.hdMetrics.centroids);

    % Expected centroid format:
    %
    %   nROIs x 2
    %   column 1 = x
    %   column 2 = y
    %
    % Transpose automatically if stored as 2 x nROIs.
    if size(centroids,2) ~= 2 && size(centroids,1) == 2
        centroids = centroids.';
    end

    if size(centroids,2) < 2
        warning(['[ANATOMY] hdMetrics.centroids does not have ', ...
            'two coordinate columns for %s.'], ...
            fishName);

    else

        centroids = centroids(:,1:2);

        nCentroidROIs = size(centroids,1);

        fprintf('[ANATOMY] roiImage size: %s\n', ...
            mat2str(size(roiImage)));

        fprintf('[ANATOMY] Centroids: %d ROIs\n', ...
            nCentroidROIs);

        %% ---------------- Validate candidate indices ----------------

        candidateIdsForPlot = candidateCellIdsUsed(:);

        validCandidateForPlot = ...
            candidateIdsForPlot >= 1 & ...
            candidateIdsForPlot <= nCentroidROIs;

        if any(~validCandidateForPlot)

            warning(['[ANATOMY] %d candidate IDs are outside ', ...
                'the centroid range 1:%d.'], ...
                sum(~validCandidateForPlot), ...
                nCentroidROIs);
        end

        candidateIdsForPlot = ...
            candidateIdsForPlot(validCandidateForPlot);

        %% ---------------- Output folder -----------------------------

       anatomyOutputFolder = fullfile( ...
            fishPath, ...
            ['PCA_PROMAX_ASSEMBLY_ANATOMY_' cellSelectionToken]);

        if ~isfolder(anatomyOutputFolder)
            mkdir(anatomyOutputFolder);
        end

        %% ---------------- Figure settings ---------------------------

        assembliesPerFigure = 9;

        nAnatomyFigures = ceil( ...
            nAssemblies / assembliesPerFigure);

        markerSizeAllCandidates = 12;
        markerSizeAssembly      = 45;

        for iAnatomyFigure = 1:nAnatomyFigures

            firstAssembly = ...
                (iAnatomyFigure-1) * assembliesPerFigure + 1;

            lastAssembly = min( ...
                iAnatomyFigure * assembliesPerFigure, ...
                nAssemblies);

            assembliesInFigure = ...
                firstAssembly:lastAssembly;

            figAnatomy = figure( ...
                'Color','w', ...
                'Name',sprintf( ...
                    '%s assemblies %d-%d', ...
                    fishName, ...
                    firstAssembly, ...
                    lastAssembly), ...
                'Position',[50 50 1500 1400]);

            anatomyLayout = tiledlayout( ...
                figAnatomy, ...
                3,3, ...
                'TileSpacing','compact', ...
                'Padding','compact');

            for iPanel = 1:assembliesPerFigure

                ax = nexttile(anatomyLayout);

                if iPanel > numel(assembliesInFigure)

                    axis(ax,'off');
                    continue
                end

                assemblyId = ...
                    assembliesInFigure(iPanel);

                memberOriginalIds = ...
                    assemblyMembersOriginalRoiIdx{assemblyId};

                memberOriginalIds = ...
                    double(memberOriginalIds(:));

                validMemberIds = ...
                    isfinite(memberOriginalIds) & ...
                    memberOriginalIds >= 1 & ...
                    memberOriginalIds <= nCentroidROIs;

                memberOriginalIds = ...
                    round(memberOriginalIds(validMemberIds));

                %% ------------ Show field of view --------------------

                showRoiImageForAssemblies(ax, roiImage);

                hold(ax,'on');

                %% ------------ All candidate cells -------------------

                if ~isempty(candidateIdsForPlot)

                    candidateXY = ...
                        centroids(candidateIdsForPlot,:);

                    validCandidateXY = ...
                        all(isfinite(candidateXY),2);

                    candidateXY = ...
                        candidateXY(validCandidateXY,:);

                    scatter( ...
                        ax, ...
                        candidateXY(:,1), ...
                        candidateXY(:,2), ...
                        markerSizeAllCandidates, ...
                        [0.65 0.65 0.65], ...
                        'filled', ...
                        'MarkerFaceAlpha',0.45, ...
                        'MarkerEdgeColor','none');
                end

                %% ------------ Assembly cells in yellow -------------

                if ~isempty(memberOriginalIds)

                    memberXY = ...
                        centroids(memberOriginalIds,:);

                    validMemberXY = ...
                        all(isfinite(memberXY),2);

                    memberXY = ...
                        memberXY(validMemberXY,:);

                    % Black border improves visibility on bright images
                    scatter( ...
                        ax, ...
                        memberXY(:,1), ...
                        memberXY(:,2), ...
                        markerSizeAssembly + 35, ...
                        'k', ...
                        'filled', ...
                        'MarkerEdgeColor','none');

                    scatter( ...
                        ax, ...
                        memberXY(:,1), ...
                        memberXY(:,2), ...
                        markerSizeAssembly, ...
                        [1 1 0], ...
                        'filled', ...
                        'MarkerEdgeColor','k', ...
                        'LineWidth',0.4);
                end

                %% ------------ Panel formatting ----------------------

                axis(ax,'image');
                axis(ax,'off');

                title( ...
                    ax, ...
                    sprintf( ...
                        'Assembly %d | n = %d', ...
                        assemblyId, ...
                        numel(memberOriginalIds)), ...
                    'Interpreter','none', ...
                    'FontWeight','bold', ...
                    'FontSize',11);

                hold(ax,'off');
            end

       sgtitle( ...
    anatomyLayout, ...
    sprintf([ ...
        '%s — PCA-promax assembly anatomy (%s)\n' ...
        'Yellow: assembly cells | Gray: all candidate cells'], ...
        fishName, ...
        cellSelectionToken), ...
    'Interpreter','none', ...
    'FontWeight','bold', ...
    'FontSize',15);

            %% ---------------- Save figure ----------------------------

           anatomyBaseName = sprintf( ...
                '%s_%s_PCA_PROMAX_assembly_anatomy_%02d_to_%02d', ...
                fishName, ...
                cellSelectionToken, ...
                firstAssembly, ...
                lastAssembly);

            anatomyPngFile = fullfile( ...
                anatomyOutputFolder, ...
                [anatomyBaseName '.png']);

            anatomyFigFile = fullfile( ...
                anatomyOutputFolder, ...
                [anatomyBaseName '.fig']);

            anatomyPdfFile = fullfile( ...
                anatomyOutputFolder, ...
                [anatomyBaseName '.pdf']);

            try
                exportgraphics( ...
                    figAnatomy, ...
                    anatomyPngFile, ...
                    'Resolution',300);
            catch
                saveas(figAnatomy, anatomyPngFile);
            end

            try
                exportgraphics( ...
                    figAnatomy, ...
                    anatomyPdfFile, ...
                    'ContentType','vector');
            catch
            end

            try
                savefig(figAnatomy, anatomyFigFile);
            catch
            end

            fprintf(['[ANATOMY] Saved assemblies %d-%d:\n', ...
                '  %s\n'], ...
                firstAssembly, ...
                lastAssembly, ...
                anatomyPngFile);

            if opts.closeFigures
                close(figAnatomy);
            end
        end
    end
end

        %% ===================== PRIMARY MEMBERSHIP ====================

        assemblyPrimaryCandidate = ...
            zeros(nCandidateCellsUsed,1);

        if nAssemblies > 0

            [maximumWeight, maximumAssembly] = ...
                max(abs(assembliesCandidate),[],2);

            assignedCells = ...
                maximumWeight > 0;

            assemblyPrimaryCandidate(assignedCells) = ...
                maximumAssembly(assignedCells);
        end

        assemblyPrimaryAllROIs = ...
            zeros(nTotalROIs,1);

        assemblyPrimaryAllROIs( ...
            candidateCellIdsUsed) = ...
            assemblyPrimaryCandidate;

        assemblyMultiplicityCandidate = ...
            sum(assemblyMembershipCandidate,2);

        assemblyMultiplicityAllROIs = ...
            zeros(nTotalROIs,1);

        assemblyMultiplicityAllROIs( ...
            candidateCellIdsUsed) = ...
            assemblyMultiplicityCandidate;

        %% ===================== ASSEMBLY DELTAF/F =====================

        % F_pca is:
        %
        %   time x candidate cell
        %
        % Convert it to:
        %
        %   candidate cell x time

        candidateDeltaFoF = F_pca.';

        assemblyMeanDeltaFoF = ...
            NaN(nAssemblies,nFramesUsed);

        assemblyMedianDeltaFoF = ...
            NaN(nAssemblies,nFramesUsed);

        assemblyWeightedDeltaFoF = ...
            NaN(nAssemblies,nFramesUsed);

        for a = 1:nAssemblies

            members = ...
                assemblyMembersCandidatePosition{a};

            if isempty(members)
                continue
            end

            memberTraces = ...
                candidateDeltaFoF(members,:);

            assemblyMeanDeltaFoF(a,:) = ...
                mean(memberTraces,1,'omitnan');

            assemblyMedianDeltaFoF(a,:) = ...
                median(memberTraces,1,'omitnan');

            weights = abs( ...
                assembliesCandidate(members,a));

            weights = weights(:);

            if sum(weights) > 0

                weights = weights ./ sum(weights);

                assemblyWeightedDeltaFoF(a,:) = ...
                    sum(memberTraces .* weights,1,'omitnan');
            end
        end

        %% ===================== METADATA ==============================

        assemblyInfo = struct();

        assemblyInfo.method = 'PCA-promax';

        assemblyInfo.groupName = groupName;
        assemblyInfo.fishName  = fishName;

        assemblyInfo.rasterFile    = rasterFile;
        assemblyInfo.hdMetricsFile = hdMetricsFile;

       assemblyInfo.cellSelectionMode = cellSelectionMode;

        if strcmpi(cellSelectionMode,'candidates')
            assemblyInfo.candidateField = ...
                'hdMetrics.r1pi.candidate_idx';
        else
            assemblyInfo.candidateField = 'all ROIs in the mask';
        end

        assemblyInfo.analysisDurationRequestedMin = ...
            analysisDurationMin;

        assemblyInfo.analysisDurationUsedMin = ...
            durationUsedMin;

        assemblyInfo.fps = fps;

        assemblyInfo.nTotalFrames = ...
            nTotalFrames;

        assemblyInfo.nFramesUsed = ...
            nFramesUsed;

        assemblyInfo.nTotalROIs = ...
            nTotalROIs;

        assemblyInfo.nCandidatesInHdMetrics = ...
            numel(candidateCellIds);

        assemblyInfo.nCandidatesUsedForPCA = ...
            nCandidateCellsUsed;

        assemblyInfo.nAssemblies = ...
            nAssemblies;

        assemblyInfo.assemblySize = ...
            assemblySize;

        assemblyInfo.note = [ ...
            'Assemblies were detected only among candidate cells. ', ...
            'Non-candidate cells were excluded before PCA. ', ...
            'Only the first 30 minutes were used.'];

        %% ===================== SAVE COMPACT OUTPUT ===================

        [~, rasterBaseName, ~] = ...
            fileparts(rasterFile);

       compactOutputFile = fullfile( ...
            fishPath, ...
            [rasterBaseName ...
            '_FIRST30MIN_' cellSelectionToken ...
            '_PCA_PROMAX_ENSEMBLES.mat']);

        if isfile(compactOutputFile) && ...
                ~opts.overwriteExisting

            fprintf(['Compact output already exists and was ', ...
                'not overwritten:\n  %s\n'], ...
                compactOutputFile);

        elseif opts.saveOutput

            save(compactOutputFile, ...
                'assembliesCandidate', ...
                'assembliesAllROIs', ...
                'assemblyMembershipCandidate', ...
                'assemblyMembershipAllROIs', ...
                'assemblyMembersCandidatePosition', ...
                'assemblyMembersOriginalRoiIdx', ...
                'assemblyPrimaryCandidate', ...
                'assemblyPrimaryAllROIs', ...
                'assemblyMultiplicityCandidate', ...
                'assemblyMultiplicityAllROIs', ...
                'assemblySize', ...
                'assemblyMeanDeltaFoF', ...
                'assemblyMedianDeltaFoF', ...
                'assemblyWeightedDeltaFoF', ...
                'engagement', ...
                'engSignificance', ...
                'candidateCellIds', ...
                'candidateCellIdsUsed', ...
                'excludedZeroVarianceIds', ...
                'frameIdx', ...
                'fps', ...
                'analysisDurationMin', ...
                'durationUsedMin', ...
                'assemblyInfo', ...
                'info', ...
                'opts', ...
                'rasterFile', ...
                'hdMetricsFile', ...
                'dataFieldForPCA', ...
                'dataFieldForEngagement', ...
                'cellSelectionMode', ...
                '-v7.3');

            fprintf('Saved compact output:\n  %s\n', ...
                compactOutputFile);
        end

        %% ===================== SAVE ORIGINAL + ASSEMBLIES ============

        rasterWithAssembliesFile = fullfile( ...
            fishPath, ...
            [rasterBaseName ...
            '_FIRST30MIN_' cellSelectionToken ...
            '_WITH_ASSEMBLIES.mat']);

        if isfile(rasterWithAssembliesFile) && ...
                ~opts.overwriteExisting

            fprintf(['Raster-like output already exists and was ', ...
                'not overwritten:\n  %s\n'], ...
                rasterWithAssembliesFile);

        elseif opts.saveOutput

            % Add assembly fields to original raster variables
            R.assembliesCandidate = ...
                assembliesCandidate;

            R.assembliesAllROIs = ...
                assembliesAllROIs;

            R.assemblyMembershipCandidate = ...
                assemblyMembershipCandidate;

            R.assemblyMembershipAllROIs = ...
                assemblyMembershipAllROIs;

            R.assemblyMembersCandidatePosition = ...
                assemblyMembersCandidatePosition;

            R.assemblyMembersOriginalRoiIdx = ...
                assemblyMembersOriginalRoiIdx;

            R.assemblyPrimaryCandidate = ...
                assemblyPrimaryCandidate;

            R.assemblyPrimaryAllROIs = ...
                assemblyPrimaryAllROIs;

            R.assemblyMultiplicityCandidate = ...
                assemblyMultiplicityCandidate;

            R.assemblyMultiplicityAllROIs = ...
                assemblyMultiplicityAllROIs;

            R.assemblySize = ...
                assemblySize;

            R.assemblyMeanDeltaFoF = ...
                assemblyMeanDeltaFoF;

            R.assemblyMedianDeltaFoF = ...
                assemblyMedianDeltaFoF;

            R.assemblyWeightedDeltaFoF = ...
                assemblyWeightedDeltaFoF;

            R.assemblyEngagement = ...
                engagement;

            R.assemblyEngagementSignificance = ...
                engSignificance;

            R.assemblyCandidateCellIds = ...
                candidateCellIds;

            R.assemblyCandidateCellIdsUsed = ...
                candidateCellIdsUsed;

            R.assemblyExcludedZeroVarianceCandidateIds = ...
                excludedZeroVarianceIds;

            R.assemblyAnalysisFrameIdx = ...
                frameIdx;

            R.assemblyAnalysisFps = ...
                fps;

            R.assemblyInfo = ...
                assemblyInfo;

            R.pcaPromaxInfo = ...
                info;

            R.pcaPromaxOptions = ...
                opts;

            R.originalRasterFile = ...
                rasterFile;

            R.assemblyHdMetricsFile = ...
                hdMetricsFile;

            save(rasterWithAssembliesFile, ...
                '-struct', ...
                'R', ...
                '-v7.3');

            fprintf('Saved raster-like output:\n  %s\n', ...
                rasterWithAssembliesFile);
        end

        %% ===================== DIAGNOSTIC FIGURE =====================

        if opts.saveDiagnostics

          diagnosticFolder = fullfile( ...
                fishPath, ...
                [rasterBaseName ...
                '_FIRST30MIN_' cellSelectionToken ...
                '_PCA_PROMAX_DIAGNOSTICS']);

            if ~isfolder(diagnosticFolder)
                mkdir(diagnosticFolder);
            end

            figH = plotPcaPromaxDiagnostics( ...
                info, ...
                assembliesCandidate, ...
                assemblySize, ...
                fishName);

           diagnosticPng = fullfile( ...
                diagnosticFolder, ...
                [fishName '_' cellSelectionToken ...
                '_PCA_PROMAX_diagnostics.png']);

            diagnosticFig = fullfile( ...
                diagnosticFolder, ...
                [fishName '_' cellSelectionToken ...
                '_PCA_PROMAX_diagnostics.fig']);

            try
                exportgraphics( ...
                    figH, ...
                    diagnosticPng, ...
                    'Resolution',200);
            catch
                saveas(figH, diagnosticPng);
            end

            try
                savefig(figH, diagnosticFig);
            catch
            end

            if opts.closeFigures
                close(figH);
            end
        end

        %% ===================== UPDATE SUMMARY ========================

        batchSummary.success(iSession) = true;

        batchSummary.outputFile(iSession) = ...
            string(compactOutputFile);

        batchSummary.nTotalFrames(iSession) = ...
            nTotalFrames;

        batchSummary.nFramesUsed(iSession) = ...
            nFramesUsed;

        batchSummary.durationUsedMin(iSession) = ...
            durationUsedMin;

        batchSummary.nTotalROIs(iSession) = ...
            nTotalROIs;

        batchSummary.nCandidateCells(iSession) = ...
            nCandidateCellsUsed;

        batchSummary.nAssemblies(iSession) = ...
            nAssemblies;

        batchSummary.lambdaThreshold(iSession) = ...
            info.lambda_TH;

        batchSummary.loadingThreshold(iSession) = ...
            info.load_TH;

        batchSummary.assemblySize{iSession} = ...
            assemblySize;

        fprintf('\nSession completed successfully.\n');

    catch ME

        batchSummary.success(iSession) = false;

        batchSummary.errorMessage(iSession) = ...
            string(ME.message);

        fprintf('\nSESSION FAILED\n');
        fprintf('%s\n', ME.message);

        warning('Failed session: %s', fishName);
    end

    clear R H F_pca_all F_eng_all F_pca F_eng ...
        assembliesCandidate assembliesAllROIs ...
        engagement engSignificance info;

    if opts.closeFigures
        close all;
    end
end

%% ========================= SAVE BATCH SUMMARY ===========================

summaryMatFile = fullfile( ...
    '/mnt/disk1/rec_results_v3', ...
    ['FIRST30MIN_' cellSelectionToken ...
    '_PCA_PROMAX_BATCH_SUMMARY.mat']);

save(summaryMatFile, ...
    'batchSummary', ...
    'opts', ...
    'analysisDurationMin', ...
    '-v7.3');

summaryCsvFile = fullfile( ...
    '/mnt/disk1/rec_results_v3', ...
    ['FIRST30MIN_' cellSelectionToken ...
    '_PCA_PROMAX_BATCH_SUMMARY.csv']);

summaryForCsv = batchSummary;

summaryForCsv.assemblySize = cellfun( ...
    @(x) strjoin(string(x),';'), ...
    summaryForCsv.assemblySize, ...
    'UniformOutput',false);

writetable(summaryForCsv, summaryCsvFile);

fprintf('\n\n============================================================\n');
fprintf('BATCH COMPLETE\n');
fprintf('============================================================\n');
fprintf('Sessions found   : %d\n', nSessions);
fprintf('Successful       : %d\n', sum(batchSummary.success));
fprintf('Failed           : %d\n', sum(~batchSummary.success));
fprintf('Total assemblies : %.0f\n', ...
    sum(batchSummary.nAssemblies,'omitnan'));

fprintf('\nSummary MAT:\n  %s\n', summaryMatFile);
fprintf('Summary CSV:\n  %s\n', summaryCsvFile);
fprintf('============================================================\n');

%% ========================================================================
% LOCAL FUNCTIONS
% ========================================================================

function fps = findImagingFps(S)
% Search for imaging fps in common fields.

    fps = NaN;

    candidateNames = {
        'fps'
        'fps_im'
        'frameRate'
        'framerate'
        'samplingRate'
    };

    if isfield(S,'params') && isstruct(S.params)

        for i = 1:numel(candidateNames)

            fieldName = candidateNames{i};

            if isfield(S.params,fieldName)

                value = double(S.params.(fieldName));

                if isscalar(value) && ...
                        isfinite(value) && ...
                        value > 0

                    fps = value;
                    return
                end
            end
        end
    end

    for i = 1:numel(candidateNames)

        fieldName = candidateNames{i};

        if isfield(S,fieldName)

            value = double(S.(fieldName));

            if isscalar(value) && ...
                    isfinite(value) && ...
                    value > 0

                fps = value;
                return
            end
        end
    end
end

function F = makeTimeByNeuronMatrix(X,fieldName)
% Convert matrix to time x neuron.
%
% Typical calcium files are ROI x time, with considerably fewer ROIs than
% time samples. In that case the matrix is transposed.

    X = double(X);
    X(~isfinite(X)) = 0;

    if ndims(X) ~= 2
        error('Field "%s" must be a 2D matrix.',fieldName);
    end

    nRows = size(X,1);
    nCols = size(X,2);

    if nRows < nCols

        F = X.';

        fprintf('%s: transposed ROI x time -> time x ROI\n', ...
            fieldName);

    else

        F = X;

        fprintf('%s: kept as time x ROI\n', ...
            fieldName);
    end
end

function [assemblies,engagement,engSignificance,info] = ...
    pcaPromaxEnsembles(F,F_forEngagement,opts)

    [nTime,nNeurons] = size(F);

    if isfield(opts,'minShiftFraction')
        minShiftFraction = opts.minShiftFraction;
    else
        minShiftFraction = 0;
    end

    fprintf('\nRunning PCA-promax on:\n');
    fprintf('  Time samples    : %d\n',nTime);
    fprintf('  Candidate cells : %d\n',nNeurons);

    %% ---------------- Z-score --------------------------------------

    Fz = localZscore(F);
    Fz(~isfinite(Fz)) = 0;

    %% ---------------- PCA ------------------------------------------

    [allPCs,eigenvalues] = ...
        localPCAFromCov(Fz);

    %% ---------------- Null model -----------------------------------

    fprintf('Circular-shift iterations: %d\n', ...
        opts.nNullIter);

    [lambdaThreshold,nullInfo] = ...
        maxNullEigenvalueCircularShift( ...
            Fz, ...
            opts);

   if strcmpi(opts.nullEigenvalueMode,'parallelAnalysis') && ...
            isfield(nullInfo,'thresholdByRank')

        % Rank-by-rank comparison. Note this keeps every rank that beats its
        % own null, rather than stopping at the first failure. On simulated
        % data the "stop at first failure" rule collapsed to zero components
        % whenever lambda_1 happened to fall below its own null, which is
        % common with strongly autocorrelated calcium traces.
        significantPCs = ...
            find(eigenvalues > nullInfo.thresholdByRank(:));

    else

        significantPCs = ...
            find(eigenvalues > lambdaThreshold);
    end

    fprintf('Null eigenvalue threshold: %.4f\n', ...
        lambdaThreshold);

    fprintf('Significant PCs: %d\n', ...
        numel(significantPCs));

    %% ---------------- No significant components --------------------

    if isempty(significantPCs)

        assemblies = zeros(nNeurons,0);
        engagement = zeros(0,nTime);
        engSignificance = ones(0,nTime);

        info = struct();

        info.eigenvals = eigenvalues;
        info.lambda_TH = lambdaThreshold;
        info.significantPCs = significantPCs;
        info.nullInfo = nullInfo;

        info.PCs_rot = zeros(nNeurons,0);
        info.neuronLoading = zeros(nNeurons,1);

        info.load_TH = opts.loadingThreshold;

        info.nROIs = nNeurons;
        info.nTime = nTime;

        return
    end

    %% ---------------- Promax ---------------------------------------

    PCs = allPCs(:,significantPCs);

    if size(PCs,2) == 1

        rotatedPCs = PCs;

    else

        rotatedPCs = localPromaxRotation( ...
            PCs, ...
            opts.promaxPower);
    end

    rotatedPCs = ...
        normalizeAndFlip(rotatedPCs);

    %% ---------------- Loading threshold ----------------------------

    zRotatedPCs = ...
        localZscore(rotatedPCs);

    neuronLoading = ...
        max(zRotatedPCs,[],2);

    loadingThreshold = ...
        opts.loadingThreshold;

    fprintf('Loading threshold: %.3f\n', ...
        loadingThreshold);

    %% ---------------- Define assemblies ----------------------------

    [assemblies,rotatedPCs] = ...
        defineAssemblies( ...
            rotatedPCs, ...
            loadingThreshold, ...
            opts.minAssemblySize);

    fprintf('Assemblies before merging: %d\n', ...
        size(assemblies,2));

    %% ---------------- Merge similar assemblies ---------------------

    if opts.mergeSimilarAssemblies && ...
            size(assemblies,2) > 1

        [assemblies,rotatedPCs] = ...
            mergeSimilarAssemblies( ...
                assemblies, ...
                rotatedPCs, ...
                loadingThreshold, ...
                opts.similarityThreshold, ...
                opts.minAssemblySize);
    end

  fprintf('Assemblies after merging: %d\n', ...
        size(assemblies,2));

    %% ---------------- Prune weak assemblies ------------------------

    pruneInfo = struct('applied',false);

    if isfield(opts,'pruneWeakAssemblies') && ...
            opts.pruneWeakAssemblies && ...
            size(assemblies,2) > 0

        [assemblies,rotatedPCs,pruneInfo] = ...
            pruneWeakAssemblies( ...
                assemblies, ...
                rotatedPCs, ...
                F_forEngagement, ...
                opts);

        fprintf('Assemblies after pruning: %d (%d removed)\n', ...
            size(assemblies,2),pruneInfo.nRemoved);
    end

    %% ---------------- Engagement -----------------------------------

    %% ---------------- Engagement -----------------------------------

    [engagement,engSignificance] = ...
        assemblyEngagement( ...
            assemblies, ...
            F_forEngagement);

    %% ---------------- Output information ---------------------------

    info = struct();

    info.eigenvals = eigenvalues;
    info.lambda_TH = lambdaThreshold;
    info.significantPCs = significantPCs;
    info.nullInfo = nullInfo;

    info.PCs_rot = rotatedPCs;
    info.neuronLoading = neuronLoading;

    info.load_TH = loadingThreshold;

    info.pruneInfo = pruneInfo;

    info.nROIs = nNeurons;
    info.nTime = nTime;
end

function [coeff,eigenvalues] = localPCAFromCov(Fz)

    nTime = size(Fz,1);

    if nTime < 2
        error('At least two time samples are required.');
    end

    covarianceMatrix = ...
        (Fz.' * Fz) ./ (nTime - 1);

    covarianceMatrix = ...
        (covarianceMatrix + covarianceMatrix.') ./ 2;

    [eigenvectors,eigenvalueVector] = ...
        eig(covarianceMatrix,'vector');

    [eigenvalues,order] = ...
        sort(real(eigenvalueVector),'descend');

    eigenvalues = ...
        max(eigenvalues,0);

    coeff = ...
        real(eigenvectors(:,order));
end

function [lambdaThreshold,nullInfo] = ...
    maxNullEigenvalueCircularShift(Fz,opts)
% Circular-shift null model for PCA eigenvalues.
%
% Returns:
%   lambdaThreshold   scalar threshold, kept for backward compatibility
%                     (diagnostic figure, info.lambda_TH)
%   nullInfo          includes thresholdByRank, the per-rank thresholds
%                     required by Horn parallel analysis
%
% Self-contained: every option it needs has a local fallback, so it does
% not depend on any edit made elsewhere in the script.

    [nTime,nNeurons] = size(Fz);

    % ---------------- Option fallbacks -------------------------------

    if isfield(opts,'minShiftFraction')
        minShiftFraction = opts.minShiftFraction;
    else
        minShiftFraction = 0;
    end

    if isfield(opts,'nullEigenvalueMode')
        nullEigenvalueMode = opts.nullEigenvalueMode;
    else
        nullEigenvalueMode = 'maxPerIteration';
    end

    % Minimum circular shift. Very small shifts barely decorrelate slow
    % calcium transients and inflate the null.
    minimumShift = max(1,round(minShiftFraction*nTime));
    maximumShift = max(minimumShift,nTime-minimumShift);

    % ---------------- Null iterations --------------------------------

    maxEigenvalueByIteration = ...
        zeros(opts.nNullIter,1);

    pooledEigenvalues = ...
        zeros(opts.nNullIter*nNeurons,1);

    shiftedMatrix = ...
        zeros(nTime,nNeurons);

    poolPosition = 1;

    for iteration = 1:opts.nNullIter

        for neuron = 1:nNeurons

            if maximumShift > minimumShift
                shiftAmount = randi([minimumShift maximumShift]);
            else
                shiftAmount = minimumShift;
            end

            shiftedMatrix(:,neuron) = ...
                circshift( ...
                    Fz(:,neuron), ...
                    shiftAmount);
        end

        [~,currentEigenvalues] = ...
            localPCAFromCov(shiftedMatrix);

        maxEigenvalueByIteration(iteration) = ...
            max(currentEigenvalues);

        currentIndices = ...
            poolPosition:( ...
            poolPosition + numel(currentEigenvalues) - 1);

        pooledEigenvalues(currentIndices) = ...
            currentEigenvalues;

        poolPosition = ...
            currentIndices(end) + 1;
    end

    pooledEigenvalues = ...
        pooledEigenvalues(1:poolPosition-1);

    % ---------------- Per-rank thresholds ----------------------------
    %
    % Each iteration wrote one contiguous block of eigenvalues already
    % sorted in descending order, so reshaping gives a rank x iteration
    % matrix. Note that maxEigenvalueByIteration is exactly row 1 of this
    % matrix, so the old scalar rule and parallel analysis agree on the
    % first component and differ only from rank 2 upwards.

    nullByRank = reshape( ...
        pooledEigenvalues, ...
        nNeurons, ...
        opts.nNullIter);

    thresholdByRank = ...
        prctile(nullByRank,opts.nullPercentile,2);

    thresholdByRank = thresholdByRank(:);

    % ---------------- Scalar threshold -------------------------------

    switch lower(nullEigenvalueMode)

        case lower('maxPerIteration')

            lambdaThreshold = prctile( ...
                maxEigenvalueByIteration, ...
                opts.nullPercentile);

        case lower('parallelAnalysis')

            lambdaThreshold = thresholdByRank(1);

        case lower('pooledEigenvalues')

            lambdaThreshold = prctile( ...
                pooledEigenvalues, ...
                opts.nullPercentile);

        otherwise

            error('Unknown nullEigenvalueMode: %s', ...
                nullEigenvalueMode);
    end

    % ---------------- Output -----------------------------------------

    nullInfo = struct();

    nullInfo.maxEigByIter = ...
        maxEigenvalueByIteration;

    nullInfo.pooledEigenvalues = ...
        pooledEigenvalues;

    nullInfo.thresholdByRank = ...
        thresholdByRank;

    nullInfo.mode = ...
        nullEigenvalueMode;

    nullInfo.percentile = ...
        opts.nullPercentile;

    nullInfo.minShiftFraction = ...
        minShiftFraction;
end

function rotatedPCs = ...
    localPromaxRotation(PCs,powerValue)

    if exist('rotatefactors','file') == 2

        try

            rotatedPCs = rotatefactors( ...
                PCs, ...
                'Method','promax', ...
                'Power',powerValue, ...
                'Normalize','off');

            return

        catch
        end
    end

    [varimaxLoadings,rotationMatrix] = ...
        localVarimax(PCs);

    targetLoadings = ...
        sign(varimaxLoadings) .* ...
        abs(varimaxLoadings).^powerValue;

    obliqueTransform = ...
        varimaxLoadings \ targetLoadings;

    columnNorms = ...
        sqrt(sum(obliqueTransform.^2,1));

    columnNorms(columnNorms == 0) = 1;

    obliqueTransform = ...
        obliqueTransform ./ columnNorms;

    rotatedPCs = ...
        PCs * rotationMatrix * obliqueTransform;
end

function [rotatedLoadings,rotationMatrix] = ...
    localVarimax(loadings)

    [nRows,nColumns] = size(loadings);

    rotationMatrix = eye(nColumns);

    previousObjective = 0;

    maxIterations = 500;
    tolerance = 1e-6;

    for iteration = 1:maxIterations

        currentLoadings = ...
            loadings * rotationMatrix;

        target = ...
            currentLoadings.^3 - ...
            (1/nRows) .* ...
            currentLoadings * ...
            diag(sum(currentLoadings.^2,1));

        [U,S,V] = ...
            svd(loadings.' * target,'econ');

        rotationMatrix = U * V.';

        objective = sum(diag(S));

        if iteration > 1 && ...
                objective < previousObjective*(1+tolerance)

            break
        end

        previousObjective = objective;
    end

    rotatedLoadings = ...
        loadings * rotationMatrix;
end

function rotatedPCs = ...
    normalizeAndFlip(rotatedPCs)

    for component = 1:size(rotatedPCs,2)

        componentNorm = ...
            norm(rotatedPCs(:,component));

        if componentNorm > 0

            rotatedPCs(:,component) = ...
                rotatedPCs(:,component) ./ componentNorm;
        end

        if max(rotatedPCs(:,component)) <= ...
                abs(min(rotatedPCs(:,component)))

            rotatedPCs(:,component) = ...
                -rotatedPCs(:,component);
        end
    end
end

function [assemblies,keptPCs] = ...
    defineAssemblies( ...
        rotatedPCs, ...
        loadingThreshold, ...
        minAssemblySize)

    zRotatedPCs = ...
        localZscore(rotatedPCs);

    zRotatedPCs(~isfinite(zRotatedPCs)) = 0;

    nNeurons = ...
        size(rotatedPCs,1);

    nComponents = ...
        size(rotatedPCs,2);

    initialAssemblies = ...
        zeros(nNeurons,nComponents);

    keepComponent = ...
        false(1,nComponents);

    for component = 1:nComponents

        selectedNeurons = find( ...
            zRotatedPCs(:,component) >= ...
            loadingThreshold);

        if numel(selectedNeurons) >= ...
                minAssemblySize

            initialAssemblies( ...
                selectedNeurons,component) = ...
                rotatedPCs( ...
                    selectedNeurons,component);

            assemblyNorm = ...
                norm(initialAssemblies(:,component));

            if assemblyNorm > 0

                initialAssemblies(:,component) = ...
                    initialAssemblies(:,component) ./ ...
                    assemblyNorm;

                keepComponent(component) = true;
            end
        end
    end

    assemblies = ...
        initialAssemblies(:,keepComponent);

    keptPCs = ...
        rotatedPCs(:,keepComponent);
end

function [assemblies,components] = ...
    mergeSimilarAssemblies( ...
        assemblies, ...
        components, ...
        loadingThreshold, ...
        similarityThreshold, ...
        minAssemblySize)

    while true

        nAssemblies = ...
            size(assemblies,2);

        if nAssemblies < 2
            break
        end

        similarityMatrix = ...
            assemblies.' * assemblies;

        similarityMatrix( ...
            triu(true(nAssemblies),0)) = 0;

        [rowIndices,columnIndices] = ...
            find( ...
                similarityMatrix > ...
                similarityThreshold);

        if isempty(rowIndices)
            break
        end

        similarityValues = ...
            similarityMatrix( ...
                sub2ind( ...
                    size(similarityMatrix), ...
                    rowIndices, ...
                    columnIndices));

        [~,maximumIndex] = ...
            max(similarityValues);

        firstAssembly = ...
            rowIndices(maximumIndex);

        secondAssembly = ...
            columnIndices(maximumIndex);

        mergedComponent = ...
            components(:,firstAssembly) + ...
            components(:,secondAssembly);

        if norm(mergedComponent) == 0

            if norm(components(:,firstAssembly)) >= ...
                    norm(components(:,secondAssembly))

                mergedComponent = ...
                    components(:,firstAssembly);

            else

                mergedComponent = ...
                    components(:,secondAssembly);
            end
        end

        mergedComponent = ...
            mergedComponent ./ ...
            max(norm(mergedComponent),eps);

        indicesToRemove = ...
            sort( ...
                [firstAssembly secondAssembly], ...
                'descend');

        components(:,indicesToRemove) = [];
        assemblies(:,indicesToRemove) = [];

        zMergedComponent = ...
            localZscore(mergedComponent);

        mergedAssembly = ...
            mergedComponent;

        mergedAssembly( ...
            zMergedComponent < loadingThreshold) = 0;

        if nnz(mergedAssembly) >= minAssemblySize && ...
                norm(mergedAssembly) > 0

            mergedAssembly = ...
                mergedAssembly ./ norm(mergedAssembly);

            components = ...
                [components mergedComponent]; %#ok<AGROW>

            assemblies = ...
                [assemblies mergedAssembly]; %#ok<AGROW>
        end
    end
end

function [engagement,significance] = ...
    assemblyEngagement(assemblies,activityMatrix)

    binaryRaster = ...
        activityMatrix > 0;

    [nTime,nNeurons] = ...
        size(binaryRaster);

    nAssemblies = ...
        size(assemblies,2);

    engagement = ...
        zeros(nAssemblies,nTime);

    significance = ...
        ones(nAssemblies,nTime);

    populationActivity = ...
        sum(binaryRaster,2);

    for assembly = 1:nAssemblies

        members = find( ...
            assemblies(:,assembly) ~= 0);

        currentAssemblySize = ...
            numel(members);

        if currentAssemblySize == 0
            continue
        end

        intersection = ...
            sum(binaryRaster(:,members),2);

        denominator = ...
            currentAssemblySize + ...
            populationActivity;

        validTimes = ...
            denominator > 0;

        engagement(assembly,validTimes) = ...
            ( ...
            2 .* intersection(validTimes) ./ ...
            denominator(validTimes) ...
            ).';

        activeTimes = ...
            find(intersection > 0);

        for iTime = 1:numel(activeTimes)

            t = activeTimes(iTime);

            significance(assembly,t) = ...
                localHypergeomUpperTail( ...
                    nNeurons, ...
                    currentAssemblySize, ...
                    populationActivity(t), ...
                    intersection(t));
        end
    end
end

function pUpper = ...
    localHypergeomUpperTail( ...
        populationSize, ...
        numberSuccesses, ...
        numberDraws, ...
        observedSuccesses)

    maximumPossible = ...
        min(numberSuccesses,numberDraws);

    minimumPossible = ...
        max( ...
            0, ...
            numberDraws - ...
            (populationSize-numberSuccesses));

    if observedSuccesses < minimumPossible

        pUpper = 1;
        return
    end

    if observedSuccesses > maximumPossible

        pUpper = 0;
        return
    end

    possibleValues = ...
        observedSuccesses:maximumPossible;

    logProbability = ...
        localLogChoose( ...
            numberSuccesses, ...
            possibleValues) + ...
        localLogChoose( ...
            populationSize-numberSuccesses, ...
            numberDraws-possibleValues) - ...
        localLogChoose( ...
            populationSize, ...
            numberDraws);

    pUpper = ...
        sum(exp(logProbability));

    pUpper = ...
        min(max(pUpper,0),1);
end

function output = localLogChoose(n,k)

    output = ...
        gammaln(n+1) - ...
        gammaln(k+1) - ...
        gammaln(n-k+1);
end

function Z = localZscore(X)

    columnMean = ...
        mean(X,1,'omitnan');

    columnStd = ...
        std(X,0,1,'omitnan');

    columnStd( ...
        columnStd == 0 | ...
        ~isfinite(columnStd)) = 1;

    Z = ...
        (X-columnMean) ./ columnStd;

    Z(~isfinite(Z)) = 0;
end

function figH = ...
    plotPcaPromaxDiagnostics( ...
        info, ...
        assemblies, ...
        assemblySize, ...
        sessionName)

    figH = figure( ...
        'Color','w', ...
        'Name',['PCA-promax diagnostics: ' sessionName], ...
        'Position',[100 100 1500 450]);

    tiledlayout( ...
        figH, ...
        1,4, ...
        'TileSpacing','compact', ...
        'Padding','compact');

    %% Eigenvalues

    nexttile;

    plot( ...
        info.eigenvals, ...
        'k.-', ...
        'LineWidth',1);

    hold on;

    yline( ...
        info.lambda_TH, ...
        'r--', ...
        'LineWidth',1.5);

    xlabel('PC index');
    ylabel('Eigenvalue');
    title('PCA and null threshold');
    box off;

    %% Loading distribution

    nexttile;

    histogram( ...
        info.neuronLoading, ...
        30);

    hold on;

    xline( ...
        info.load_TH, ...
        'r--', ...
        'LineWidth',1.5);

    xlabel('Maximum z-loading');
    ylabel('Candidate cells');

    title(sprintf( ...
        'Loading threshold = %.2f', ...
        info.load_TH));

    box off;

    %% Membership matrix

    nexttile;

    if isempty(assemblies)

        imagesc( ...
            zeros(info.nROIs,1));

    else

        imagesc( ...
            assemblies ~= 0);
    end

    xlabel('Assembly');
    ylabel('Candidate cell');

    title(sprintf( ...
        'Membership, A = %d', ...
        size(assemblies,2)));

    colormap(gca,gray);
    box off;

    %% Assembly sizes

    nexttile;

    if isempty(assemblySize)

        axis off;

        text( ...
            0.5,0.5, ...
            'No assemblies detected', ...
            'HorizontalAlignment','center');

    else

        bar(assemblySize);

        xlabel('Assembly');
        ylabel('Candidate cells');
        title('Assembly size');
        box off;
    end

    sgtitle( ...
        sessionName, ...
        'Interpreter','none');
end

function showRoiImageForAssemblies(ax, roiImage)
% Display hdMetrics.roiImage while handling common image formats.
%
% Supported:
%   height x width
%   height x width x 3
%   height x width x nPlanes
%
% For a multi-plane non-RGB image, a maximum-intensity projection is used.

    imageToShow = roiImage;

    imageToShow = double(imageToShow);

    imageToShow(~isfinite(imageToShow)) = 0;

    %% ---------------- Determine image format ------------------------

    if ndims(imageToShow) == 2

        % Already a grayscale 2D image

    elseif ndims(imageToShow) == 3 && ...
            size(imageToShow,3) == 3

        % RGB image: preserve the three channels

    elseif ndims(imageToShow) == 3

        % Multi-plane grayscale image
        imageToShow = max(imageToShow,[],3);

    else

        % Squeeze singleton dimensions
        imageToShow = squeeze(imageToShow);

        if ndims(imageToShow) > 2
            imageToShow = max(imageToShow,[],3);
        end
    end

    %% ---------------- Robust intensity normalization ----------------

    if ndims(imageToShow) == 2

        finiteValues = imageToShow(isfinite(imageToShow));

        if isempty(finiteValues)

            imageToShow = zeros(size(imageToShow));

        else

            lowerLimit = prctile(finiteValues,1);
            upperLimit = prctile(finiteValues,99.5);

            if upperLimit <= lowerLimit
                lowerLimit = min(finiteValues);
                upperLimit = max(finiteValues);
            end

            if upperLimit > lowerLimit

                imageToShow = ...
                    (imageToShow-lowerLimit) ./ ...
                    (upperLimit-lowerLimit);

            else

                imageToShow = zeros(size(imageToShow));
            end

            imageToShow = min(max(imageToShow,0),1);
        end

        imagesc(ax,imageToShow);
        colormap(ax,gray);

    else

        % Normalize each RGB channel independently
        for channel = 1:3

            currentChannel = imageToShow(:,:,channel);
            finiteValues = currentChannel(isfinite(currentChannel));

            if isempty(finiteValues)
                imageToShow(:,:,channel) = 0;
                continue
            end

            lowerLimit = prctile(finiteValues,1);
            upperLimit = prctile(finiteValues,99.5);

            if upperLimit > lowerLimit

                currentChannel = ...
                    (currentChannel-lowerLimit) ./ ...
                    (upperLimit-lowerLimit);

            else

                currentChannel = zeros(size(currentChannel));
            end

            imageToShow(:,:,channel) = ...
                min(max(currentChannel,0),1);
        end

        image(ax,imageToShow);
    end

    set(ax,'YDir','reverse');

    axis(ax,'image');
    axis(ax,'off');
end

function [assemblies,components,pruneInfo] = ...
    pruneWeakAssemblies(assemblies,components,activityMatrix,opts)
% Remove assemblies that are no more coherent than random cell sets.
%
% Two statistics per assembly, each with its own null:
%
%   synchrony   max over bins of (member cells active in the bin) divided by
%               assembly size. A group that never fires together scores low
%               even if its cells are individually busy.
%
%   coherence   mean pairwise correlation between member cells in the binary
%               raster.
%
% DEVIATION FROM ROMANO: the null here draws random cell sets of exactly the
% same size as the assembly being tested, separately for each assembly.
% Romano draws one random assembly index per iteration, so a single null
% distribution is shared across assemblies of different sizes. Both
% statistics depend on assembly size, so a shared null penalises large
% assemblies and favours small ones. Size matching removes that bias.

    binaryRaster = activityMatrix > 0;

    [nTime,nCells] = size(binaryRaster);

    nAssemblies = size(assemblies,2);

    pruneInfo = struct( ...
        'applied',true, ...
        'nBefore',nAssemblies, ...
        'nRemoved',0, ...
        'synchrony',nan(nAssemblies,1), ...
        'coherence',nan(nAssemblies,1), ...
        'synchronyThreshold',nan(nAssemblies,1), ...
        'coherenceThreshold',nan(nAssemblies,1), ...
        'removed',false(nAssemblies,1));

    if nAssemblies == 0 || nTime < 2
        return
    end

    binFrames = max(1,round(opts.pruneSyncBinFrames));
    binEdges = 1:binFrames:nTime;

    % ---------------- Observed statistics ----------------------------

    memberSets = cell(nAssemblies,1);

    for a = 1:nAssemblies

        members = find(assemblies(:,a) ~= 0);
        memberSets{a} = members;

        if numel(members) < 2
            continue
        end

        pruneInfo.synchrony(a) = ...
            maxSyncFraction(binaryRaster,members,binEdges);

        pruneInfo.coherence(a) = ...
            meanPairwiseCorrelation(binaryRaster,members);
    end

    % ---------------- Size-matched nulls -----------------------------

    for a = 1:nAssemblies

        nMembers = numel(memberSets{a});

        if nMembers < 2 || ~isfinite(pruneInfo.synchrony(a))
            continue
        end

        nullSync = nan(opts.pruneNullIterations,1);
        nullCoherence = nan(opts.pruneNullIterations,1);

        for iNull = 1:opts.pruneNullIterations

            pick = randperm(nCells,nMembers);

            nullSync(iNull) = ...
                maxSyncFraction(binaryRaster,pick,binEdges);

            nullCoherence(iNull) = ...
                meanPairwiseCorrelation(binaryRaster,pick);
        end

        pruneInfo.synchronyThreshold(a) = max( ...
            prctile(nullSync,95), ...
            opts.pruneMinSyncFraction);

        pruneInfo.coherenceThreshold(a) = ...
            prctile(nullCoherence,95);
    end

    % ---------------- Decide ------------------------------------------

    failsSync = pruneInfo.synchrony < pruneInfo.synchronyThreshold;
    failsCoherence = pruneInfo.coherence < pruneInfo.coherenceThreshold;

    % Assemblies too small to score at all are removed.
    unscored = ~isfinite(pruneInfo.synchrony) | ...
        ~isfinite(pruneInfo.coherence);

    weak = failsSync | failsCoherence | unscored;

    if isfield(opts,'pruneRelativeZ') && opts.pruneRelativeZ && ...
            nAssemblies >= 3

        zSync = zscoreSafe(pruneInfo.synchrony);
        zCoherence = zscoreSafe(pruneInfo.coherence);

        combined = zscoreSafe((zSync+zCoherence)/2);

        weak = weak | combined < -1;
    end

    pruneInfo.removed = weak;
    pruneInfo.nRemoved = sum(weak);

    assemblies(:,weak) = [];
    components(:,weak) = [];
end

function value = maxSyncFraction(binaryRaster,members,binEdges)
% Largest fraction of the member set active within one bin.

    value = NaN;

    if isempty(members)
        return
    end

    nBins = numel(binEdges)-1;

    if nBins < 1
        return
    end

    counts = zeros(nBins,1);

    for iBin = 1:nBins

        window = binEdges(iBin):(binEdges(iBin+1)-1);

        counts(iBin) = ...
            sum(any(binaryRaster(window,members),1));
    end

    value = max(counts)/numel(members);
end

function value = meanPairwiseCorrelation(binaryRaster,members)
% Mean off-diagonal correlation among the member cells.

    value = NaN;

    if numel(members) < 2
        return
    end

    block = double(binaryRaster(:,members));

    keep = std(block,0,1) > 0;

    if sum(keep) < 2
        return
    end

    C = corr(block(:,keep));

    upper = triu(true(sum(keep)),1);

    value = mean(C(upper),'omitnan');
end

function z = zscoreSafe(x)

    x = x(:);

    mu = mean(x,'omitnan');
    sigma = std(x,'omitnan');

    if ~isfinite(sigma) || sigma == 0
        z = zeros(size(x));
    else
        z = (x-mu)/sigma;
    end
end