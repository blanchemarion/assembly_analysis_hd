% Fish-level detection controls. Run this script from any working directory.
% Surface comparisons are separate two-group Kruskal-Wallis p-values,
% without correction. KW_p is the three-morph omnibus.
% SEM = sample SD / sqrt(number of finite fish values); n < 2 gives NaN.
% Reads existing results only; writes one summary CSV in outputs/assemblies.
baseDir = fullfile(fileparts(mfilename('fullpath')), 'outputs', 'assemblies');
folders = [ ...
 "z1p5_m0p75_morph_comparison_raw_candidates_merge080_K2_Ron0p5_Roff0p5_MI0p5_a0p05_n1000_dur0p5_safe10_seed1"
 "z1p75_m0p6_morph_comparison_raw_candidates_thr175_K2_Ron0p5_Roff0p5_MI0p5_a0p05_n1000_dur0p5_safe10_seed1"
 "z1p25_m0p6_morph_comparison_raw_candidates_thr125_K2_Ron0p5_Roff0p5_MI0p5_a0p05_n1000_dur0p5_safe10_seed1"];
metrics = ["meanAssemblySpikeRatePerMinuteOutsideActivationEpisodes", ...
    "medianEpisodeDurationFrames", "medianWithinAssemblyPhaseR"];
groups = ["surface", "molino", "pachon"];
summary = table('Size', [9 8], ...
    'VariableTypes', [repmat({'string'},1,5), repmat({'double'},1,3)], ...
    'VariableNames', {'configuration','metric','Surface_mean_SEM', ...
    'Molino_mean_SEM','Pachon_mean_SEM','KW_p','Surface_Molino_p','Surface_Pachon_p'});
summary{:,6:8} = NaN;
for c = 1:numel(folders)
    folder = fullfile(baseDir, folders(c));
    if ~isfolder(folder) && c == 1 && isfolder(folder + "1")
        warning('Requested seed folder absent; using %s.', folder + "1");
        folders(c) = folders(c) + "1";
        folder = folder + "1";
    end
    fish = readResult(folder, 'fish_level_metrics.csv', 'fishLevel');
    kw = readResult(folder, 'morph_omnibus_kw_tests.csv', 'morphOmnibusKWTests');
    for m = 1:numel(metrics)
        row = (c-1)*numel(metrics)+m;
        metric = metrics(m);
        summary.configuration(row) = folders(c);
        summary.metric(row) = metric;
        summary{row,3:5} = repmat("NaN ± NaN",1,3);
        if ~all(ismember(["morph","fish",metric], string(fish.Properties.VariableNames)))
            warning('Missing metric or fish-level data: %s / %s.', folders(c), metric);
            continue
        end
        g = lower(strtrim(string(fish.morph)));
        g = replace(g, "pachón", "pachon");
        % Refuse repeated fish observations rather than silently overweight them.
        keys = table(g, string(fish.fish));
        if height(unique(keys,'rows')) ~= height(keys)
            warning('Duplicate fish identifiers in %s; skipping %s.', folders(c), metric);
            continue
        end
        y = double(fish.(metric));
        valid = isfinite(y) & ismember(g,groups);
        n = zeros(1,3);
        for j = 1:3
            v = y(valid & g == groups(j));
            n(j) = numel(v);
            if isempty(v)
                warning('No finite %s fish values: %s / %s.', groups(j), folders(c), metric);
                continue
            end
            se = NaN;
            if n(j) > 1, se = std(v,0)/sqrt(n(j)); end
            summary{row,j+2} = string(sprintf('%.6g ± %.6g',mean(v),se));
        end
        if any(n == 0), continue, end
        % Only reuse an omnibus result with df=2 (the saved df=1 tests
        % compare two morphs and cannot serve as the three-morph omnibus).
        if all(ismember({'metric','degreesOfFreedom','pOmnibus'},kw.Properties.VariableNames))
            hit = find(string(kw.metric)==metric & kw.degreesOfFreedom==2,1);
            if ~isempty(hit), summary.KW_p(row) = kw.pOmnibus(hit); end
        end
        if ~isfinite(summary.KW_p(row))
            try
                summary.KW_p(row) = kruskalwallis(y(valid),cellstr(g(valid)),'off');
            catch err
                warning('Cannot compute KW for %s / %s: %s',folders(c),metric,err.message);
            end
        end
        for j = 2:3
            if all(ismember({'metric','group1','group2','degreesOfFreedom','pOmnibus'},kw.Properties.VariableNames))
                hit = find(string(kw.metric)==metric & kw.degreesOfFreedom==1 & ...
                    ((string(kw.group1)=="surface" & string(kw.group2)==groups(j)) | ...
                     (string(kw.group2)=="surface" & string(kw.group1)==groups(j))),1);
                if ~isempty(hit), summary{row,j+5} = kw.pOmnibus(hit); end
            end
            if ~isfinite(summary{row,j+5}) && n(1)>=3 && n(j)>=3
                pair = valid & (g=="surface" | g==groups(j));
                try
                    summary{row,j+5} = kruskalwallis(y(pair),cellstr(g(pair)),'off');
                catch err
                    warning('Cannot compute Surface-%s KW for %s / %s: %s', ...
                        groups(j),folders(c),metric,err.message);
                end
            elseif ~isfinite(summary{row,j+5})
                warning('Insufficient fish for Surface-%s KW: %s / %s.',groups(j),folders(c),metric);
            end
        end
    end
end
outputFile = fullfile(baseDir,'assembly_detection_controls_summary.csv');
writetable(summary,outputFile,'Encoding','UTF-8');
fprintf('Saved %d rows to %s\nSurface comparisons: two-group Kruskal-Wallis (uncorrected).\n',height(summary),outputFile);

function T = readResult(folder, csvName, field)
% Prefer lightweight CSVs; fall back to the matching fish-level report table.
T = table();
csvFile = fullfile(folder,csvName);
try
    if isfile(csvFile)
        T = readtable(csvFile,'TextType','string');
        return
    end
    matFile = fullfile(folder,'assembly_morph_comparison_report.mat');
    if isfile(matFile)
        S = load(matFile,'report');
        if isfield(S,'report') && isfield(S.report,field) && istable(S.report.(field))
            T = S.report.(field);
        end
    end
catch err
    warning('Could not read %s: %s',csvFile,err.message);
end
end
