# Assembly detection and analysis

This folder contains the PCA–Promax assembly-detection pipeline and the analyses that compare detected assemblies across surface, Molino, and Pachón fish. Run MATLAB from the repository or add this folder to the MATLAB path:

```matlab
addpath('assembly_analysis_hd')
```

## 1. Detect assemblies

The main entry point is `Assem_detection.m`:

```matlab
results = Assem_detection();
```

It finds eligible recordings under `/home/blanche/data/{surface,molino,pachon}/rec*`, applies `QC_ENTRY_KEPT_files.csv`, uses candidate HD cells, and analyzes the first 30 minutes by default. Results and the session manifest are written under `assembly_analysis_hd/data_processed/assemblies`; figures are written under `assembly_analysis_hd/outputs/assemblies`.

Detection defaults are defined in `assembly_default_options.m`. The chosen parameter set was chosen to be:
- `nNullIter`:1000
- `nullPercentile`:95
- `loadingThreshold`: 1.75
- `minAssemblySize`:2
- `similarityThreshold`:0.60
- `randomSeed`:1

The principal parameters to review are:

- Data and selection: `dataRoot`, `morphs`, `analysisDurationMin`, and `maxSessionsPerMorph`.
- PCA/null detection: `nNullIter`, `nullPercentile`, `minShiftFraction`, and `promaxPower`.
- Membership: `loadingThreshold` and `minAssemblySize`.
- Post-processing: `mergeSimilarAssemblies`, `similarityThreshold`, `pruneWeakAssemblies`, and the `prune*` options.
- Output/reproducibility: `analysisLabel`, `randomSeed`, `makeFigures`, and `overwriteExisting`.

`analysisLabel` should be unique when changing detection settings. `Assem_detection.m` also reads a legacy `data_processed/assemblies/validation/selected_configuration.json` if that file still exists; explicit name/value arguments take precedence over it.

## 2. Analyze detected assemblies

Pass the returned `results` structure to each analysis. Most analysis functions can also be called with `[]`, in which case they load the saved `session_manifest.mat`.

```matlab
properties = compare_assembly_properties_across_morphs(results);
triggered  = compare_assembly_triggered_averages_across_morphs(results);
spatial    = compare_assembly_spatial_compactness_across_morphs(results);
corrs      = compare_assembly_correlations_across_morphs(results);
```

Each function's adjustable parameters are the `addParameter` entries near the top of its file and can be supplied as name/value pairs.

- `compare_assembly_properties_across_morphs.m`: assembly size, activity, episode, phase, behavioral-class, compactness, CCG, and transition summaries. Tune activation options, `nCompactnessPermutations`, `ccgMaxLagFrames`, `transitionLagFrames`.
- `compare_assembly_triggered_averages_across_morphs.m`: onset-aligned assembly/member activity and pseudo-onset controls. Tune `relativeFrames`, `isolatedWindows`, `controlWindowCount`, and `controlRandomSeed`.
- `compare_assembly_spatial_compactness_across_morphs.m`: Hansen-style anatomical compactness. Tune `distanceBinWidthUm`, `maxDistanceUm`, `sensitivityBinWidthsUm`, `compactEntropyCutoff`, and bootstrap counts.
- `compare_assembly_correlations_across_morphs.m`: correlations as a function of preferred-phase distance, onset-shuffle controls, and turn-related analyses. Tune `phaseDistanceThresholdDeg`, `phaseDistanceBinWidthDeg`, `maxLagSeconds`, `turnWindowSeconds`, shuffle/bootstrap counts, and behavior-time settings.

The properties, triggered-average, and correlation analyses share activation defaults from `assembly_activation_defaults.m`. If activation criteria are overridden, pass the same values to all three analyses. The shared options are `minActiveMembers`, `onRecruitmentFraction`, `offRecruitmentFraction`, `minMatchingIndex`, `activationNullAlpha`, `nNullShuffles`, `minDurationSeconds`, `minNullShiftSeconds`, and `activationRandomSeed`.

Finally, `summarize_assembly_detection_controls.m` reads completed comparison outputs and creates a compact control summary. It is a script rather than a function; edit the `folders` and `metrics` arrays near its top before running it:

```matlab
run('1_Assemblies/summarize_assembly_detection_controls.m')
```
