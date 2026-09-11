function results = Assem_detection(varargin)
%ASSEM_DETECTION Non-interactive batch PCA-Promax assembly detection.
% Scans /home/blanche/data/{surface,molino,pachon}/rec* and writes only to
% data_processed/assemblies and outputs/assemblies in this project.
rootDir=fileparts(mfilename('fullpath'));
opts=assembly_default_options(rootDir);
lockedPath=fullfile(opts.processedRoot,'validation','selected_configuration.json');
if isfile(lockedPath)
 locked=jsondecode(fileread(lockedPath)); names={'loadingThreshold','minAssemblySize','nullPercentile'};
 for i=1:numel(names),if isfield(locked,names{i}),opts.(names{i})=locked.(names{i});end,end
end
opts=override(opts,varargin{:});
fprintf("\n=== Assembly detection started %s ===\n",char(datetime("now","Format","yyyy-MM-dd HH:mm:ss")));
fprintf("Output: %s\n",fullfile(opts.processedRoot,opts.analysisLabel)); drawnow;
results=run_assembly_batch(opts);
fprintf("=== Assembly detection finished %s ===\n\n",char(datetime("now","Format","yyyy-MM-dd HH:mm:ss"))); drawnow;
end
function o=override(o,varargin),if mod(numel(varargin),2),error('Options must be name/value pairs.');end,for i=1:2:numel(varargin),if ~isfield(o,varargin{i}),error('Unknown option: %s',varargin{i});end,o.(varargin{i})=varargin{i+1};end,end
