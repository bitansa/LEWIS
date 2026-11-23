% =============================================================================
% run_pwling_batch.m   (Pairwise-LiNGAM, stage 2 of 3)
%
% Runs the MATLAB function pwling of Hyvarinen and Smith (2013) on every
% dataset listed in manifest.csv, and writes the likelihood-ratio matrices.
% For each dataset it calls
%
%     LR = pwling(X', 1)
%
% pwling takes the variables in the rows, so the n x p dataset is
% transposed. Method 1 is the general likelihood ratio of the paper, which
% does not assume a particular error distribution. LR(i, j) > 0 means the
% direction i -> j. Stage 3 orients the edges with this rule.
%
% Stage 1 (R) writes the datasets to inputs/ and the list to manifest.csv:
%     01_motivation/pairwise_lingam/export_motivation_for_matlab.R
%     02_comparison/pairwise_lingam/export_datasets_for_matlab.R
% Stage 3 (R) reads outputs/ and evaluates the estimates:
%     01_motivation/pairwise_lingam/import_motivation_pairwise.R
%     02_comparison/pairwise_lingam/import_pairwise_results.R
%
% THE CODE pwling IS NOT PART OF THIS REPOSITORY. It is the code of
%     A. Hyvarinen and S. M. Smith (2013). Pairwise likelihood ratios for
%     estimation of non-Gaussian structural equation models. Journal of
%     Machine Learning Research, 14, 111-152.
% Download page: https://www.cs.helsinki.fi/u/ahyvarin/code/pwcausal/
% File:          https://www.cs.helsinki.fi/u/ahyvarin/code/pwcausal/pwling.zip
% Unpack the file into a folder of your choice. Then set the environment
% variable PWLING_PATH to the folder that contains pwling.m before you
% start MATLAB, for example
%     export PWLING_PATH=/path/to/pwling
% You can also set a MATLAB variable PWLING_PATH before you call this
% script.
%
% HOW TO RUN. Make the folder that holds manifest.csv the current folder,
% put this folder on the MATLAB path, and call the script by name. From a
% terminal, in the folder of the study (for example
% reproducibility/02_comparison/pairwise_lingam):
%
%     matlab -batch "addpath('../../common/pairwise_lingam'), run_pwling_batch"
%
% or, in an open MATLAB session:
%
%     cd <reproducibility>/02_comparison/pairwise_lingam
%     addpath('<reproducibility>/common/pairwise_lingam')
%     run_pwling_batch
%
% Do not call it as "run path/to/run_pwling_batch.m". The command run
% changes the current folder to the folder of the script, and manifest.csv
% is then not found. Instead, you can set a MATLAB variable DATA_DIR to the
% folder that holds manifest.csv before you call the script.
%
% PARALLEL RUN. With the Parallel Computing Toolbox, the datasets run on
% NUM_WORKERS workers (default: all cores but one). Set a MATLAB variable
% NUM_WORKERS = 1 before the call for a serial run. Without the toolbox,
% the loop runs serially.
%
% The script can be restarted. It skips every dataset whose output file
% exists already. It stops with an error at the end if a dataset failed.
%
% Columns of manifest.csv:
%   experiment, cell_index, rep_id, graph_type, p, n, rho, effect_scale,
%   X_file, LR_file
% Output files in outputs/:
%   <LR_file>      the p x p matrix LR, comma separated, no header
%   runtimes.csv   experiment, cell_index, rep_id, LR_file, runtime_s, status
%                  (one row per run of a dataset, appended)
% =============================================================================

% ---- Folder with manifest.csv, inputs/ and outputs/ -------------------------
% The folder is DATA_DIR if the user set it, and the current folder
% otherwise. The script does not assign DATA_DIR itself, so a later call in
% the same MATLAB session from another folder uses that folder.
if exist('DATA_DIR', 'var') && ~isempty(DATA_DIR)
    data_dir = DATA_DIR;
else
    data_dir = pwd;
end
fprintf('Data folder: %s\n', data_dir);

% ---- Location of pwling ------------------------------------------------------
% Read PWLING_PATH from the environment unless it is set in MATLAB already.
% Then put the folder on the MATLAB path and check that pwling.m is found.
if ~exist('PWLING_PATH', 'var') || isempty(PWLING_PATH)
    PWLING_PATH = getenv('PWLING_PATH');
end
if isempty(PWLING_PATH)
    error(['Set PWLING_PATH to the folder that contains pwling.m ', ...
           '(environment variable or MATLAB variable).']);
end
addpath(PWLING_PATH);
if exist('pwling', 'file') ~= 2
    error('pwling.m was not found in PWLING_PATH = %s.', PWLING_PATH);
end

% ---- Parallel workers --------------------------------------------------------
% One core is left for the operating system. Set NUM_WORKERS = 1 for a
% serial run.
if ~exist('NUM_WORKERS', 'var') || isempty(NUM_WORKERS)
    NUM_WORKERS = max(1, feature('numcores') - 1);
end
if NUM_WORKERS > 1 && exist('gcp', 'file') == 2
    if isempty(gcp('nocreate'))
        parpool('local', NUM_WORKERS);
    end
    fprintf('Parallel pool with %d workers.\n', NUM_WORKERS);
else
    fprintf('Serial run.\n');
end

% ---- Read the manifest -------------------------------------------------------
% Every column is read with a fixed type, so that the file names stay text.
manifest_path = fullfile(data_dir, 'manifest.csv');
if ~isfile(manifest_path)
    error(['manifest.csv was not found in %s. Make that folder the ', ...
           'current folder and call the script by name, or set DATA_DIR.'], ...
          data_dir);
end
opts = delimitedTextImportOptions('NumVariables', 10, ...
    'Delimiter', ',', ...
    'DataLines', [2, Inf], ...
    'ExtraColumnsRule', 'ignore', ...
    'VariableNames', {'experiment','cell_index','rep_id','graph_type', ...
                      'p','n','rho','effect_scale','X_file','LR_file'}, ...
    'VariableTypes', {'string','double','double','string', ...
                      'double','double','double','double', ...
                      'string','string'});
M = readtable(manifest_path, opts);
fprintf('Manifest with %d datasets.\n', height(M));

% ---- Input and output folders ------------------------------------------------
in_dir = fullfile(data_dir, 'inputs');
out_dir = fullfile(data_dir, 'outputs');
if ~isfolder(in_dir)
    error('The folder inputs/ was not found in %s. Run stage 1 first.', data_dir);
end
% Create the folder of the outputs.
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% ---- Log of the run times ----------------------------------------------------
% LR_file is the key that stage 3 uses to match a run time to a dataset.
runtime_log = fullfile(out_dir, 'runtimes.csv');
if ~isfile(runtime_log)
    fid = fopen(runtime_log, 'w');
    fprintf(fid, 'experiment,cell_index,rep_id,LR_file,runtime_s,status\n');
    fclose(fid);
end

% ---- Columns of the manifest as plain arrays (needed by parfor) --------------
N = height(M);
expts = string(M.experiment);
ks = M.cell_index;
rs = M.rep_id;
LR_fns = string(M.LR_file);
X_fns = string(M.X_file);
status_arr = strings(N, 1);
runtime_arr = nan(N, 1);
errmsg_arr = strings(N, 1);

% ---- Main loop ---------------------------------------------------------------
% Each dataset is read, transposed, passed to pwling, and its LR matrix is
% written. An error is kept with its message and reported after the loop.
t_batch = tic;
parfor row = 1:N
    LR_out = fullfile(out_dir, char(LR_fns(row)));
    % Skip a dataset that is done already.
    if isfile(LR_out)
        status_arr(row) = "DONE";
        continue
    end
    X_fn = fullfile(in_dir, char(X_fns(row)));
    if ~isfile(X_fn)
        status_arr(row) = "MISSING";
        errmsg_arr(row) = "input file not found";
        continue
    end
    try
        % pwling takes the variables in the rows.
        X = csvread(X_fn);
        t0 = tic;
        LR = pwling(X', 1);
        rt = toc(t0);
        csvwrite(LR_out, LR);
        status_arr(row) = "OK";
        runtime_arr(row) = rt;
    catch err
        status_arr(row) = "ERROR";
        errmsg_arr(row) = string(err.message);
    end
end

% ---- Append the run times after the loop -------------------------------------
% A parfor loop cannot write to one file safely, so the log is written here.
fid = fopen(runtime_log, 'a');
for row = 1:N
    st = status_arr(row);
    if st == "OK"
        fprintf(fid, '%s,%d,%d,%s,%.6f,OK\n', ...
                expts(row), ks(row), rs(row), LR_fns(row), runtime_arr(row));
    elseif st == "ERROR" || st == "MISSING"
        fprintf(fid, '%s,%d,%d,%s,NA,%s\n', ...
                expts(row), ks(row), rs(row), LR_fns(row), st);
    end
end
fclose(fid);

% ---- Report ------------------------------------------------------------------
% Print every failed dataset and stop with an error if there is any.
n_done = sum(status_arr == "OK") + sum(status_arr == "DONE");
err_idx = find(status_arr == "ERROR" | status_arr == "MISSING");
for i = 1:numel(err_idx)
    ir = err_idx(i);
    fprintf('  [row %d] %s: %s\n', ir, status_arr(ir), errmsg_arr(ir));
end
fprintf('Done in %.1f s. %d datasets done, %d failed.\n', ...
        toc(t_batch), n_done, numel(err_idx));
if ~isempty(err_idx)
    error('%d datasets failed. Fix the cause and run the script again.', ...
          numel(err_idx));
end
