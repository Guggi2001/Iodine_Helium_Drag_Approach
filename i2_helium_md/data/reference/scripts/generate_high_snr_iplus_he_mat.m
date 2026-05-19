% Generate the averaged I+He high-SNR VMI reference MAT used by the Python
% paper_v2 pipeline.
%
% Replaces the previous personal "generate_high_snr_matfile_data" workflow,
% which accidentally kept only the first file because of a debug-leftover
% line `filenumbers = filenumbers(1);`. Every shipped MAT inspected on
% 2026-05-19 had n=1 (e.g. ressumI2HeNI^+He.mat -> filenumbers = [43702],
% from the 17.10.24 session).
%
% The reference dataset is now aligned to simulation_image.m
% (legacy_matlab_repository/single_pulse_simulation/HeDFT_comparison/
% simulation_image.m:153) and to the paper_cov pipeline:
%
%   FILENUMBERS = [45668, 45662, 45667, 45686]   % I+He droplet, low doping,
%                                                % 300 mW, 03.12.24
%   VM_CENTER   = [524.5297, 380.8430]           % hardcoded center from
%                                                % simulation_image.m:159
%   VELOCITY_FACTOR = 8.6178                     % canonical VMI calibration
%
% The only intended improvement over simulation_image.m is optional global
% background subtraction. Default is OFF (NaN), because the previous global
% bg file 43655 is from the 17.10.24 session and a same-session no-gate
% bg frame from 03.12.24 has not been identified. Set GLOBAL_BG_FILENUMBER
% once a 03.12.24 bg frame is known.
%
% Scope: I+He droplet channel only (set 2). Other channels can be added by
% duplicating the SET block.
%
% Output (in-place replacement, MAT v7):
%   i2_helium_md/old_data/ressumI2HeNI^+He.mat
%       res_sum     : averaged plot_processed_VMI struct (1/N sum over N files)
%       filenumbers : 1xN array of the file IDs that went into res_sum
%
% Diagnostic figure:
%   i2_helium_md/old_data/ressumI2HeNI^+He.fig
%
% Requirements:
%   - setup_VMI_path_office_flir on the MATLAB path
%     (defines path_to_processed_data and puts the VMI toolbox on the path)
%   - process_raw_VMI, plot_processed_VMI, add_processed_data,
%     subtract_processed_data, multiply_processed_data
%   - autocenter_from_extended_data (only used when INSPECT_CENTERS_ONLY = true)
%   - colorcet (only used for the diagnostic figure)
%
% Two modes via the INSPECT_CENTERS_ONLY flag in USER SETTINGS:
%   true  -> report per-file auto-detected centers, plot them, then return.
%            Use this once to pick VM_CENTER.
%   false -> use the hardcoded VM_CENTER, average, save .mat + .fig.

clear; close all;

setup_VMI_path_office_flir;

% ===================== USER SETTINGS =====================
INSPECT_CENTERS_ONLY = false;

% Hardcode after running the inspection pass once.
% This is the simulation_image.m / paper_cov hardcoded center.
VM_CENTER = [524.5297, 380.8430];

% Global no-gate-voltage background subtracted from every measurement.
% Set to NaN to disable (default). 43655 is the 17.10.24 bg frame; do NOT
% use it for the 03.12.24 FILENUMBERS below. Set to a same-session bg
% frame from 03.12.24 once one has been identified.
GLOBAL_BG_FILENUMBER = NaN;

% Velocity calibration. Aligned to the paper_cov pipeline value 8.6178
% (cov exporter line 81) and the export_paper_v2_reference_data.m vf_single.
% The previous personal script had 10.0995, which gave a different absolute
% velocity axis.
VELOCITY_FACTOR = 8.6178;

% If true, re-run process_raw_VMI on every file (also the global bg).
% If false, only re-process files that have no VMIdata_*.mat cached.
% Re-process whenever the chosen center changes, since process_raw_VMI
% bakes the center into the cached .mat.
REPROCESS_RAW = true;

% Output location. Matches the existing shipped file path so the Python
% paper_v2 export reads the new average without any path change.
script_dir = fileparts(mfilename('fullpath'));
OUT_DIR    = ['T:\github synchronized\VMI_matlab\matfile_data_scripts\A_state_paper_figures_single_pulse\high_snr'];
OUT_NAME   = 'ressumI2HeNI^+He';   % .mat and .fig share this stem

% Measurement set: I+He droplet, low doping, 300 mW, 03.12.24.
% Identical to the simulation_image.m reference set (line 153) and a
% superset of the paper_cov triplet (45668, 45662, 45667) by the addition
% of 45686.
SET_NAME    = 'I_2He_N:I^+He';
FILENUMBERS = [45668, 45662, 45667, 45686];   % 4 IDs
% =========================================================

if ~exist(OUT_DIR, 'dir')
    mkdir(OUT_DIR);
end

% --------- center inspection mode ---------
if INSPECT_CENTERS_ONLY
    if exist('autocenter_from_extended_data', 'file') ~= 2
        error(['autocenter_from_extended_data is not on the MATLAB path. ', ...
            'Add the legacy VMI toolbox before running the inspection pass.']);
    end
    fprintf('Per-file auto-detected centers (%s, %d files):\n', ...
            SET_NAME, numel(FILENUMBERS));
    centers = zeros(numel(FILENUMBERS), 2);
    for k = 1:numel(FILENUMBERS)
        fn_k = FILENUMBERS(k);
        c_k  = autocenter_from_extended_data(fn_k);
        centers(k, :) = c_k(:).';
        fprintf('  fn=%d  center=[%.4f, %.4f]  diff_from_VM_CENTER=[%+.3f, %+.3f]\n', ...
                fn_k, c_k(1), c_k(2), ...
                c_k(1) - VM_CENTER(1), c_k(2) - VM_CENTER(2));
    end
    fprintf('\nSummary:\n');
    fprintf('  median = [%.4f, %.4f]\n', median(centers(:,1)), median(centers(:,2)));
    fprintf('  mean   = [%.4f, %.4f]\n', mean(centers(:,1)),   mean(centers(:,2)));
    fprintf('  std    = [%.4f, %.4f]\n', std(centers(:,1)),    std(centers(:,2)));
    fprintf('  max abs deviation from VM_CENTER = [%.3f, %.3f]\n', ...
            max(abs(centers(:,1) - VM_CENTER(1))), ...
            max(abs(centers(:,2) - VM_CENTER(2))));

    figure('Name', sprintf('Auto-centers: %s', SET_NAME));
    scatter(centers(:,1), centers(:,2), 80, 1:numel(FILENUMBERS), 'filled');
    hold on;
    scatter(VM_CENTER(1), VM_CENTER(2), 200, 'k', 'x', 'LineWidth', 2);
    for k = 1:numel(FILENUMBERS)
        text(centers(k,1) + 0.2, centers(k,2) + 0.2, ...
             sprintf('%d', FILENUMBERS(k)), 'FontSize', 8);
    end
    xlabel('center_x / px');
    ylabel('center_y / px');
    title(sprintf('Per-file auto-detected centers (%s, x = hardcoded VM\\_CENTER)', ...
                  SET_NAME), 'Interpreter', 'tex');
    axis equal; grid on; colorbar;
    return;
end

% --------- averaging mode ---------
fprintf('Averaging %d I+He files at hardcoded center [%.4f, %.4f]\n', ...
        numel(FILENUMBERS), VM_CENTER(1), VM_CENTER(2));

% Ensure every raw file (including the global bg) has been processed at
% the chosen center.
required = FILENUMBERS;
if ~isnan(GLOBAL_BG_FILENUMBER)
    required = [required, GLOBAL_BG_FILENUMBER];
end

if REPROCESS_RAW
    unprocessed = required;
else
    unprocessed = [];
    for fn = required
        cache_path = fullfile(path_to_processed_data, ...
                              sprintf('VMIdata_%d.mat', fn));
        if ~exist(cache_path, 'file')
            unprocessed = [unprocessed, fn]; %#ok<AGROW>
        end
    end
end
if ~isempty(unprocessed)
    fprintf('Running process_raw_VMI on %d file(s)...\n', numel(unprocessed));
    process_raw_VMI(unprocessed, VM_CENTER);
end

% Load global bg once.
if ~isnan(GLOBAL_BG_FILENUMBER)
    res_global_bg = plot_processed_VMI( ...
        GLOBAL_BG_FILENUMBER, 0, VM_CENTER, true);
end

% Average over all files.
res_sum = [];
for k = 1:numel(FILENUMBERS)
    fn_k = FILENUMBERS(k);
    res_k = plot_processed_VMI(fn_k, 0, VM_CENTER, true);
    if ~isnan(GLOBAL_BG_FILENUMBER)
        res_k = subtract_processed_data(res_k, res_global_bg);
    end
    if k == 1
        res_sum = res_k;
    else
        res_sum = add_processed_data(res_sum, res_k);
    end
end
res_sum = multiply_processed_data(res_sum, 1 / numel(FILENUMBERS));

% Preserve original variable names so existing consumers
% (export_paper_v2_reference_data.m) do not need to change.
filenumbers = FILENUMBERS; %#ok<NASGU>

out_mat = fullfile(OUT_DIR, [OUT_NAME, '.mat']);
save(out_mat, 'res_sum', 'filenumbers');
fprintf('Saved averaged res_sum (n=%d) to %s\n', ...
        numel(FILENUMBERS), out_mat);

% Diagnostic figure.
fig = figure('Name', sprintf('%s averaged (n=%d)', SET_NAME, numel(FILENUMBERS)));
surf((res_sum.Y - res_sum.image_center_y) * VELOCITY_FACTOR, ...
     (res_sum.X - res_sum.image_center_x) * VELOCITY_FACTOR, ...
     res_sum.image);
view(90, 90);
xlabel('v_x / m/s');
ylabel('v_y / m/s');
xlim([-3000, 3000]);
ylim([-3000, 3000]);
if exist('colorcet', 'file') == 2
    colormap(colorcet('L08'));
end
pbaspect([1, 1, 1]);
title(sprintf('%s averaged (n=%d, vf=%.4f, bg=%d)', ...
              SET_NAME, numel(FILENUMBERS), VELOCITY_FACTOR, ...
              GLOBAL_BG_FILENUMBER), 'Interpreter', 'tex');

out_fig = fullfile(OUT_DIR, [OUT_NAME, '.fig']);
savefig(fig, out_fig);
fprintf('Saved diagnostic figure to %s\n', out_fig);
