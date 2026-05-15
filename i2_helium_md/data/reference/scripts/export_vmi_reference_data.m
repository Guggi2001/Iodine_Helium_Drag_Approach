% export_vmi_reference_data.m
%
% Exports the Abel-inverted, image-smoothed, mass-corrected VMI summary
% references used by scripts/post_processing/plot_experimental_comparison.py
% and the consolidated scripts/post_processing/plot_run_summary.py.
%
% This pipeline is intentionally different from the paper_v2/v3/v4 exporters:
% it averages four I+He raw measurements, runs movmean image smoothing, does
% Abel inversion on both channels, and applies a sqrt(127/131) mass correction
% to the I+He radial axis. The output files therefore are not directly
% comparable to the paper_v2/v3/v4 radial exports.
%
% Outputs, written under data/reference/vmi_summary/:
% - vmi_iplus_he.csv:  v_mps,signal_arb
% - vmi_iplus_gas.csv: v_mps,signal_arb
clear; close all;

% NOTE: Ensure your VMI toolbox is on the MATLAB path before running this!

vf_single = 8.6178;
mass_correction_factor = sqrt(127/131);

script_dir = fileparts(mfilename('fullpath'));
out_dir = fullfile(script_dir, '..', 'vmi_summary');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

fprintf('Processing I+He Droplet data...\n');
%% 1. Process I+He Measurements
I_plus_He_measurements = [45668, 45662, 45667, 45686];
for k = 1:length(I_plus_He_measurements)
    fn = I_plus_He_measurements(k);
    if k == 1
        res_Iplus_He = plot_processed_VMI(fn, true, [524.5297 380.8430], true);
    else
        res_temp = plot_processed_VMI(fn, true, [524.5297 380.8430], true);
        res_Iplus_He = add_processed_data(res_Iplus_He, res_temp);
    end
end
res_Iplus_He = multiply_processed_data(res_Iplus_He, 1/length(I_plus_He_measurements));
res_Iplus_He.image(res_Iplus_He.image < 0) = 0;
res_Iplus_He.image = movmean(res_Iplus_He.image, 3, 1);
res_Iplus_He.image = movmean(res_Iplus_He.image, 3, 2);
res_Iplus_He = abel_invert_processed_VMI(res_Iplus_He);

% Scale arrays as done in the original plot. Output is m/s (the canonical
% on-disk velocity unit for all reference CSVs). Python loaders convert to
% A/ps internally.
v_he_mps = res_Iplus_He.r * vf_single * mass_correction_factor;
signal_he = movmean(res_Iplus_He.radial_distribution, 1);

% Export
T_he = table(v_he_mps(:), signal_he(:), 'VariableNames', {'v_mps', 'signal_arb'});
writetable(T_he, fullfile(out_dir, 'vmi_iplus_he.csv'));
fprintf('Saved vmi_iplus_he.csv\n\n');


fprintf('Processing I+ Gas Phase data...\n');
%% 2. Process Gas Phase I+
gas_measurement = 43632;
res_Iplus_gas = plot_processed_VMI(gas_measurement, 1, [482.9299 392.4866], true);
res_Iplus_gas = abel_invert_processed_VMI(res_Iplus_gas);

% Scale arrays
v_gas_mps = res_Iplus_gas.r * vf_single;
signal_gas = res_Iplus_gas.radial_distribution;

% Export
T_gas = table(v_gas_mps(:), signal_gas(:), 'VariableNames', {'v_mps', 'signal_arb'});
writetable(T_gas, fullfile(out_dir, 'vmi_iplus_gas.csv'));
fprintf('Saved vmi_iplus_gas.csv\n');