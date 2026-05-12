% export_paper_v3_reference_data.m
%
% Export experimental references used by
% legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper_v3.m
% into small CSV files that the Python port can load.
%
% Provenance:
% - Source script: post_process_single_pulse_paper_v3.m, active non-effusive
%   droplet branch, lines 114-222.
% - I+ droplet reference measurement: 43567 (600 mW), center [524.5297 380.8430].
% - I+ gas reference measurement: 43568 (600 mW), center [524.5297 380.8430].
% - I+He measurement: 43563 (300 mW), later replaced by high-SNR MAT data.
% - High-SNR MAT input:
%   T:\github synchronized\VMI_matlab\matfile_data_scripts\A_state_paper_figures_single_pulse\high_snr\ressumI2HeNI^+He
% - Timescan input: mean_timescan_2d_VMI([296:297], false,
%   [524.5297 380.8430], false), then columns with t > 150 ps.
% - Velocity factors: velocity_factor for I+He high-SNR data, 5.636 for
%   timescan data.
% - Angular threshold: VMIN_ANGULAR_DISTR = 0 m/s.
%
% Requirements:
% - The legacy VMI MATLAB toolbox must be on the MATLAB path.
% - The high-SNR MAT file path above must exist on the MATLAB machine.
%
% Outputs, written to data/reference/ when this script is run from this folder:
% - paper_v3_iplus_he_radial.csv: v_mps,signal_arb
% - paper_v3_timescan_radial.csv: v_mps,signal_t_001,...
% - paper_v3_iplus_he_phi.csv: phi_rad,signal_arb

clear; close all;

VMIN_ANGULAR_DISTR = 0;
VM_center_Iplus_He = [509.3664 387.6409];
VM_center_shared = [524.5297 380.8430];
I_plus_from_drop_reference_measurement = 43567; %#ok<NASGU>
gas_phase_reference_measurement = 43568; %#ok<NASGU>
I_plus_He_from_drop_reference_measurement = 43563; %#ok<NASGU>
vf_timescan = 5.636;

fprintf('Exporting paper v3 I+He references...\n');

global plot_processed_with_ROI
plot_processed_with_ROI = false;

% Keep these calls next to the provenance even though the active v3 path uses
% the high-SNR MAT file for res_Iplus_He after loading/subtracting.
res_Iplus_from_He = plot_processed_VMI(I_plus_from_drop_reference_measurement, true, VM_center_shared, true);
res_Iplus_gas = plot_processed_VMI(gas_phase_reference_measurement, true, VM_center_shared, true);
res_Iplus_from_He = subtract_processed_data(res_Iplus_from_He, res_Iplus_gas); %#ok<NASGU>
res_Iplus_He = plot_processed_VMI(I_plus_He_from_drop_reference_measurement, true, VM_center_Iplus_He, true); %#ok<NASGU>

data_in = load('T:\github synchronized\VMI_matlab\matfile_data_scripts\A_state_paper_figures_single_pulse\high_snr\ressumI2HeNI^+He');
res_Iplus_He = data_in.res_sum;

global velocity_factor
vf = velocity_factor;

global GS_bleach_correction
GS_bleach_correction = 1;
res2 = mean_timescan_2d_VMI([296:297], false, VM_center_shared, false);

out_dir = fileparts(mfilename('fullpath'));
out_dir = fullfile(out_dir, '..');

% Top panel I+He radial trace: v3 line 178.
v_he_mps = res_Iplus_He.r(:) * vf;
signal_he = res_Iplus_He.radial_distribution(:);
T_he = table(v_he_mps, signal_he, ...
    'VariableNames', {'v_mps', 'signal_arb'});
writetable(T_he, fullfile(out_dir, 'paper_v3_iplus_he_radial.csv'));
fprintf('Saved paper_v3_iplus_he_radial.csv\n');

% Top panel timescan traces: v3 lines 169-173.
timescan_mask = res2.t > 150;
timescan_signal = res2.data(:, timescan_mask);
v_timescan_mps = res2.r(:) * vf_timescan;
T_ts = array2table([v_timescan_mps, timescan_signal]);
signal_names = compose("signal_t_%03d", 1:size(timescan_signal, 2));
T_ts.Properties.VariableNames = [{'v_mps'}, cellstr(signal_names)];
writetable(T_ts, fullfile(out_dir, 'paper_v3_timescan_radial.csv'));
fprintf('Saved paper_v3_timescan_radial.csv\n');

% Bottom panel angular trace: v3 lines 195-201.
b_r = res_Iplus_He.r * vf > VMIN_ANGULAR_DISTR;
y_phi = mean(res_Iplus_He.image_polar(:, b_r), 2);
T_phi = table(res_Iplus_He.phi(:), y_phi(:), ...
    'VariableNames', {'phi_rad', 'signal_arb'});
writetable(T_phi, fullfile(out_dir, 'paper_v3_iplus_he_phi.csv'));
fprintf('Saved paper_v3_iplus_he_phi.csv\n');
