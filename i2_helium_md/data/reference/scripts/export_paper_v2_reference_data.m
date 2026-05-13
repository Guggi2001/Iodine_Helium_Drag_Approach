% Export reference data for the Python paper-v2 I+He comparison port.
%
% Legacy source:
%   legacy_matlab_repository/single_pulse_simulation/
%       post_process_single_pulse_paper_IplusHe_comparison.m
%
% Scope:
%   Active non-effusive "comparison IHe+" branch only. This exporter uses
%   MATLAB/VMI toolbox processing as the source of truth for experimental
%   data and writes small Python-consumable references under
%   data/reference/paper_v2/.
%
% Experimental radial references actually plotted by the MATLAB radial panel:
%   43562: I+ gas, 300 mW, center [482.9299 392.4866]
%   high-SNR res_sum: I+He, loaded from the MAT file below
% Extra I+He power-scan radial references exported for comparison:
%   43556: I+He, 160 mW, center [524.5297 380.8430]
%   43569: I+He, 600 mW, center [524.5297 380.8430]
%
% 2-D VMI image reference:
%   The MATLAB figure's experimental image panel uses the processed
%   high-SNR res_sum loaded from:
%     T:\github synchronized\VMI_matlab\matfile_data_scripts\
%       A_state_paper_figures_single_pulse\high_snr\ressumI2HeNI^+He
%   and plots:
%     surf((res.Y-res.image_center_y)*vf_single/100, ...
%          (res.X-res.image_center_x)*vf_single/100, res.image)
%   The exported Python reference keeps the same processed intensity but uses
%   Matplotlib-ready axis names: vx_Aps is the plot x-grid and vy_Aps is the
%   plot y-grid for pcolormesh(vx_Aps, vy_Aps, intensity).
%
% Output conventions:
%   Radial CSV columns: v_Aps,signal_arb
%   Phi CSV columns:    phi_rad,signal_arb
%   Image MAT fields:  vx_Aps, vy_Aps, intensity
%                      vx_Aps and vy_Aps are full calibrated coordinate grids,
%                      normalized for Matplotlib pcolormesh axes.
%   Image metadata:    JSON sidecar with source, units, center, and scaling
%
% Requirements:
%   - The legacy VMI MATLAB toolbox on the MATLAB path.
%   - plot_processed_VMI available.
%   - colorcet is optional, only used by the original plot.
%   - No MATLAB Python bridge is required.

clear; close all;

fprintf('Exporting paper v2 I+He comparison references...\n');

if exist('plot_processed_VMI', 'file') ~= 2
    error(['plot_processed_VMI is not on the MATLAB path. ', ...
        'Add the legacy VMI toolbox before running this exporter.']);
end

script_dir = fileparts(mfilename('fullpath'));
out_dir = fullfile(script_dir, '..', 'paper_v2');
image_dir = fullfile(out_dir, 'images');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
if ~exist(image_dir, 'dir')
    mkdir(image_dir);
end

vf_single = 10.0995; % same paper-v2 script scale; radial/image axes in A/ps use /100
common_center = [524.5297 380.8430];
gas_phase_reference_measurement = 43562;
gas_center = [482.9299 392.4866];

res_Iplus_gas = plot_processed_VMI(gas_phase_reference_measurement, true, gas_center, true);
v_Aps = res_Iplus_gas.r(:) * vf_single / 100;
signal_arb = res_Iplus_gas.radial_distribution(:);
writetable( ...
    table(v_Aps, signal_arb, 'VariableNames', {'v_Aps', 'signal_arb'}), ...
    fullfile(out_dir, 'iplus_gas_300mw_43562_radial.csv'));
fprintf('Exported I+ gas radial graph curve from measurement %d\n', gas_phase_reference_measurement);

extra_ihe_specs = {
    43556, 'iplus_he_160mw_43556_radial.csv', 160;
    43569, 'iplus_he_600mw_43569_radial.csv', 600;
};

for i = 1:size(extra_ihe_specs, 1)
    measurement_id = extra_ihe_specs{i, 1};
    filename = extra_ihe_specs{i, 2};
    power_mw = extra_ihe_specs{i, 3};

    res_extra_ihe = plot_processed_VMI(measurement_id, true, common_center, true);
    v_Aps = res_extra_ihe.r(:) * vf_single / 100;
    signal_arb = res_extra_ihe.radial_distribution(:);
    writetable( ...
        table(v_Aps, signal_arb, 'VariableNames', {'v_Aps', 'signal_arb'}), ...
        fullfile(out_dir, filename));
    fprintf('Exported optional I+He radial curve from measurement %d (%d mW)\n', ...
        measurement_id, power_mw);
end

% Export the processed high-SNR image used by the MATLAB image panel.
high_snr_path = ['T:\github synchronized\VMI_matlab\matfile_data_scripts\', ...
    'A_state_paper_figures_single_pulse\high_snr\ressumI2HeNI^+He'];
high_snr_load_path = high_snr_path;
if exist(high_snr_load_path, 'file') ~= 2 && exist([high_snr_path, '.mat'], 'file') == 2
    high_snr_load_path = [high_snr_path, '.mat'];
end
if exist(high_snr_load_path, 'file') ~= 2
    error('High-SNR MAT file not found: %s', high_snr_path);
end
data_in = load(high_snr_load_path);
if ~isfield(data_in, 'res_sum')
    error('High-SNR MAT file does not contain expected variable res_sum: %s', high_snr_load_path);
end
res = data_in.res_sum;

v_Aps = res.r(:) * vf_single / 100;
signal_arb = res.radial_distribution(:);
writetable( ...
    table(v_Aps, signal_arb, 'VariableNames', {'v_Aps', 'signal_arb'}), ...
    fullfile(out_dir, 'iplus_he_high_snr_radial.csv'));
fprintf('Exported high-SNR I+He radial graph curve.\n');

vx_Aps = (res.X - res.image_center_x) * vf_single / 100;
vy_Aps = (res.Y - res.image_center_y) * vf_single / 100;
intensity = res.image;

b_r = res.r * vf_single > 0;
phi_rad = res.phi(:);
signal_arb = mean(res.image_polar(:, b_r), 2);
signal_arb = signal_arb(:) / max(signal_arb(:));
writetable( ...
    table(phi_rad, signal_arb, 'VariableNames', {'phi_rad', 'signal_arb'}), ...
    fullfile(out_dir, 'iplus_he_high_snr_phi.csv'));

save(fullfile(image_dir, 'iplus_he_high_snr_vmi_image.mat'), ...
    'vx_Aps', 'vy_Aps', 'intensity');

fid = fopen(fullfile(image_dir, 'iplus_he_high_snr_vmi_image.json'), 'w');
if fid < 0
    error('Could not open metadata JSON for writing.');
end
fprintf(fid, '{\n');
fprintf(fid, '  "legacy_script": "post_process_single_pulse_paper_IplusHe_comparison.m",\n');
fprintf(fid, '  "source_mat_file": "%s",\n', strrep(high_snr_load_path, '\', '\\'));
fprintf(fid, '  "channel": "I+He high-SNR processed VMI",\n');
fprintf(fid, '  "vf_single": %.15g,\n', vf_single);
fprintf(fid, '  "units": "A/ps",\n');
fprintf(fid, '  "axis_equations": {\n');
fprintf(fid, '    "vx_Aps": "(res.X - res.image_center_x) * vf_single / 100",\n');
fprintf(fid, '    "vy_Aps": "(res.Y - res.image_center_y) * vf_single / 100"\n');
fprintf(fid, '  },\n');
fprintf(fid, '  "image_center_x": %.15g,\n', res.image_center_x);
fprintf(fid, '  "image_center_y": %.15g,\n', res.image_center_y);
fprintf(fid, '  "output_fields": ["vx_Aps", "vy_Aps", "intensity"],\n');
fprintf(fid, '  "external_requirements": ["legacy VMI MATLAB toolbox", "high-SNR res_sum MAT file"]\n');
fprintf(fid, '}\n');
fclose(fid);

preview = figure('Visible', 'off');
surf(vx_Aps, vy_Aps, intensity, 'EdgeColor', 'none');
view(90, 90);
pbaspect([1, 1, 1]);
xlim([-35, 35]);
ylim([-35, 35]);
xlabel('v_x / A/ps');
ylabel('v_y / A/ps');
colorbar;
print(preview, fullfile(image_dir, 'iplus_he_high_snr_vmi_image_preview.png'), '-dpng', '-r150');
close(preview);

fprintf('Exported paper-v2 high-SNR 2-D VMI image reference.\n');
