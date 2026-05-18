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
%   The same high-SNR directory also provides a separate I+He2 res_sum:
%       A_state_paper_figures_single_pulse\high_snr\ressumI2HeNI^+He2
%   and plots:
%     surf((res.Y-res.image_center_y)*vf_single/100, ...
%          (res.X-res.image_center_x)*vf_single/100, res.image)
%   The exported Python reference keeps the same processed intensity but uses
%   Matplotlib-ready axis names: vx_Aps is the plot x-grid and vy_Aps is the
%   plot y-grid for pcolormesh(vx_Aps, vy_Aps, intensity).
%
% Output conventions:
%   Radial CSV columns:       v_mps,signal_arb
%   Phi CSV columns:          phi_rad,signal_arb
%   Cartesian image MAT:      vx_mps, vy_mps, intensity (Matplotlib pcolormesh layout)
%   Polar image MAT:          phi_rad, v_radius_mps, intensity_polar
%                             (rows = phi, cols = v_radius; mirrors res.image_polar)
%   Image metadata:           JSON sidecar per image with source, units, center,
%                             and scaling
%
% Velocity unit: m/s (canonical on-disk format for all reference CSVs).
% The Python loaders convert to A/ps internally so the documented
% A/ps binning conventions stay valid.
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

vf_single = 8.6178; % canonical VMI toolbox calibration (m/s per radial pixel).
                    % Earlier versions of this script used 10.0995, which is the
                    % outdated calibration baked into the original paper-v2 MATLAB
                    % source; the corrected value matches paper_v4 and the
                    % vmi_summary exporter.
common_center = [524.5297 380.8430];
gas_phase_reference_measurement = 43562;
gas_center = [482.9299 392.4866];

res_Iplus_gas = plot_processed_VMI(gas_phase_reference_measurement, true, gas_center, true);
v_mps = res_Iplus_gas.r(:) * vf_single;
signal_arb = res_Iplus_gas.radial_distribution(:);
writetable( ...
    table(v_mps, signal_arb, 'VariableNames', {'v_mps', 'signal_arb'}), ...
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
    v_mps = res_extra_ihe.r(:) * vf_single;
    signal_arb = res_extra_ihe.radial_distribution(:);
    writetable( ...
        table(v_mps, signal_arb, 'VariableNames', {'v_mps', 'signal_arb'}), ...
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

v_mps = res.r(:) * vf_single;
signal_arb = res.radial_distribution(:);
writetable( ...
    table(v_mps, signal_arb, 'VariableNames', {'v_mps', 'signal_arb'}), ...
    fullfile(out_dir, 'iplus_he_high_snr_radial.csv'));
fprintf('Exported high-SNR I+He radial graph curve.\n');

vx_mps = (res.X - res.image_center_x) * vf_single;
vy_mps = (res.Y - res.image_center_y) * vf_single;
intensity = res.image;

b_r = res.r * vf_single > 0;
phi_rad = res.phi(:);
signal_arb = mean(res.image_polar(:, b_r), 2);
signal_arb = signal_arb(:) / max(signal_arb(:));
writetable( ...
    table(phi_rad, signal_arb, 'VariableNames', {'phi_rad', 'signal_arb'}), ...
    fullfile(out_dir, 'iplus_he_high_snr_phi.csv'));

% Full 2-D polar image (phi rows, v_radius cols) used for the side-by-side
% polar comparison figure in plot_paper_v2.py. Keeps the existing 1-D phi
% CSV intact; the polar image is a separate, richer reference.
v_radius_mps = res.r(:) * vf_single;
phi_rad = res.phi(:);
intensity_polar = res.image_polar;
save(fullfile(image_dir, 'iplus_he_high_snr_vmi_polar_image.mat'), ...
    'phi_rad', 'v_radius_mps', 'intensity_polar');

fid = fopen(fullfile(image_dir, 'iplus_he_high_snr_vmi_polar_image.json'), 'w');
if fid < 0
    error('Could not open polar metadata JSON for writing.');
end
fprintf(fid, '{\n');
fprintf(fid, '  "legacy_script": "post_process_single_pulse_paper_IplusHe_comparison.m",\n');
fprintf(fid, '  "source_mat_file": "%s",\n', strrep(high_snr_load_path, '\', '\\'));
fprintf(fid, '  "channel": "I+He high-SNR processed VMI (polar)",\n');
fprintf(fid, '  "vf_single": %.15g,\n', vf_single);
fprintf(fid, '  "units": "m/s",\n');
fprintf(fid, '  "axis_equations": {\n');
fprintf(fid, '    "phi_rad": "res.phi",\n');
fprintf(fid, '    "v_radius_mps": "res.r * vf_single"\n');
fprintf(fid, '  },\n');
fprintf(fid, '  "image_center_x": %.15g,\n', res.image_center_x);
fprintf(fid, '  "image_center_y": %.15g,\n', res.image_center_y);
fprintf(fid, '  "matrix_layout": "rows=phi, cols=v_radius",\n');
fprintf(fid, '  "output_fields": ["phi_rad", "v_radius_mps", "intensity_polar"],\n');
fprintf(fid, '  "external_requirements": ["legacy VMI MATLAB toolbox", "high-SNR res_sum MAT file"]\n');
fprintf(fid, '}\n');
fclose(fid);

preview_polar = figure('Visible', 'off');
pcolor(phi_rad, v_radius_mps, intensity_polar.');
shading flat;
xlabel('phi / rad');
ylabel('v / m/s');
xlim([0, 2 * pi]);
ylim([0, max(v_radius_mps)]);
colorbar;
print(preview_polar, fullfile(image_dir, 'iplus_he_high_snr_vmi_polar_image_preview.png'), '-dpng', '-r150');
close(preview_polar);

save(fullfile(image_dir, 'iplus_he_high_snr_vmi_image.mat'), ...
    'vx_mps', 'vy_mps', 'intensity');

fid = fopen(fullfile(image_dir, 'iplus_he_high_snr_vmi_image.json'), 'w');
if fid < 0
    error('Could not open metadata JSON for writing.');
end
fprintf(fid, '{\n');
fprintf(fid, '  "legacy_script": "post_process_single_pulse_paper_IplusHe_comparison.m",\n');
fprintf(fid, '  "source_mat_file": "%s",\n', strrep(high_snr_load_path, '\', '\\'));
fprintf(fid, '  "channel": "I+He high-SNR processed VMI",\n');
fprintf(fid, '  "vf_single": %.15g,\n', vf_single);
fprintf(fid, '  "units": "m/s",\n');
fprintf(fid, '  "axis_equations": {\n');
fprintf(fid, '    "vx_mps": "(res.X - res.image_center_x) * vf_single",\n');
fprintf(fid, '    "vy_mps": "(res.Y - res.image_center_y) * vf_single"\n');
fprintf(fid, '  },\n');
fprintf(fid, '  "image_center_x": %.15g,\n', res.image_center_x);
fprintf(fid, '  "image_center_y": %.15g,\n', res.image_center_y);
fprintf(fid, '  "output_fields": ["vx_mps", "vy_mps", "intensity"],\n');
fprintf(fid, '  "external_requirements": ["legacy VMI MATLAB toolbox", "high-SNR res_sum MAT file"]\n');
fprintf(fid, '}\n');
fclose(fid);

preview = figure('Visible', 'off');
surf(vx_mps, vy_mps, intensity, 'EdgeColor', 'none');
view(90, 90);
pbaspect([1, 1, 1]);
xlim([-3500, 3500]);
ylim([-3500, 3500]);
xlabel('v_x / m/s');
ylabel('v_y / m/s');
colorbar;
print(preview, fullfile(image_dir, 'iplus_he_high_snr_vmi_image_preview.png'), '-dpng', '-r150');
close(preview);

% Export the analogous high-SNR I+He2 data for separate Mass 135 comparisons.
high_snr_he2_path = ['T:\github synchronized\VMI_matlab\matfile_data_scripts\', ...
    'A_state_paper_figures_single_pulse\high_snr\ressumI2HeNI^+He2'];
high_snr_he2_load_path = high_snr_he2_path;
if exist(high_snr_he2_load_path, 'file') ~= 2 && exist([high_snr_he2_path, '.mat'], 'file') == 2
    high_snr_he2_load_path = [high_snr_he2_path, '.mat'];
end
if exist(high_snr_he2_load_path, 'file') ~= 2
    error('High-SNR I+He2 MAT file not found: %s', high_snr_he2_path);
end
data_in = load(high_snr_he2_load_path);
if ~isfield(data_in, 'res_sum')
    error('High-SNR I+He2 MAT file does not contain expected variable res_sum: %s', high_snr_he2_load_path);
end
res = data_in.res_sum;

v_mps = res.r(:) * vf_single;
signal_arb = res.radial_distribution(:);
writetable( ...
    table(v_mps, signal_arb, 'VariableNames', {'v_mps', 'signal_arb'}), ...
    fullfile(out_dir, 'iplus_he2_high_snr_radial.csv'));
fprintf('Exported high-SNR I+He2 radial graph curve.\n');

b_r = res.r * vf_single > 0;
phi_rad = res.phi(:);
signal_arb = mean(res.image_polar(:, b_r), 2);
signal_arb = signal_arb(:) / max(signal_arb(:));
writetable( ...
    table(phi_rad, signal_arb, 'VariableNames', {'phi_rad', 'signal_arb'}), ...
    fullfile(out_dir, 'iplus_he2_high_snr_phi.csv'));

v_radius_mps = res.r(:) * vf_single;
phi_rad = res.phi(:);
intensity_polar = res.image_polar;
save(fullfile(image_dir, 'iplus_he2_high_snr_vmi_polar_image.mat'), ...
    'phi_rad', 'v_radius_mps', 'intensity_polar');

fid = fopen(fullfile(image_dir, 'iplus_he2_high_snr_vmi_polar_image.json'), 'w');
if fid < 0
    error('Could not open I+He2 polar metadata JSON for writing.');
end
fprintf(fid, '{\n');
fprintf(fid, '  "legacy_script": "post_process_single_pulse_paper_IplusHe_comparison.m",\n');
fprintf(fid, '  "source_mat_file": "%s",\n', strrep(high_snr_he2_load_path, '\', '\\'));
fprintf(fid, '  "channel": "I+He2 high-SNR processed VMI (polar)",\n');
fprintf(fid, '  "vf_single": %.15g,\n', vf_single);
fprintf(fid, '  "units": "m/s",\n');
fprintf(fid, '  "axis_equations": {\n');
fprintf(fid, '    "phi_rad": "res.phi",\n');
fprintf(fid, '    "v_radius_mps": "res.r * vf_single"\n');
fprintf(fid, '  },\n');
fprintf(fid, '  "image_center_x": %.15g,\n', res.image_center_x);
fprintf(fid, '  "image_center_y": %.15g,\n', res.image_center_y);
fprintf(fid, '  "matrix_layout": "rows=phi, cols=v_radius",\n');
fprintf(fid, '  "output_fields": ["phi_rad", "v_radius_mps", "intensity_polar"],\n');
fprintf(fid, '  "external_requirements": ["legacy VMI MATLAB toolbox", "high-SNR res_sum MAT file"]\n');
fprintf(fid, '}\n');
fclose(fid);

preview_polar = figure('Visible', 'off');
pcolor(phi_rad, v_radius_mps, intensity_polar.');
shading flat;
xlabel('phi / rad');
ylabel('v / m/s');
xlim([0, 2 * pi]);
ylim([0, max(v_radius_mps)]);
colorbar;
print(preview_polar, fullfile(image_dir, 'iplus_he2_high_snr_vmi_polar_image_preview.png'), '-dpng', '-r150');
close(preview_polar);

vx_mps = (res.X - res.image_center_x) * vf_single;
vy_mps = (res.Y - res.image_center_y) * vf_single;
intensity = res.image;
save(fullfile(image_dir, 'iplus_he2_high_snr_vmi_image.mat'), ...
    'vx_mps', 'vy_mps', 'intensity');

fid = fopen(fullfile(image_dir, 'iplus_he2_high_snr_vmi_image.json'), 'w');
if fid < 0
    error('Could not open I+He2 metadata JSON for writing.');
end
fprintf(fid, '{\n');
fprintf(fid, '  "legacy_script": "post_process_single_pulse_paper_IplusHe_comparison.m",\n');
fprintf(fid, '  "source_mat_file": "%s",\n', strrep(high_snr_he2_load_path, '\', '\\'));
fprintf(fid, '  "channel": "I+He2 high-SNR processed VMI",\n');
fprintf(fid, '  "vf_single": %.15g,\n', vf_single);
fprintf(fid, '  "units": "m/s",\n');
fprintf(fid, '  "axis_equations": {\n');
fprintf(fid, '    "vx_mps": "(res.X - res.image_center_x) * vf_single",\n');
fprintf(fid, '    "vy_mps": "(res.Y - res.image_center_y) * vf_single"\n');
fprintf(fid, '  },\n');
fprintf(fid, '  "image_center_x": %.15g,\n', res.image_center_x);
fprintf(fid, '  "image_center_y": %.15g,\n', res.image_center_y);
fprintf(fid, '  "output_fields": ["vx_mps", "vy_mps", "intensity"],\n');
fprintf(fid, '  "external_requirements": ["legacy VMI MATLAB toolbox", "high-SNR res_sum MAT file"]\n');
fprintf(fid, '}\n');
fclose(fid);

preview = figure('Visible', 'off');
surf(vx_mps, vy_mps, intensity, 'EdgeColor', 'none');
view(90, 90);
pbaspect([1, 1, 1]);
xlim([-3500, 3500]);
ylim([-3500, 3500]);
xlabel('v_x / m/s');
ylabel('v_y / m/s');
colorbar;
print(preview, fullfile(image_dir, 'iplus_he2_high_snr_vmi_image_preview.png'), '-dpng', '-r150');
close(preview);

fprintf('Exported paper-v2 high-SNR I+He2 2-D VMI image reference.\n');

fprintf('Exported paper-v2 high-SNR 2-D VMI image reference.\n');
