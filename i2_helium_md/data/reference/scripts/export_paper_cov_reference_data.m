% Export reference data for the Python paper-cov port.
%
% Legacy source:
%   legacy_matlab_repository/single_pulse_simulation/
%       post_process_single_pulse_paper_IplusHe_comparison_cov.m
%
% Scope:
%   Active non-effusive "comparison IHe+ covariance" branch only. The
%   experimental pair-covariance matrices are computed by the legacy VMI
%   toolbox function generate_VMI_covariance_matrices and exported as a
%   small, inspectable reference under data/reference/paper_cov/.
%   The Python script plot_paper_cov.py loads this reference and pairs
%   it with a simulated counterpart computed from the Python ion
%   checkpoint.
%
% Experimental measurements (I+He droplet, low doping, 300 mW, 03.12.24):
%   45668, 45662, 45667
%
% Output fields (canonical on-disk units, m/s for velocity axes):
%   cov_angular            : 2-D (N_theta x N_theta), diagonal already zeroed
%   cov_radial             : 2-D (N_v x N_v), diagonal already zeroed AND
%                            2 x 2 movmean already applied along each axis
%   theta_centers_rad      : 1-D length N_theta (radians)
%   velocity_centers_mps   : 1-D length N_v (m/s, mass-corrected)
%
% Companion CSV (separate file, same out_dir):
%   iplus_he_phi.csv       : phi_rad, signal_arb -- 1-D phi(angle)
%                            distribution literal-port of
%                            _cov.m lines 187-195 applied to the average
%                            of plot_processed_VMI over the same three
%                            measurement IDs as the covariance.
%
% Companion JSON sidecar (same basename, .json suffix) documents:
%   - source MATLAB function and arguments
%   - measurement IDs and center
%   - velocity factor and mass correction
%   - axis units and processing steps applied
%   - companion files (the phi CSV) and their provenance
%   - export timestamp
%
% Requirements:
%   - The legacy VMI MATLAB toolbox on the MATLAB path.
%   - generate_VMI_covariance_matrices available.
%   - autocenter_from_extended_data available.
%   - plot_processed_VMI, add_processed_data, multiply_processed_data
%     available (needed for the phi CSV companion file).
%   - physical_constants.m on the path (defines velocity_factor and
%     mass_correction_factor).
%   - No MATLAB Python bridge is required.

clear; close all;

fprintf('Exporting paper-cov experimental I+He pair-covariance reference...\n');

if exist('generate_VMI_covariance_matrices', 'file') ~= 2
    error(['generate_VMI_covariance_matrices is not on the MATLAB path. ', ...
        'Add the legacy VMI toolbox before running this exporter.']);
end
if exist('autocenter_from_extended_data', 'file') ~= 2
    error(['autocenter_from_extended_data is not on the MATLAB path. ', ...
        'Add the legacy VMI toolbox before running this exporter.']);
end
if exist('plot_processed_VMI', 'file') ~= 2
    error(['plot_processed_VMI is not on the MATLAB path. ', ...
        'Add the legacy VMI toolbox before running this exporter.']);
end
if exist('add_processed_data', 'file') ~= 2
    error(['add_processed_data is not on the MATLAB path. ', ...
        'Add the legacy VMI toolbox before running this exporter.']);
end
if exist('multiply_processed_data', 'file') ~= 2
    error(['multiply_processed_data is not on the MATLAB path. ', ...
        'Add the legacy VMI toolbox before running this exporter.']);
end
if exist('physical_constants', 'file') ~= 2
    error(['physical_constants.m is not on the MATLAB path. ', ...
        'Add the legacy single-pulse repository before running this exporter.']);
end

run physical_constants.m;
velocity_factor =  8.6178;
if ~exist('velocity_factor', 'var')
    error('physical_constants.m did not define velocity_factor.');
end
mass_correction_factor= sqrt((127) /(131));
if ~exist('mass_correction_factor', 'var')
    error('physical_constants.m did not define mass_correction_factor.');
end

script_dir = fileparts(mfilename('fullpath'));
out_dir = fullfile(script_dir, '..', 'paper_cov');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% Literal arguments from _cov.m (line 362).
I_plus_He_from_drop_reference_measurements = [45668, 45662, 45667];
apply_angular_filter = false;
event_filter = true;
theta_target = pi;
theta_range = 40 / 180 * pi;

fprintf('Auto-centering on I+He droplet measurements...\n');
center = autocenter_from_extended_data(I_plus_He_from_drop_reference_measurements);
fprintf('  center = [%.6f, %.6f]\n', center(1), center(2));

fprintf('Computing experimental covariance matrices...\n');
result = generate_VMI_covariance_matrices( ...
    I_plus_He_from_drop_reference_measurements, ...
    [0, 600], ...
    center, ...
    [90, 90], ...
    apply_angular_filter, ...
    event_filter, ...
    theta_target, ...
    theta_range);

% Apply the legacy post-processing in order (lines 367, 385, 388-389 of _cov.m).
cov_angular = result.cov_angular - diag(diag(result.cov_angular));
cov_radial = result.cov_radial - diag(diag(result.cov_radial));
cov_radial = movmean(cov_radial, 2, 1);
cov_radial = movmean(cov_radial, 2, 2);

theta_centers_rad = result.theta(:);

% Convert radial axis to m/s and apply mass correction (line 391 of _cov.m).
% MATLAB pixel -> A/ps is "* velocity_factor / 100"; on-disk we keep m/s so we
% drop the /100 (the Python loader divides by 100 to get A/ps).
velocity_centers_mps = result.r(:) * velocity_factor * mass_correction_factor;

mat_path = fullfile(out_dir, 'iplus_he_covariance.mat');
save(mat_path, ...
    'cov_angular', 'cov_radial', 'theta_centers_rad', 'velocity_centers_mps');
fprintf('Saved covariance matrices to %s\n', mat_path);

% Companion 1-D phi(angle) reference. Literal port of _cov.m lines 100-205
% so the Python paper_cov_phi_distribution.png figure overlays the same
% experimental shape MATLAB plots. The same three I+He droplet IDs are
% reused. The phi pipeline uses the HARDCODED center [524.5297, 380.8430]
% from _cov.m line 100, NOT the auto-detected center used for the
% covariance pipeline above. A different center pivots the polar
% transform and rotates the resulting angular distribution, so reusing
% the auto-detected center here produces a phi curve that does not match
% the live MATLAB _cov.m figure.
phi_center = [524.5297, 380.8430];
fprintf('Building averaged res_Iplus_He from plot_processed_VMI (center = [%.4f, %.4f])...\n', ...
    phi_center(1), phi_center(2));
for k = 1:length(I_plus_He_from_drop_reference_measurements)
    fn = I_plus_He_from_drop_reference_measurements(k);
    if k == 1
        res_Iplus_He = plot_processed_VMI(fn, true, phi_center, true);
    else
        res_Iplus_He = add_processed_data( ...
            res_Iplus_He, plot_processed_VMI(fn, true, phi_center, true));
    end
end
res_Iplus_He = multiply_processed_data( ...
    res_Iplus_He, 1 / length(I_plus_He_from_drop_reference_measurements));

% _cov.m lines 187-195 with VMIN_ANGULAR_DISTR = 0 (b_r is all-true).
b_r        = res_Iplus_He.r * velocity_factor > 0;
phi_rad    = res_Iplus_He.phi(:);
signal_arb = mean(res_Iplus_He.image_polar(:, b_r), 2);
signal_arb = signal_arb(:) / max(signal_arb(:));

phi_csv_path = fullfile(out_dir, 'iplus_he_phi.csv');
writetable( ...
    table(phi_rad, signal_arb, 'VariableNames', {'phi_rad', 'signal_arb'}), ...
    phi_csv_path);
fprintf('Saved phi(angle) reference to %s\n', phi_csv_path);

% Sidecar JSON metadata.
json_path = fullfile(out_dir, 'iplus_he_covariance.json');
fid = fopen(json_path, 'w');
if fid < 0
    error('Could not open metadata JSON for writing: %s', json_path);
end
fprintf(fid, '{\n');
fprintf(fid, '  "legacy_script": "post_process_single_pulse_paper_IplusHe_comparison_cov.m",\n');
fprintf(fid, '  "exporter": "data/reference/scripts/export_paper_cov_reference_data.m",\n');
fprintf(fid, '  "channel": "I+He droplet pair covariance",\n');
fprintf(fid, '  "measurement_ids": [45668, 45662, 45667],\n');
fprintf(fid, '  "doping": "low",\n');
fprintf(fid, '  "probe_power_mW": 300,\n');
fprintf(fid, '  "measurement_date": "2024-12-03",\n');
fprintf(fid, '  "center_pixels_covariance": [%.15g, %.15g],\n', center(1), center(2));
fprintf(fid, '  "center_pixels_phi": [%.15g, %.15g],\n', phi_center(1), phi_center(2));
fprintf(fid, '  "velocity_factor": %.15g,\n', velocity_factor);
fprintf(fid, '  "mass_correction_factor": %.15g,\n', mass_correction_factor);
fprintf(fid, '  "covariance_arguments": {\n');
fprintf(fid, '    "velocity_range_px": [0, 600],\n');
fprintf(fid, '    "bins": [90, 90],\n');
fprintf(fid, '    "apply_angular_filter": false,\n');
fprintf(fid, '    "event_filter": true,\n');
fprintf(fid, '    "theta_target_rad": %.15g,\n', theta_target);
fprintf(fid, '    "theta_range_rad": %.15g\n', theta_range);
fprintf(fid, '  },\n');
fprintf(fid, '  "processing_steps": [\n');
fprintf(fid, '    "cov_angular -= diag(diag(cov_angular))",\n');
fprintf(fid, '    "cov_radial -= diag(diag(cov_radial))",\n');
fprintf(fid, '    "cov_radial = movmean(cov_radial, 2, 1)",\n');
fprintf(fid, '    "cov_radial = movmean(cov_radial, 2, 2)",\n');
fprintf(fid, '    "res_Iplus_He = mean(plot_processed_VMI(ids, true, [524.5297 380.8430], true))",\n');
fprintf(fid, '    "phi_signal = mean(res_Iplus_He.image_polar(:, res_Iplus_He.r * vf > 0), 2) / max(...)"\n');
fprintf(fid, '  ],\n');
fprintf(fid, '  "companion_files": [\n');
fprintf(fid, '    {\n');
fprintf(fid, '      "filename": "iplus_he_phi.csv",\n');
fprintf(fid, '      "purpose": "experimental I+He phi(angle) distribution overlay for paper_cov_phi_distribution.png",\n');
fprintf(fid, '      "columns": ["phi_rad", "signal_arb"],\n');
fprintf(fid, '      "source": "mean of plot_processed_VMI for the same 3 measurement IDs at hardcoded center [524.5297, 380.8430] (matches _cov.m line 100), polar mean over all radii"\n');
fprintf(fid, '    }\n');
fprintf(fid, '  ],\n');
fprintf(fid, '  "axis_units": {\n');
fprintf(fid, '    "theta_centers_rad": "radian (range [-pi, pi])",\n');
fprintf(fid, '    "velocity_centers_mps": "m/s, mass-corrected by mass_correction_factor"\n');
fprintf(fid, '  },\n');
fprintf(fid, '  "axis_equations": {\n');
fprintf(fid, '    "theta_centers_rad": "result.theta",\n');
fprintf(fid, '    "velocity_centers_mps": "result.r * velocity_factor * mass_correction_factor"\n');
fprintf(fid, '  },\n');
fprintf(fid, '  "external_requirements": ["legacy VMI MATLAB toolbox"],\n');
fprintf(fid, '  "export_timestamp_utc": "%s"\n', datestr(datetime('now', 'TimeZone', 'UTC'), 'yyyy-mm-ddTHH:MM:SSZ'));
fprintf(fid, '}\n');
fclose(fid);
fprintf('Saved metadata sidecar to %s\n', json_path);

% Preview PNGs for quick visual inspection (not loaded by Python).
preview_angular = figure('Visible', 'off');
imagesc(theta_centers_rad, theta_centers_rad, cov_angular);
axis equal tight;
xlabel('theta / rad');
ylabel('theta / rad');
colorbar;
title('Experimental angular pair covariance (diag zeroed)');
print(preview_angular, fullfile(out_dir, 'iplus_he_cov_angular_preview.png'), '-dpng', '-r150');
close(preview_angular);

preview_radial = figure('Visible', 'off');
imagesc(velocity_centers_mps / 100, velocity_centers_mps / 100, cov_radial);
axis equal tight;
xlabel('v / A/ps');
ylabel('v / A/ps');
colorbar;
title('Experimental radial pair covariance (diag zeroed + movmean)');
print(preview_radial, fullfile(out_dir, 'iplus_he_cov_radial_preview.png'), '-dpng', '-r150');
close(preview_radial);

fprintf('Exported paper-cov experimental pair-covariance reference.\n');
