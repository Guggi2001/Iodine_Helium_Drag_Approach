% export_paper_v4_reference_data.m
%
% Export experimental radial references used by
% legacy_matlab_repository/single_pulse_simulation/post_process_single_pulse_paper_v4.m
% into small CSV files that the Python v4 port can load.
%
% Provenance:
% - Source script: post_process_single_pulse_paper_v4.m, active non-effusive
%   droplet branch, measurement declarations on lines 26-34.
% - I+ gas measurement 43555: 160 mW, center [524.5297 380.8430].
% - I+ droplet measurement 43554: 160 mW, center [524.5297 380.8430].
% - I+ gas measurement 43568: 600 mW, center [524.5297 380.8430].
% - I+ droplet measurement 43567: 600 mW, center [524.5297 380.8430].
% - I+He measurement 43556: 160 mW, center [509.3664 387.6409].
% - I+He measurement 43563: 300 mW, center [509.3664 387.6409].
% - Velocity factor: vf_single = 8.6178 from v4 line 171.
%
% Requirements:
% - The legacy VMI MATLAB toolbox must be on the MATLAB path.
% - plot_processed_VMI must return r and radial_distribution fields.
%
% Outputs, written under data/reference/paper_v4/:
% - iplus_gas_160mw_43555_radial.csv: v_mps,signal_arb
% - iplus_drop_160mw_43554_radial.csv: v_mps,signal_arb
% - iplus_gas_600mw_43568_radial.csv: v_mps,signal_arb
% - iplus_drop_600mw_43567_radial.csv: v_mps,signal_arb
% - iplus_he_160mw_43556_radial.csv: v_mps,signal_arb
% - iplus_he_300mw_43563_radial.csv: v_mps,signal_arb

clear; close all;

fprintf('Exporting paper v4 radial references...\n');

global plot_processed_with_ROI
plot_processed_with_ROI = false;

vf_single = 8.6178;
center_shared = [524.5297 380.8430];
center_ihe = [509.3664 387.6409];

script_dir = fileparts(mfilename('fullpath'));
out_dir = fullfile(script_dir, '..', 'paper_v4');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

refs = {
    43555, 'iplus_gas_160mw_43555_radial.csv', center_shared;
    43554, 'iplus_drop_160mw_43554_radial.csv', center_shared;
    43568, 'iplus_gas_600mw_43568_radial.csv', center_shared;
    43567, 'iplus_drop_600mw_43567_radial.csv', center_shared;
    43556, 'iplus_he_160mw_43556_radial.csv', center_ihe;
    43563, 'iplus_he_300mw_43563_radial.csv', center_ihe;
};

for k = 1:size(refs, 1)
    measurement_id = refs{k, 1};
    filename = refs{k, 2};
    center = refs{k, 3};

    res = plot_processed_VMI(measurement_id, true, center, true);
    v_mps = res.r(:) * vf_single;
    signal_arb = res.radial_distribution(:);
    T = table(v_mps, signal_arb, ...
        'VariableNames', {'v_mps', 'signal_arb'});
    writetable(T, fullfile(out_dir, filename));
    fprintf('Saved %s\n', filename);
end
