% export_vmi_reference_data.m
clear; close all;

% NOTE: Ensure your VMI toolbox is on the MATLAB path before running this!

vf_single = 8.6178;
mass_correction_factor = sqrt(127/131);

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

% Scale arrays as done in the original plot
r_he_Aps = res_Iplus_He.r * vf_single / 100 * mass_correction_factor;
signal_he = movmean(res_Iplus_He.radial_distribution, 1); 

% Export
T_he = table(r_he_Aps(:), signal_he(:), 'VariableNames', {'v_Aps', 'signal_arb'});
writetable(T_he, 'vmi_iplus_he.csv');
fprintf('Saved vmi_iplus_he.csv\n\n');


fprintf('Processing I+ Gas Phase data...\n');
%% 2. Process Gas Phase I+
gas_measurement = 43632;
res_Iplus_gas = plot_processed_VMI(gas_measurement, 1, [482.9299 392.4866], true);
res_Iplus_gas = abel_invert_processed_VMI(res_Iplus_gas);

% Scale arrays
r_gas_Aps = res_Iplus_gas.r * vf_single / 100;
signal_gas = res_Iplus_gas.radial_distribution;

% Export
T_gas = table(r_gas_Aps(:), signal_gas(:), 'VariableNames', {'v_Aps', 'signal_arb'});
writetable(T_gas, 'vmi_iplus_gas.csv');
fprintf('Saved vmi_iplus_gas.csv\n');