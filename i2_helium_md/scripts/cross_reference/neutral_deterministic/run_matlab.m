% run_matlab.m
% ---------------------------------------------------------------------------
% MATLAB equivalent of run_python.py.
%
% Runs a deterministic 100-step neutral-propagation simulation on a single
% I2 molecule with the bond aligned along x and starting at rest at 3.0 A
% separation. Collisions are disabled (cross section = 0). Pure leapfrog
% dynamics over the Morse pair potential + droplet solvation potential.
%
% Outputs:
%   matlab_trajectory.csv  -- with the same columns as python_trajectory.csv
%
% Run with:
%   matlab -batch "run_matlab"
% or interactively from this directory:
%   run_matlab
%
% Adjust MATLAB_LEGACY_DIR below to point to your copy of the legacy
% Iodine_Helium_Simulation/ directory.
% ---------------------------------------------------------------------------

clear; close all;



% --- globals required by the legacy code ---
global Xdip_active partner_interaction he_direction_scattering scatter_strength
Xdip_active = false;             % disable the Gaussian dip in X-state Morse
partner_interaction = true;      % include I-I Morse pair potential
he_direction_scattering = false; % don't scatter velocity directions
scatter_strength = 0;

% --- simulation parameters (must match run_python.py) ---
NUM_STEPS_TO_RUN = 100;
DT_PS = 0.01;
BOND_LENGTH_A = 3.0;
DROPLET_RADIUS_A = 27.97;
NUM_MOLECULES = 1;

% --- physical constants (MATLAB legacy values) ---
u_amu = 1.66053907e-27;       % kg
eV_J = 1.602e-19;             % J  (note: 4 sig figs, matches MATLAB)

% --- mass per atom (2N = 2 atoms, both iodine) ---
mass_amu = 127.0;
m = mass_amu * u_amu * ones(2, 1);    % column vector, kg

% --- droplet potential (steepness, binding energy in eV) ---
% These come from the production cfg single_pulse_N2000:
potential_steepness = 14.2;          % unitless (for atom-droplet)
binding_energy_I_atom = 318.43 * 1.380649e-23 / eV_J;   % K -> eV via k_B/eV
% MATLAB defines droplet_potential_atom = @(r) droplet_potential(...)
droplet_potential_atom = @(r) droplet_potential([potential_steepness, binding_energy_I_atom, 0], r);
h_grad = 0.0001;
droplet_force = @(x) (droplet_potential_atom(x + h_grad) - droplet_potential_atom(x)) / h_grad;

% --- initial geometry (atoms placed by hand; no rng) ---
half = BOND_LENGTH_A / 2.0;
x0 = [+half; -half];   % atom 1 at +half, atom 2 at -half
y0 = [0; 0];
z0 = [0; 0];
vx0 = [0; 0];
vy0 = [0; 0];
vz0 = [0; 0];
droplet_radii = DROPLET_RADIUS_A * ones(2, 1);

% --- assemble frog_step_neutral input struct ---
frog_step_crate = struct;
frog_step_crate.m = m;
frog_step_crate.droplet_radii = droplet_radii;
frog_step_crate.droplet_force = droplet_force;
frog_step_crate.he_direction_scattering = he_direction_scattering;

% --- output buffer ---
NUM_OUT_ROWS = NUM_STEPS_TO_RUN + 1;
out = zeros(NUM_OUT_ROWS, 16);
% columns:
% 1=t, 2..7=positions (x1,y1,z1,x2,y2,z2), 8..13=velocities, 14=E_kin, 15=E_pot, 16=E_total

% --- compute t=0 partner Morse via add_partner_interaction (energy only) ---
%
% NOTE: The legacy MATLAB code in vmi_sim_3d_neutral_propa_HeDFT_mimic.m
% has a bug here: line 476 sets E_pot(:,1) to droplet only, while line 885
% (subsequent steps) uses droplet + half Morse. We use the FIXED formula
% (droplet + half Morse at t=0) to match our Python port. See README.md
% "A note on the t=0 E_pot fix" for details.
[~, ~, ~, E_pot_partner_0] = add_partner_interaction(x0, y0, z0, m);
[Ek0, Ep0, Et0] = compute_energies(x0, y0, z0, vx0, vy0, vz0, m, droplet_radii, ...
                                     droplet_potential_atom, eV_J, E_pot_partner_0);
out(1, :) = [0.0, x0(1), y0(1), z0(1), x0(2), y0(2), z0(2), ...
             vx0(1), vy0(1), vz0(1), vx0(2), vy0(2), vz0(2), Ek0, Ep0, Et0];

% --- step loop ---
x = x0; y = y0; z = z0;
vx = vx0; vy = vy0; vz = vz0;

for t_idx = 1:NUM_STEPS_TO_RUN
    [x_new, y_new, z_new, vx_new, vy_new, vz_new, E_pot_partner] = ...
        frog_step_neutral(x, y, z, vx, vy, vz, frog_step_crate, DT_PS);

    [Ek, Ep, Et] = compute_energies(x_new, y_new, z_new, vx_new, vy_new, vz_new, ...
                                     m, droplet_radii, droplet_potential_atom, eV_J, ...
                                     E_pot_partner);

    t_ps = t_idx * DT_PS;
    out(t_idx + 1, :) = [t_ps, x_new(1), y_new(1), z_new(1), x_new(2), y_new(2), z_new(2), ...
                          vx_new(1), vy_new(1), vz_new(1), vx_new(2), vy_new(2), vz_new(2), ...
                          Ek, Ep, Et];

    x = x_new; y = y_new; z = z_new;
    vx = vx_new; vy = vy_new; vz = vz_new;
end

% --- write CSV ---
header_str = ['t_ps,x1_A,y1_A,z1_A,x2_A,y2_A,z2_A,', ...
              'vx1_Aps,vy1_Aps,vz1_Aps,vx2_Aps,vy2_Aps,vz2_Aps,', ...
              'E_kin_eV,E_pot_eV,E_total_eV'];
% Define where you want to save it (replace with your actual path)
save_directory = 'T:\github synchronized\Iodine_Helium_Drag_Approach\Drag_function\i2_helium_md\scripts\matlab_comparison_test'; 

% Combine the folder path and the file name safely
full_file_path = fullfile(save_directory, 'matlab_trajectory.csv');
fid = fopen(full_file_path, 'w');
fprintf(fid, '%s\n', header_str);
for ii = 1:size(out, 1)
    fprintf(fid, '%.16e', out(ii, 1));
    for k = 2:size(out, 2)
        fprintf(fid, ',%.16e', out(ii, k));
    end
    fprintf(fid, '\n');
end
fclose(fid);

% --- summary ---
fprintf('Wrote matlab_trajectory.csv with %d rows.\n', size(out, 1));
fprintf('\nFinal state:\n');
fprintf('  t = %.4f ps\n', out(end, 1));
fprintf('  atom 1 x = %.6f A   (started at %.4f)\n', out(end, 2), BOND_LENGTH_A/2);
fprintf('  atom 2 x = %.6f A   (started at %.4f)\n', out(end, 5), -BOND_LENGTH_A/2);
fprintf('  bond length = %.6f A\n', out(end, 2) - out(end, 5));
fprintf('  E_kin = %.6f eV\n', out(end, 14));
fprintf('  E_pot = %.6f eV\n', out(end, 15));
fprintf('  E_total = %.6f eV\n', out(end, 16));
fprintf('\nE_total drift: %.4f%%\n', ...
    abs(out(end, 16) - out(1, 16)) / abs(out(1, 16)) * 100);


% ---------------------------------------------------------------------------
% Local functions (must be at end of script in MATLAB R2016b+)
% ---------------------------------------------------------------------------
function [Ek, Ep, Et] = compute_energies(x, y, z, vx, vy, vz, m, droplet_radii, ...
                                          droplet_potential_atom, eV_J, E_pot_partner)
    v_sq = vx.^2 + vy.^2 + vz.^2;
    Ek_per_atom = 0.5 * m .* (v_sq * 100^2) / eV_J;   % A/ps -> m/s factor 100
    r = sqrt(x.^2 + y.^2 + z.^2);
    E_drop = droplet_potential_atom(r - droplet_radii);
    % per-pair Morse E_pot_partner is length N=1; tile to 2N=2 and split half/half
    E_partner_per_atom = [E_pot_partner; E_pot_partner] / 2.0;
    Ep_per_atom = E_drop + E_partner_per_atom;
    Ek = sum(Ek_per_atom);
    Ep = sum(Ep_per_atom);
    Et = Ek + Ep;
end
