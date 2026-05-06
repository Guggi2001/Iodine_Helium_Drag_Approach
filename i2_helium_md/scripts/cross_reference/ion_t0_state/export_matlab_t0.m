% export_matlab_t0.m
%
% Self-contained inline reimplementation of the ion-init block of
% vmi_sim_3d_ion_propa.m (lines 133-294), with all stochastic ion
% physics disabled and NO propagation step taken. Reads inputs.json
% and writes matlab_t0.csv.
%
% This script deliberately does NOT call run_simulation, inputfile_*,
% or any function from legacy_matlab_repository/. The legacy ion-init
% formulas at t=0 are reproduced inline -- including the documented
% bugs -- so that the resulting CSV reflects exactly what the legacy
% code produces.
%
% Run:
%     matlab -batch "cd('scripts/cross_reference/ion_t0_state'); export_matlab_t0"
% or, from the MATLAB GUI with the working directory at this script:
%     export_matlab_t0
%
% Output: matlab_t0.csv (12 rows, columns: quantity, atom_0, atom_1).

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

% --- Read shared inputs --------------------------------------------------
inputs = jsondecode(fileread(fullfile(script_dir, 'inputs.json')));

% --- MATLAB legacy physical constants (verbatim from physical_constants.m)
eV = 1.602e-19;          % Joule          -- rounded relative to CODATA
u  = 1.66053907e-27;     % kg             -- rounded relative to CODATA

% --- Hand-crafted neutral end-state, 2N layout (N=1) ---------------------
num_molecules = inputs.num_molecules;
num_particles = 2 * num_molecules;

x_ci  = [inputs.atom_0.x_A;    inputs.atom_1.x_A];
y_ci  = [inputs.atom_0.y_A;    inputs.atom_1.y_A];
z_ci  = [inputs.atom_0.z_A;    inputs.atom_1.z_A];
vx_ci = [inputs.atom_0.vx_Aps; inputs.atom_1.vx_Aps];
vy_ci = [inputs.atom_0.vy_Aps; inputs.atom_1.vy_Aps];
vz_ci = [inputs.atom_0.vz_Aps; inputs.atom_1.vz_Aps];

mass_i        = repmat(inputs.mass_amu * u,         num_particles, 1);
droplet_radii = repmat(inputs.droplet_radius_A,     num_particles, 1);

% Single-charge ionization is disabled -> all atoms carry +1.
charge_i = ones(num_particles, 1);

time_i_t0 = inputs.time_ps_at_t0;

% --- Droplet ion potential (vmi_sim_3d_ion_propa.m line 92, 101) ---------
%   beta = [potential_steepness, binding_energy, offset]
droplet_potential = @(beta, x) ((erf((x - beta(3))/beta(1))*1 + 1)/2) * beta(2);
beta_drop = [
    inputs.cfg_flags.potential_steepness, ...
    inputs.cfg_flags.binding_energy_I_ion_eV, ...
    0.0 ...
];
droplet_potential_ion = @(r) droplet_potential(beta_drop, r);

% --- Legacy t=0 bookkeeping (vmi_sim_3d_ion_propa.m lines 284-294) -------
% Reproduced VERBATIM, including the two known bugs:
%   * E_kin uses (vx^2 + vy^2)^2 instead of vx^2 + vy^2 + vz^2,
%   * E_pot uses sqrt(x^2 + y^2) (2D) instead of 3D, and omits the
%     partner-Coulomb contribution.
% These rows are tagged INTENTIONAL_FIX_* by the comparison script.
E_kin_ion_t0    = mass_i .* (vx_ci.^2 + vy_ci.^2).^2 / 2 / eV;
E_pot_ion_t0    = droplet_potential_ion( sqrt(x_ci.^2 + y_ci.^2) - droplet_radii );
E_dissip_ion_t0 = zeros(num_particles, 1);

% No propagation step is taken.

% --- Write CSV (long-format: quantity, atom_0, atom_1) -------------------
rows = {
    'x_A',              x_ci(1),            x_ci(2);
    'y_A',              y_ci(1),            y_ci(2);
    'z_A',              z_ci(1),            z_ci(2);
    'vx_Aps',           vx_ci(1),           vx_ci(2);
    'vy_Aps',           vy_ci(1),           vy_ci(2);
    'vz_Aps',           vz_ci(1),           vz_ci(2);
    'mass_kg',          mass_i(1),          mass_i(2);
    'droplet_radius_A', droplet_radii(1),   droplet_radii(2);
    'time_ps',          time_i_t0,          time_i_t0;
    'E_kin_eV',         E_kin_ion_t0(1),    E_kin_ion_t0(2);
    'E_pot_eV',         E_pot_ion_t0(1),    E_pot_ion_t0(2);
    'E_dissip_eV',      E_dissip_ion_t0(1), E_dissip_ion_t0(2);
};

out_csv = fullfile(script_dir, 'matlab_t0.csv');
fid = fopen(out_csv, 'w');
fprintf(fid, 'quantity,atom_0,atom_1\n');
for i = 1:size(rows, 1)
    fprintf(fid, '%s,%.16e,%.16e\n', rows{i,1}, rows{i,2}, rows{i,3});
end
fclose(fid);

fprintf('Wrote %s with %d quantities.\n', out_csv, size(rows,1));
for i = 1:size(rows, 1)
    fprintf('  %-18s  atom_0 = %+.6e   atom_1 = %+.6e\n', ...
        rows{i,1}, rows{i,2}, rows{i,3});
end
