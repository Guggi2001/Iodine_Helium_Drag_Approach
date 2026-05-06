% export_matlab_multistep.m
%
% Inline legacy-equivalent MATLAB reference for the multi-step
% deterministic ion-propagation cross-reference (CLAUDE.md validation
% targets 3 and 4: several deterministic ion steps with collisions
% disabled, plus energy bookkeeping).
%
% This script reimplements the relevant pieces of frog_step_ion.m,
% add_partner_interaction_ion.m, ion_interaction_potential.m, and
% droplet_potential.m INLINE so it has no dependency on the legacy
% repository's path or globals. The numerical update order, finite-
% difference force forms, and unit-conversion literals are taken
% verbatim from the legacy code so any constant-rounding effect
% remains observable in the comparison.
%
% Stochastic ion physics is disabled at the configuration level:
% cross-section = 0 (no collisions), mass_attach_probability = 0
% (no helium attachment), sigma_dependent_on_v = false,
% ion_scatter_angle_std_deg = 0, additional_droplet_charges = 0.
%
% --- t=0 bookkeeping -----------------------------------------------------
% The legacy MATLAB ion script has two known t=0 bookkeeping bugs:
%   * vmi_sim_3d_ion_propa.m:289 -- E_kin uses (vx^2+vy^2)^2 (missing vz,
%     squared); we use the CORRECTED form here so the comparison can
%     focus on the propagation logic instead of the t=0 bugs.
%   * vmi_sim_3d_ion_propa.m:291 -- E_pot uses 2D radial and omits the
%     partner Coulomb; we use the CORRECTED form (3D radial + half-
%     pair Coulomb), again to match the Python intentional fix.
%
% Bug-by-bug t=0 verification was already covered by the prior task
% (scripts/cross_reference/ion_t0_state).
%
% Run:
%   matlab -batch "cd('scripts/cross_reference/ion_multistep_no_collision'); export_matlab_multistep"
%
% Output: matlab_multistep.csv (one row per (step, atom) pair).

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

% --- Read shared inputs --------------------------------------------------
inputs = jsondecode(fileread(fullfile(script_dir, 'inputs.json')));

% --- Legacy MATLAB physical constants (verbatim from physical_constants.m)
eV = 1.602e-19;            % J         -- rounded vs CODATA
u  = 1.66053907e-27;       % kg        -- rounded vs CODATA

% --- Numerical inputs ----------------------------------------------------
num_molecules  = inputs.num_molecules;
num_particles  = 2 * num_molecules;
num_steps      = inputs.num_steps;
dt             = inputs.dt_ion_ps;
mass_per_atom  = inputs.mass_amu * u;
R              = inputs.droplet_radius_A;
S              = inputs.cfg_flags.potential_steepness;
B              = inputs.cfg_flags.binding_energy_I_ion_eV;
E_coul_scale   = inputs.cfg_flags.E_coulomb_scale;

x  = [inputs.atom_0.x_A;    inputs.atom_1.x_A];
y  = [inputs.atom_0.y_A;    inputs.atom_1.y_A];
z  = [inputs.atom_0.z_A;    inputs.atom_1.z_A];
vx = [inputs.atom_0.vx_Aps; inputs.atom_1.vx_Aps];
vy = [inputs.atom_0.vy_Aps; inputs.atom_1.vy_Aps];
vz = [inputs.atom_0.vz_Aps; inputs.atom_1.vz_Aps];

mass_i        = repmat(mass_per_atom, num_particles, 1);
droplet_radii = repmat(R,             num_particles, 1);
charge_i      = ones(num_particles, 1);   % all +1, single ionization OFF

% --- Inline droplet potential and finite-difference force ---------------
% droplet_potential.m + frog_step_ion.m line 105
h_drop = 1e-6;
drop_pot   = @(r) ((erf((r - 0)/S) + 1)/2) * B;
drop_force = @(r) (drop_pot(r + h_drop) - drop_pot(r))/h_drop;     % +dU/dr (eV/A)

% --- Inline ion-ion partner Coulomb potential and FD force --------------
% ion_interaction_potential.m + add_partner_interaction_ion.m
h_pair = 1e-4;
ion_pair_U = @(r) E_coul_scale * 1.0 * 1.0 * 14.39964548 ./ r;     % eV at separation r (q1=q2=1)

% --- Output buffer -------------------------------------------------------
% One row per (step, atom). Columns:
%   step, t_ps, atom, x_A, y_A, z_A, vx_Aps, vy_Aps, vz_Aps,
%   mass_kg, E_kin_eV, E_pot_eV, E_dissip_eV
out_rows = cell(0, 13);

% --- t=0 energy bookkeeping (CORRECTED — see header) --------------------
[E_drop_t0, E_part_t0_per_atom] = deal_energies(x, y, z, droplet_radii, drop_pot, ion_pair_U);
v_sq_t0 = vx.^2 + vy.^2 + vz.^2;
E_kin_t0    = mass_i .* (v_sq_t0 .* (100.^2)) / 2 / eV;
E_pot_t0    = E_drop_t0 + E_part_t0_per_atom;
E_dissip_t0 = zeros(num_particles, 1);
out_rows = append_step(out_rows, 0, 0.0, x, y, z, vx, vy, vz, mass_i, E_kin_t0, E_pot_t0, E_dissip_t0);

% --- Time integration: velocity-Verlet, droplet + partner Coulomb -------
t_ps = 0.0;
for step = 1:num_steps
    % accelerations at current positions
    [ax0, ay0, az0] = ion_acceleration( ...
        x, y, z, mass_i, droplet_radii, charge_i, ...
        drop_force, ion_pair_U, h_pair, eV, u);

    % drift to new positions
    x1 = x + dt*vx + 0.5*ax0*dt^2;
    y1 = y + dt*vy + 0.5*ay0*dt^2;
    z1 = z + dt*vz + 0.5*az0*dt^2;

    % accelerations at new positions
    [ax1, ay1, az1] = ion_acceleration( ...
        x1, y1, z1, mass_i, droplet_radii, charge_i, ...
        drop_force, ion_pair_U, h_pair, eV, u);

    % kick velocities
    vx1 = vx + 0.5*(ax0 + ax1)*dt;
    vy1 = vy + 0.5*(ay0 + ay1)*dt;
    vz1 = vz + 0.5*(az0 + az1)*dt;

    % advance scalars/state
    x = x1; y = y1; z = z1;
    vx = vx1; vy = vy1; vz = vz1;
    t_ps = t_ps + dt;

    % per-step energies (matches vmi_sim_3d_ion_propa.m lines 761, 765)
    [E_drop, E_part_per_atom] = deal_energies(x, y, z, droplet_radii, drop_pot, ion_pair_U);
    v_sq = vx.^2 + vy.^2 + vz.^2;
    E_kin    = mass_i .* (v_sq .* (100.^2)) / 2 / eV;
    E_pot    = E_drop + E_part_per_atom;
    E_dissip = zeros(num_particles, 1);     % stays 0: collisions disabled

    out_rows = append_step(out_rows, step, t_ps, x, y, z, vx, vy, vz, mass_i, E_kin, E_pot, E_dissip);
end

% --- Write CSV -----------------------------------------------------------
out_csv = fullfile(script_dir, 'matlab_multistep.csv');
fid = fopen(out_csv, 'w');
fprintf(fid, ['step,t_ps,atom,', ...
              'x_A,y_A,z_A,', ...
              'vx_Aps,vy_Aps,vz_Aps,', ...
              'mass_kg,E_kin_eV,E_pot_eV,E_dissip_eV\n']);
for i = 1:size(out_rows, 1)
    r = out_rows(i, :);
    fprintf(fid, '%d,%.16e,%d,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n', ...
        r{1}, r{2}, r{3}, r{4}, r{5}, r{6}, r{7}, r{8}, r{9}, r{10}, r{11}, r{12}, r{13});
end
fclose(fid);

fprintf('Wrote %s with %d rows (%d steps x %d atoms).\n', ...
    out_csv, size(out_rows, 1), num_steps + 1, num_particles);


% =========================================================================
% Local helper functions
% =========================================================================
function [E_drop, E_part_per_atom] = deal_energies(x, y, z, droplet_radii, drop_pot, ion_pair_U)
    % Per-atom droplet + half-pair Coulomb energies (eV).
    % Mirrors line 765 of vmi_sim_3d_ion_propa.m and the t=0 corrected
    % form. r is 3D (correct), and partner Coulomb is split half/half
    % between the two atoms of each molecule.
    n = size(x, 1);
    half = n/2;
    r_atom = sqrt(x.^2 + y.^2 + z.^2);
    E_drop = drop_pot(r_atom - droplet_radii);

    % pair separation, per molecule (length N)
    dx = x(1:half) - x(half+1:end);
    dy = y(1:half) - y(half+1:end);
    dz = z(1:half) - z(half+1:end);
    dr = sqrt(dx.^2 + dy.^2 + dz.^2);

    U_pair = ion_pair_U(dr);                         % per pair (length N)
    E_part_per_atom = [U_pair; U_pair] / 2.0;        % per atom (length 2N)
end

function [ax, ay, az] = ion_acceleration(x, y, z, mass, droplet_radii, charge, ...
                                         drop_force, ion_pair_U, h_pair, eV, u)
    % Combined acceleration: droplet + partner Coulomb.
    n = size(x, 1);
    half = n/2;

    % --- droplet contribution (frog_step_ion.m lines 32-51, 117-137) ---
    r0_atom = sqrt(x.^2 + y.^2 + z.^2);
    depth   = r0_atom - droplet_radii;
    F_drop  = drop_force(depth) * 1.602e-9;          % eV/A -> N (legacy literal)
    a_drop  = -F_drop ./ mass;                       % m/s^2
    a_drop  = a_drop * 1e-14;                        % -> A/ps^2 (legacy literal)
    rx = x ./ r0_atom; ry = y ./ r0_atom; rz = z ./ r0_atom;
    ax = a_drop .* rx;
    ay = a_drop .* ry;
    az = a_drop .* rz;

    % --- partner Coulomb contribution (add_partner_interaction_ion.m) -
    dx = x(1:half) - x(half+1:end);
    dy = y(1:half) - y(half+1:end);
    dz = z(1:half) - z(half+1:end);
    dr = sqrt(dx.^2 + dy.^2 + dz.^2);
    dr_unit = [dx ./ dr, dy ./ dr, dz ./ dr];        % length N x 3

    F_pair = (ion_pair_U(dr) - ion_pair_U(dr + h_pair)) / h_pair;   % eV/A
    F_full = [F_pair; F_pair];                                       % length 2N
    a_pair = F_full ./ (mass / u) * 9648.53322;                      % legacy literal

    a_vec  = a_pair .* [dr_unit; -dr_unit];          % atom 1: +; atom 2: -
    ax = ax + a_vec(:, 1);
    ay = ay + a_vec(:, 2);
    az = az + a_vec(:, 3);
end

function rows = append_step(rows, step, t_ps, x, y, z, vx, vy, vz, mass, E_kin, E_pot, E_dissip)
    n = size(x, 1);
    for atom = 1:n
        rows(end+1, :) = { ...
            step, t_ps, atom-1, ...
            x(atom), y(atom), z(atom), ...
            vx(atom), vy(atom), vz(atom), ...
            mass(atom), E_kin(atom), E_pot(atom), E_dissip(atom) ...
        };
    end
end
