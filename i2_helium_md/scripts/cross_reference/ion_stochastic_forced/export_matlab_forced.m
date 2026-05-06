% export_matlab_forced.m
%
% Inline legacy-equivalent MATLAB reference for the stochastic forced-
% event ion-driver cross-reference (CLAUDE.md validation target 5).
%
% Reimplements the relevant pieces of vmi_sim_3d_ion_propa.m
% (lines ~300-770) and frog_step_ion.m INLINE so the script has no
% dependency on the legacy repo path or globals. Uses the legacy-rounded
% physical constants and the literal force-conversion factors from the
% legacy code so the constants-rounding effect remains observable.
%
% Forced-event design: cross-section sigma_0 = 1e6 A^2 with
% sigma_dependent_on_v = true and exponent = -2 makes the Mode-3
% probability dr * sigma * rho_droplet >> 1, so trial < p_scatter is
% always satisfied and a collision fires every step (after the first,
% which sets the previous-step distance). mass_attach_probability = 1
% makes every collision an attachment. Event counts and mass history
% are deterministic on both sides; post-collision velocities still
% depend on the RNG stream and are NOT compared cross-language.
%
% Bookkeeping note: this script implements the corrected t=0 energies
% (3D radial + half-pair Coulomb + full vz). The t=0 legacy bugs are
% out of scope here -- they were validated by ion_t0_state.
%
% Run:
%   matlab -batch "cd('scripts/cross_reference/ion_stochastic_forced'); export_matlab_forced"
%
% Output: matlab_forced.csv (one row per (step, atom)).

script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir)
    script_dir = pwd;
end

% --- Read shared inputs --------------------------------------------------
inputs = jsondecode(fileread(fullfile(script_dir, 'inputs.json')));

% --- Legacy MATLAB physical constants (verbatim from physical_constants.m)
eV = 1.602e-19;
u  = 1.66053907e-27;

% --- Numerical inputs ----------------------------------------------------
num_molecules = inputs.num_molecules;
num_particles = 2 * num_molecules;
num_steps     = inputs.num_steps;
dt            = inputs.dt_ion_ps;
mass_per_atom = inputs.mass_amu * u;
R             = inputs.droplet_radius_A;
S             = inputs.cfg_flags.potential_steepness;
B             = inputs.cfg_flags.binding_energy_I_ion_eV;
E_coul_scale  = inputs.cfg_flags.E_coulomb_scale;

sigma0_const                   = inputs.cfg_flags.geometric_scattering_crosssection_Iplus;
sigma_dependent_on_v           = inputs.cfg_flags.sigma_dependent_on_v;
sigma_ion_exponent             = inputs.cfg_flags.sigma_ion_exponent;
mass_attach_probability        = inputs.cfg_flags.mass_attach_probability;
scatter_mass_ion               = inputs.cfg_flags.scatter_mass_ion_amu;
ion_scatter_angle_std_deg      = inputs.cfg_flags.ion_scatter_angle_std_deg;
v_limit_m_per_s                = inputs.cfg_flags.v_limit_m_per_s;

% Landau cutoff E_min (matches MATLAB run_simulation.m):
%   E_min = (127 u) * v_limit^2 / 2 / eV  with v_limit in A/ps
v_limit_Aps = v_limit_m_per_s / 100.0;
E_min = (127.0 * u) * (v_limit_Aps * 100.0)^2 / 2 / eV;

% Helium droplet density inside the droplet (legacy convention).
bulk_density_helium = 0.0219;
density_droplet     = 0.8 * bulk_density_helium;

x  = [inputs.atom_0.x_A;    inputs.atom_1.x_A];
y  = [inputs.atom_0.y_A;    inputs.atom_1.y_A];
z  = [inputs.atom_0.z_A;    inputs.atom_1.z_A];
vx = [inputs.atom_0.vx_Aps; inputs.atom_1.vx_Aps];
vy = [inputs.atom_0.vy_Aps; inputs.atom_1.vy_Aps];
vz = [inputs.atom_0.vz_Aps; inputs.atom_1.vz_Aps];

mass_i        = repmat(mass_per_atom, num_particles, 1);
droplet_radii = repmat(R,             num_particles, 1);
charge_i      = ones(num_particles, 1);

% --- Inline droplet potential and finite-difference force ---------------
h_drop = 1e-6;
drop_pot   = @(r) ((erf((r - 0)/S) + 1)/2) * B;
drop_force = @(r) (drop_pot(r + h_drop) - drop_pot(r))/h_drop;

% --- Inline ion-ion partner Coulomb potential ---------------------------
h_pair = 1e-4;
ion_pair_U = @(r) E_coul_scale * 1.0 * 1.0 * 14.39964548 ./ r;

% --- Output buffer --------------------------------------------------------
% Columns: step, t_ps, atom, x, y, z, vx, vy, vz, mass_kg,
%          E_kin, E_pot, E_dissip, E_mass_attach_defect,
%          number_of_collisions, b_collision, b_attach, sigma_used, depth
out_rows = cell(0, 19);

% --- t=0 corrected bookkeeping -----------------------------------------
[E_drop_t0, E_part_t0_per_atom] = deal_energies(x, y, z, droplet_radii, drop_pot, ion_pair_U);
v_sq_t0 = vx.^2 + vy.^2 + vz.^2;
E_kin     = mass_i .* (v_sq_t0 .* (100.^2)) / 2 / eV;
E_pot     = E_drop_t0 + E_part_t0_per_atom;
E_dissip  = zeros(num_particles, 1);
E_defect  = zeros(num_particles, 1);
n_coll    = zeros(num_particles, 1);

% sigma_used at t=0 (recorded for cross-language analytic check)
v0_speed = sqrt(v_sq_t0);
if sigma_dependent_on_v
    sigma_used_t0 = sigma0_const * v0_speed.^sigma_ion_exponent;
else
    sigma_used_t0 = repmat(sigma0_const, num_particles, 1);
end
depth_t0 = sqrt(x.^2 + y.^2 + z.^2) - droplet_radii;

out_rows = append_step(out_rows, 0, 0.0, x, y, z, vx, vy, vz, mass_i, ...
    E_kin, E_pot, E_dissip, E_defect, n_coll, ...
    zeros(num_particles,1), zeros(num_particles,1), sigma_used_t0, depth_t0);

% --- Time integration with collisions/attachments -----------------------
prev_dist = nan(num_particles, 1);   % NaN -> first step has no collisions
t_ps = 0.0;

for step = 1:num_steps
    % ---- pre-step bookkeeping for sigma_used (records the value the
    %      driver feeds to sample_collision_events: sigma(v_pre)) -----
    v_pre_speed = sqrt(vx.^2 + vy.^2 + vz.^2);
    if sigma_dependent_on_v
        sigma_used = sigma0_const * v_pre_speed.^sigma_ion_exponent;
    else
        sigma_used = repmat(sigma0_const, num_particles, 1);
    end

    % ---- velocity-Verlet leapfrog (frog_step_ion equivalent) ----------
    [ax0, ay0, az0] = ion_acceleration(x, y, z, mass_i, droplet_radii, ...
        charge_i, drop_force, ion_pair_U, h_pair, eV, u);
    x1 = x + dt*vx + 0.5*ax0*dt^2;
    y1 = y + dt*vy + 0.5*ay0*dt^2;
    z1 = z + dt*vz + 0.5*az0*dt^2;
    [ax1, ay1, az1] = ion_acceleration(x1, y1, z1, mass_i, droplet_radii, ...
        charge_i, drop_force, ion_pair_U, h_pair, eV, u);
    vx1 = vx + 0.5*(ax0 + ax1)*dt;
    vy1 = vy + 0.5*(ay0 + ay1)*dt;
    vz1 = vz + 0.5*(az0 + az1)*dt;

    % ---- depth at new position --------------------------------------
    r1    = sqrt(x1.^2 + y1.^2 + z1.^2);
    depth = r1 - droplet_radii;

    % ---- post-leapfrog speed and E0 (matches MATLAB line 384) -------
    v_post_leap_sq    = vx1.^2 + vy1.^2 + vz1.^2;
    v_post_leap_speed = sqrt(v_post_leap_sq);
    E0 = mass_i .* v_post_leap_sq * (100^2) / 2 / eV;

    % ---- v-dependent cross section AT post-leapfrog speed (matches
    %      vmi_sim_3d_ion_propa.m:444 which uses the post-leapfrog v).
    if sigma_dependent_on_v
        sigma = sigma0_const * v_post_leap_speed.^sigma_ion_exponent;
    else
        sigma = repmat(sigma0_const, num_particles, 1);
    end

    % ---- Mode-3 collision sampling ---------------------------------
    if any(isnan(prev_dist))
        b_collision = false(num_particles, 1);
    else
        trial_random_number = rand(num_particles, 1);
        p_scatter = prev_dist .* sigma .* density_droplet;
        b_collision = (trial_random_number < p_scatter) & (depth < 0);
        % Landau cutoff
        b_collision = b_collision & ~(E0 < E_min);
    end

    % ---- collision kinematics (only if any collisions) -------------
    [vx_after, vy_after, vz_after, dE] = apply_collision_inline( ...
        vx1, vy1, vz1, mass_i, b_collision, ...
        scatter_mass_ion, ion_scatter_angle_std_deg, u, eV);

    % ---- mass attachment (after collision) -------------------------
    mass_attach_trial = rand(num_particles, 1);
    b_attach = (mass_attach_trial < mass_attach_probability) & b_collision;
    mass_old = mass_i;
    mass_i   = mass_i + b_attach * 4 * u;
    mass_diff = mass_i - mass_old;

    % ---- accumulators (E_kin uses NEW mass; defect uses post-collision v)
    v_post_sq    = vx_after.^2 + vy_after.^2 + vz_after.^2;
    E_kin    = mass_i .* v_post_sq * (100^2) / 2 / eV;

    [E_drop_step, E_part_step] = deal_energies(x1, y1, z1, droplet_radii, drop_pot, ion_pair_U);
    E_pot    = E_drop_step + E_part_step;
    E_dissip = E_dissip + dE;
    E_defect = E_defect - 0.5 * mass_diff .* v_post_sq * (100^2) / eV;
    n_coll   = n_coll + double(b_collision);

    % ---- per-atom step distance for next-step collision sampler ----
    prev_dist = sqrt((x1 - x).^2 + (y1 - y).^2 + (z1 - z).^2);

    % ---- advance state ---------------------------------------------
    x  = x1;  y  = y1;  z  = z1;
    vx = vx_after; vy = vy_after; vz = vz_after;
    t_ps = t_ps + dt;

    out_rows = append_step(out_rows, step, t_ps, x, y, z, vx, vy, vz, mass_i, ...
        E_kin, E_pot, E_dissip, E_defect, n_coll, ...
        double(b_collision), double(b_attach), sigma_used, depth);
end

% --- Write CSV -----------------------------------------------------------
out_csv = fullfile(script_dir, 'matlab_forced.csv');
fid = fopen(out_csv, 'w');
fprintf(fid, ['step,t_ps,atom,', ...
              'x_A,y_A,z_A,vx_Aps,vy_Aps,vz_Aps,', ...
              'mass_kg,E_kin_eV,E_pot_eV,E_dissip_eV,E_mass_attach_defect_eV,', ...
              'number_of_collisions,b_collision,b_attach,sigma_used_A2,depth_A\n']);
for i = 1:size(out_rows, 1)
    r = out_rows(i, :);
    fprintf(fid, ['%d,%.16e,%d,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e,', ...
                  '%.16e,%.16e,%.16e,%.16e,%.16e,', ...
                  '%d,%d,%d,%.16e,%.16e\n'], ...
        r{1}, r{2}, r{3}, r{4}, r{5}, r{6}, r{7}, r{8}, r{9}, ...
        r{10}, r{11}, r{12}, r{13}, r{14}, ...
        r{15}, r{16}, r{17}, r{18}, r{19});
end
fclose(fid);

fprintf('Wrote %s with %d rows (%d steps x %d atoms).\n', ...
    out_csv, size(out_rows, 1), num_steps + 1, num_particles);


% =========================================================================
% Local helpers (same as multistep script + collision kinematics)
% =========================================================================
function [E_drop, E_part_per_atom] = deal_energies(x, y, z, droplet_radii, drop_pot, ion_pair_U)
    n = size(x, 1);
    half = n/2;
    r_atom = sqrt(x.^2 + y.^2 + z.^2);
    E_drop = drop_pot(r_atom - droplet_radii);
    dx = x(1:half) - x(half+1:end);
    dy = y(1:half) - y(half+1:end);
    dz = z(1:half) - z(half+1:end);
    dr = sqrt(dx.^2 + dy.^2 + dz.^2);
    U_pair = ion_pair_U(dr);
    E_part_per_atom = [U_pair; U_pair] / 2.0;
end

function [ax, ay, az] = ion_acceleration(x, y, z, mass, droplet_radii, charge, ...
                                         drop_force, ion_pair_U, h_pair, eV, u)
    n = size(x, 1);
    half = n/2;
    r0_atom = sqrt(x.^2 + y.^2 + z.^2);
    depth = r0_atom - droplet_radii;
    F_drop = drop_force(depth) * 1.602e-9;
    a_drop = -F_drop ./ mass;
    a_drop = a_drop * 1e-14;
    rx = x ./ r0_atom; ry = y ./ r0_atom; rz = z ./ r0_atom;
    ax = a_drop .* rx;
    ay = a_drop .* ry;
    az = a_drop .* rz;

    dx = x(1:half) - x(half+1:end);
    dy = y(1:half) - y(half+1:end);
    dz = z(1:half) - z(half+1:end);
    dr = sqrt(dx.^2 + dy.^2 + dz.^2);
    dr_unit = [dx ./ dr, dy ./ dr, dz ./ dr];
    F_pair = (ion_pair_U(dr) - ion_pair_U(dr + h_pair)) / h_pair;
    F_full = [F_pair; F_pair];
    a_pair = F_full ./ (mass / u) * 9648.53322;
    a_vec = a_pair .* [dr_unit; -dr_unit];
    ax = ax + a_vec(:, 1);
    ay = ay + a_vec(:, 2);
    az = az + a_vec(:, 3);
end

function [vx_new, vy_new, vz_new, dE] = apply_collision_inline( ...
        vx, vy, vz, mass, b_collision, scatter_mass_ion, scatter_angle_std_deg, u, eV)
    % Inline reimplementation of the lab-frame elastic-scattering kinematics
    % from vmi_sim_3d_ion_propa.m:507-678. RHO is the per-atom mass ratio.
    n = size(vx, 1);
    v_speed = sqrt(vx.^2 + vy.^2 + vz.^2);
    safe = v_speed > 0;
    v_unit_x = zeros(n, 1); v_unit_y = zeros(n, 1); v_unit_z = zeros(n, 1);
    v_unit_x(safe) = vx(safe) ./ v_speed(safe);
    v_unit_y(safe) = vy(safe) ./ v_speed(safe);
    v_unit_z(safe) = vz(safe) ./ v_speed(safe);
    v_unit_x(~safe) = 1.0;

    E0 = 0.5 * mass .* (v_speed * 100).^2 / eV;

    % Impact parameter b/R via inverse CDF of 2 b/R^2 -> b/R = sqrt(u).
    impact_parameter_norm = sqrt(rand(n, 1));
    COSTHETA = 2 * impact_parameter_norm.^2 - 1;
    SINTHETA = sqrt(max(0, 1 - COSTHETA.^2));
    COSTHETA(~b_collision) = 1;
    SINTHETA(~b_collision) = 0;

    RHO = (mass / u) / scatter_mass_ion;
    E1  = E0 .* (1 + 2*RHO.*COSTHETA + RHO.^2) ./ (1 + RHO).^2;
    dE  = E0 - E1;
    dE(~b_collision) = 0;

    arg = 1 + 2*RHO.*COSTHETA + RHO.^2;
    arg = max(arg, 0);
    den = sqrt(arg);
    safe_den = den > 0;
    COStheta_lab = ones(n, 1);
    COStheta_lab(safe_den) = (COSTHETA(safe_den) + RHO(safe_den)) ./ den(safe_den);
    COStheta_lab = max(-1, min(1, COStheta_lab));
    SINtheta_lab = sqrt(max(0, 1 - COStheta_lab.^2));

    if scatter_angle_std_deg > 0
        idx = find(b_collision);
        if ~isempty(idx)
            std_rad = scatter_angle_std_deg * pi / 180;
            theta = acos(COStheta_lab(idx)) + randn(numel(idx), 1) * std_rad;
            COStheta_lab(idx) = cos(theta);
            SINtheta_lab(idx) = sqrt(max(0, 1 - COStheta_lab(idx).^2));
        end
    end

    % Random orthonormal basis in plane perpendicular to v.
    ref = rand(n, 3) - 0.5;
    ref_norm = sqrt(sum(ref.^2, 2));
    ref_norm(ref_norm == 0) = 1;
    ref = ref ./ ref_norm;

    v_unit = [v_unit_x, v_unit_y, v_unit_z];
    n1 = cross(v_unit, ref);
    n1_norm = sqrt(sum(n1.^2, 2));
    n1_norm(n1_norm == 0) = 1;
    n1 = n1 ./ n1_norm;
    n2 = cross(v_unit, n1);

    % MATLAB legacy azimuth convention.
    COSBETA = (rand(n, 1) - 0.5) * 2;
    SINBETA = sqrt(max(0, 1 - COSBETA.^2));

    v_new_speed = sqrt(2 * E1 * eV ./ mass) / 100;
    v_par  = v_new_speed .* COStheta_lab;
    v_perp = v_new_speed .* SINtheta_lab;

    new_v = v_unit .* v_par + n1 .* (COSBETA .* v_perp) + n2 .* (SINBETA .* v_perp);
    new_v(~b_collision, 1) = vx(~b_collision);
    new_v(~b_collision, 2) = vy(~b_collision);
    new_v(~b_collision, 3) = vz(~b_collision);

    vx_new = new_v(:, 1);
    vy_new = new_v(:, 2);
    vz_new = new_v(:, 3);
end

function rows = append_step(rows, step, t_ps, x, y, z, vx, vy, vz, mass, ...
                            E_kin, E_pot, E_dissip, E_defect, n_coll, ...
                            b_collision, b_attach, sigma_used, depth)
    n = size(x, 1);
    for atom = 1:n
        rows(end+1, :) = { ...
            step, t_ps, atom-1, ...
            x(atom), y(atom), z(atom), ...
            vx(atom), vy(atom), vz(atom), ...
            mass(atom), E_kin(atom), E_pot(atom), E_dissip(atom), E_defect(atom), ...
            n_coll(atom), b_collision(atom), b_attach(atom), sigma_used(atom), depth(atom) ...
        };
    end
end
