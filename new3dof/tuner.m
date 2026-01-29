clear; clc;

%% ===== USER CONSTANTS =====
ell  = 0.225;          % [m] TVC moment arm
D    = 0.079;            % [m] body diameter
Lref = D;              % reference length
S    = pi*D^2/4;       % reference area
Iyy  = .007272;        % [kg*m^2] constant pitch inertia

out_file = 'pd_gainsched_thrust_constantIyy.mat';

% PIDTUNE settings
wc_min = 2;            % [rad/s]
wc_max = 8;            % [rad/s]
designFocus = 'reference-tracking';

%% ===== LOAD NOMINAL TRAJECTORY =====
load('nominaltraj.mat');  % contains 'out' struct

%% ===== Aerodynamic table (Mach -> Cm_alpha) =====
% Replace with your actual data
Mach_tbl = [0.00 0.02 0.04 0.06 0.08];
Cma_tbl  = [-0.005 -0.007 -0.009 -0.011 -0.013];

%% ===== Create thrust grid =====
Tmin = floor(min(out.T.Data)/10)*10;
Tmax = ceil(max(out.T.Data)/10)*10;

if Tmin == Tmax
    warning('Thrust range is a single value. Expanding grid slightly.');
    T_grid = Tmin + (0:1);  % 2 points for interpolation
else
    T_grid = (Tmin:10:Tmax).';
end

%% ===== Robust nearest index helper =====
function idx = nearestIndexSafe(x, xq)
    if isempty(x)
        warning('nearestIndexSafe: input array is empty. Returning 1.');
        idx = 1;
        return;
    end
    [~, idx] = min(abs(x(:) - xq)); % ensures scalar
end

%% ===== Allocate schedule =====
schedule = struct([]);
keep = false(numel(T_grid),1);

%% ===== Loop over thrust grid =====
for i = 1:numel(T_grid)
    Tq = T_grid(i);

    % Nearest trajectory point
    idx_traj = nearestIndexSafe(out.T.Data, Tq);

    % Extract operating point
    qbar_op = out.qbar.Data(idx_traj);
    T_op    = out.T.Data(idx_traj);
    Mach_op = out.Mach.Data(idx_traj);

    % Skip invalid points
    if ~isfinite(qbar_op) || ~isfinite(T_op)
        warning('Skipping T = %.2f N: invalid trajectory data', Tq);
        continue;
    end

    % Interpolate Cm_alpha
    if numel(Mach_tbl) < 2
        Cm_alpha_op = Cma_tbl(1);  % fallback if only 1 point
    else
        Cm_alpha_op = interp1(Mach_tbl, Cma_tbl, Mach_op, 'linear', 'extrap');
    end

    % Compute plant gains
    k_aero = (qbar_op * S * Lref / Iyy) * Cm_alpha_op;
    k_tvc  = (T_op * ell / Iyy);

    if ~isfinite(k_aero) || ~isfinite(k_tvc) || k_tvc <= 0
        warning('Skipping T = %.2f N: invalid gains', Tq);
        continue;
    end

    % Linear plant state-space
    A = [0 1;
         k_aero 0];
    B = [0;
         k_tvc];
    Gtheta = ss(A,B,[1 0],0);

    % Crossover frequency (constant)
    wc = (wc_min + wc_max)/2;

    opts = pidtuneOptions('DesignFocus', designFocus, 'CrossoverFrequency', wc);

    try
        Cpd = pidtune(Gtheta,'PD',opts);
    catch ME
        warning('pidtune failed at T = %.2f N: %s', Tq, ME.message);
        continue;
    end

    % Store schedule
    schedule(i).T  = Tq;
    schedule(i).Kp = Cpd.Kp;
    schedule(i).Kd = Cpd.Kd;
    keep(i) = true;
end

% Remove invalid points
schedule = schedule(keep);

%% ===== Save schedule =====
save(out_file,'schedule');
fprintf('Saved PD gain schedule over thrust (constant Iyy): %d points\n', numel(schedule));