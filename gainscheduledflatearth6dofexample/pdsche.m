%% Gain-scheduled PD (pidtune) vs Mach (step = 0.1)
clear; clc;

load("nominaltraj.mat");   % provides 'out'
load("aero.mat");          % provides 'aero_data'

%% ===== USER CONSTANTS =====
ell  = 3.3317;           % [m] TVC moment arm
D    = 0.4;              % [m] body diameter
S    = pi*(D^2)/4;       % [m^2]
Lref = D;                % [m] reference length used in Cm_alpha definition

out_file = "pd_gainsched.mat";
% ===========================

%% ===== Pull trajectory signals =====
tM = out.mach.Time(:);
M  = out.mach.Data(:);

tq   = out.q_star.Time(:);
qbar = out.q_star.Data(:);

tT = out.thrust.Time(:);
T  = out.thrust.Data(:);

tI = out.Iyy.Time(:);
Idata = out.Iyy.Data;

% Sample time (still computed, but only saved if you want it)
dt = diff(tM);
Ts = median(dt(~isnan(dt) & isfinite(dt) & dt > 0));
if isempty(Ts) || Ts <= 0
    error("Could not determine a valid sample time Ts from out.mach.Time.");
end

%% ===== Build Mach -> Cm_alpha map (average over Alpha) =====
Mach_tbl = aero_data.Mach(:);
Cma_tbl  = aero_data.Cm_alpha(:);

[Mach_u, ~, g] = unique(Mach_tbl);
Cma_u = accumarray(g, Cma_tbl, [], @mean);

%% ===== Mach grid (step 0.1) over trajectory =====
Mmin = floor(min(M)/0.1)*0.1;
Mmax = ceil(max(M)/0.1)*0.1;
Mach_grid = (Mmin:0.1:Mmax).';

%% ===== Extract Iyy(t) vector robustly =====
sz = size(Idata);
if numel(sz) == 3 && all(sz(1:2) == [3 3])
    Iyy_traj = squeeze(Idata(2,2,:));
elseif numel(sz) == 3 && all(sz(2:3) == [3 3])
    Iyy_traj = squeeze(Idata(:,2,2));
else
    error('Unexpected inertia Data size. Got size = [%s].', num2str(sz));
end
Iyy_traj = Iyy_traj(:);

%% ===== PIDTUNE SETTINGS =====
wc_min = 2;    % [rad/s]
wc_max = 8;    % [rad/s]
designFocus = "reference-tracking";   % or "disturbance-rejection"

%% ===== Allocate output schedule (ONLY what you need) =====
schedule = struct([]);
schedule(numel(Mach_grid)).Mach_grid = []; % preallocate struct array

keep = false(numel(Mach_grid),1);

for i = 1:numel(Mach_grid)
    Mq = Mach_grid(i);

    % nearest trajectory operating point at this Mach
    idx_traj = nearestIndex(M, Mq);
    t_op = tM(idx_traj);
    M_op = M(idx_traj);

    % aligned qbar, thrust, Iyy
    qbar_op = qbar(nearestIndex(tq, t_op));
    T_op    = T(nearestIndex(tT, t_op));
    Iyy_op  = Iyy_traj(nearestIndex(tI, t_op));

    if ~isfinite(qbar_op) || ~isfinite(T_op) || ~isfinite(Iyy_op) || Iyy_op <= 0
        continue;
    end

    % Cm_alpha from aero table (nearest Mach)
    idx_aero = nearestIndex(Mach_u, M_op);
    Cm_alpha_op = Cma_u(idx_aero);

    % linear model terms
    k_aero = (qbar_op * S * Lref / Iyy_op) * Cm_alpha_op;
    k_tvc  = (T_op * ell / Iyy_op);

    if ~isfinite(k_aero) || ~isfinite(k_tvc) || k_tvc <= 0
        continue;
    end

    % continuous plant u -> theta
    A = [0      1;
         k_aero 0];
    B = [0;
         k_tvc];
    Gtheta = ss(A, B, [1 0], 0);

    % schedule crossover frequency wc linearly across Mach range
    frac = (Mq - Mach_grid(1)) / max(eps, (Mach_grid(end) - Mach_grid(1)));
    wc = wc_min + (wc_max - wc_min) * frac;
    wc = max(min(wc, wc_max), wc_min);

    % model-based PD autotune
    opts = pidtuneOptions('DesignFocus', designFocus, 'CrossoverFrequency', wc);

    try
        Cpd = pidtune(Gtheta, "PD", opts);
    catch
        continue;
    end

    % store ONLY essentials
    schedule(i).Mach_grid = Mq;
    schedule(i).Kp = Cpd.Kp;
    schedule(i).Kd = Cpd.Kd;
    schedule(i).Ts = Ts;     % optional; remove if you truly don't want it

    keep(i) = true;
end

% drop failed points
schedule = schedule(keep);

if isempty(schedule)
    error("No valid schedule points produced. Check signals and units.");
end

save(out_file, "schedule");   % only schedule saved

fprintf("Saved PD gain schedule: %s\n", out_file);
fprintf("Points: %d, Ts = %.6g s, Mach range %.2f to %.2f (step 0.1)\n", ...
        numel(schedule), Ts, Mach_grid(1), Mach_grid(end));

%% ===== Local helper =====
function idx = nearestIndex(x, x_query)
    [~, idx] = min(abs(x - x_query));
end
