clear; clc;

%% ===== USER CONSTANTS =====
ell  = 0.225;          % [m] TVC moment arm
Iyy  = 0.007272;       % [kg*m^2] constant pitch inertia

out_file = 'pd_gainsched_thrust_constantIyy.mat';

% PIDTUNE settings
wc_min = 2;            % [rad/s]
wc_max = 8;            % [rad/s]
designFocus = 'reference-tracking';

%% ===== LOAD NOMINAL TRAJECTORY =====
load('nominaltraj.mat');  % contains 'out' struct

%% ===== Create thrust grid =====
Tmin = floor(min(out.T.Data)/10)*10;
Tmax = ceil(max(out.T.Data)/10)*10;

if Tmin == Tmax
    warning('Thrust range is a single value. Expanding grid slightly.');
    T_grid = Tmin + (0:1);  % 2 points for interpolation
else
    T_grid = linspace(Tmin, Tmax, 5).';
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
    T_op = out.T.Data(idx_traj);

    % Skip invalid points
    if ~isfinite(T_op)
        warning('Skipping T = %.2f N: invalid trajectory data', Tq);
        continue;
    end

    % ===== Compute TVC gain only =====
    k_tvc = (T_op * ell / Iyy);

    if ~isfinite(k_tvc) || k_tvc <= 0
        warning('Skipping T = %.2f N: invalid TVC gain', Tq);
        continue;
    end

    % ===== Double Integrator Plant =====
    %   theta_ddot = k_tvc * delta
    %
    %   x1 = theta
    %   x2 = theta_dot

    A = [0 1;
         0 0];

    B = [0;
         k_tvc];

    Gtheta = ss(A,B,[1 0],0);

    % ===== Crossover frequency =====
    wc = (wc_min + wc_max)/2;

    opts = pidtuneOptions('DesignFocus', designFocus, ...
                          'CrossoverFrequency', wc);

    try
        Cpd = pidtune(Gtheta,'PD',opts);
    catch ME
        warning('pidtune failed at T = %.2f N: %s', Tq, ME.message);
        continue;
    end

    % ===== Store schedule =====
    schedule(i).T  = Tq;
    schedule(i).Kp = Cpd.Kp;
    schedule(i).Kd = Cpd.Kd;
    keep(i) = true;
end

% Remove invalid points
schedule = schedule(keep);

%% ===== Save schedule =====
save(out_file,'schedule');
fprintf('Saved PD gain schedule over thrust (double integrator model): %d points\n', numel(schedule));