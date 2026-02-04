% --- Aero CSV -> MAT processing (adds Cm_alpha) ---
clc
clear 

filename = 'spaceshot.csv';
data = readtable(filename, "VariableNamingRule","modify");

% ---- CONSTANTS (edit these) ----
cgx_m    = 4.2683;   % CG from nose datum [m]
d_ref_m  = 0.40;     % reference length for Cm (usually body diameter) [m]
cp_units = "in";     % "in" if CP column is inches, "m" if already meters
out_mat  = 'aero.mat';
% --------------------------------

% Inspect actual MATLAB column names (useful for debugging)
vnames = string(data.Properties.VariableNames);

% Robustly find CNalpha column (handles trailing underscores, etc.)
cn_candidates = vnames(contains(lower(vnames), "cnalpha") & contains(lower(vnames), "perrad"));
if isempty(cn_candidates)
    cn_candidates = vnames(contains(lower(vnames), "cnalpha"));
end
if isempty(cn_candidates)
    error("Could not find CNalpha column. Available columns:\n%s", strjoin(vnames, ", "));
end
CNalpha_name = cn_candidates(1);

% Pull required columns
Mach    = data{:, 'Mach'};
Alpha   = data{:, 'Alpha'};
CNalpha = data{:, CNalpha_name};
CP_raw  = data{:, 'CP'};

% Optional (keep if you want them in the MAT)
hasCD = any(vnames == "CD");
hasCL = any(vnames == "CL");
if hasCD, CD = data{:, 'CD'}; end
if hasCL, CL = data{:, 'CL'}; end

% Convert CP to meters
switch lower(cp_units)
    case "in"
        CP_m = CP_raw * 0.0254;
    case "m"
        CP_m = CP_raw;
    otherwise
        error("cp_units must be ""in"" or ""m"".");
end

% Compute Cm_alpha about CG:
%   Cm_alpha = -((CP - CG)/d_ref) * CNalpha
rx_m     = CP_m - cgx_m;
Cm_alpha = -(rx_m ./ d_ref_m) .* CNalpha;

% Build aero_data struct
aero_data.Mach     = Mach;
aero_data.Alpha    = Alpha;
if hasCD, aero_data.CD = CD; end
if hasCL, aero_data.CL = CL; end
aero_data.CP       = CP_m;        % [m] from nose
aero_data.CNalpha  = CNalpha;     % [1/rad]
aero_data.Cm_alpha = Cm_alpha;    % [1/rad]

save(out_mat, 'aero_data');
