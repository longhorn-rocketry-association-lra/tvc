filename = 'TVCaero.csv';
data = readtable(filename);

filtered = data(abs(data.Mach - 0.01) < 1e-6, :);

aero_map = filtered(:, {'Alpha', 'CD', 'CL', 'CP'});

aero_data.Alpha = aero_map.Alpha;
aero_data.CD = aero_map.CD;
aero_data.CL = aero_map.CL;
aero_data.CP = aero_map.CP;

save('alpha_interp.mat', 'aero_data');