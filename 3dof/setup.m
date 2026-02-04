mass = 0.5 %kg
mdot = poly(1); %polyfit fuel mdot




motor = table2array(readtable('motor.txt'));

thrust_poly = poly(1); %polyfit thrust


earth_radius = 6378e3; %m 
earth_mass = 5.9722 * 10 ^ (24); %kg
G = 6.6743 * 10 ^ (-11); 


save("params.mat");

data = readmatrix('i55_thrust.csv');
t_data = data(:,1); % Time values
T_data = data(:,2); % Thrust values

% Save variables to a .mat file
save('thrust_data.mat', 't_data', 'T_data');
disp('Thrust data saved to thrust_data.mat');