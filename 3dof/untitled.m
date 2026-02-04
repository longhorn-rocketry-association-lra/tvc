data = readmatrix('F10_thrust.csv');  
t_data = data(:,1);  % Time values from CSV
T_data = data(:,2);  % Thrust values from CSV

% Generate interpolation points
t_interp = linspace(min(t_data), max(t_data), 100);  

% Perform PCHIP interpolation
T_pchip = interp1(t_data, T_data, t_interp, 'pchip'); 

% Plot results
figure;
plot(t_data, T_data, 'ro', 'MarkerFaceColor', 'r', 'DisplayName', 'Original Data'); % Original data points
hold on;
plot(t_interp, T_pchip, 'b-', 'LineWidth', 1.5, 'DisplayName', 'PCHIP Interpolation'); % Interpolated curve
hold off;

% Labels and title
xlabel('Time (s)');
ylabel('Thrust (N)');
title('Thrust Curve Using PCHIP Interpolation');
legend('Location', 'Best');
grid on;

save('thrust_data1.mat', 't_data', 'T_data');