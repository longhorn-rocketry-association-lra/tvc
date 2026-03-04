load('pd_gainsched_thrust_constantIyy.mat','schedule');

T_array  = [schedule.T];   % numeric 1xN
Kp_array = [schedule.Kp];  % numeric 1xN
Kd_array = [schedule.Kd];  % numeric 1xN

T_array = [0; 8; 10; 24];
Kp_array = [0.15, 0.1, 0.9, 0.8];
Kd_array = [0.05, 0.05, 0.05, 0.05];
save('Scheduledata.mat', 'T_array', 'Kp_array', 'Kd_array');