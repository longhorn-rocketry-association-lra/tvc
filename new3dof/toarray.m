load('pd_gainsched_thrust_constantIyy.mat','schedule');

T_array  = [schedule.T];   % numeric 1xN
Kp_array = [schedule.Kp];  % numeric 1xN
Kd_array = [schedule.Kd];  % numeric 1xN