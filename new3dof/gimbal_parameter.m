%% Gimbal actuator parameters

wn = 15;
zeta = 0.8;
max_degree = 7;

%% Transfer Function
Ys = [wn^2];
Us = [1 2*zeta*wn wn^2];