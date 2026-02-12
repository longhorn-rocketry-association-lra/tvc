%% Gimbal actuator parameters

wn = 60;       % natural frequency [rad/s]
zeta = 0.8;    % damping ratio
max_degree = 7;

%% Transfer Function
Ys = [wn^2];
Us = [1 2*zeta*wn wn^2];
 