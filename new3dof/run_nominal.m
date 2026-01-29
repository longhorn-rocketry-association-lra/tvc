init_rocket
simOut = sim('simu');
logsout = simOut.logsout;

theta = logsout.get('theta').Values;
q     = logsout.get('q').Values;
Mach  = logsout.get('Mach').Values;
T     = logsout.get('Thrust').Values;
Iyy   = logsout.get('Iyy').Values;
qbar  = logsout.get('qbar').Values;


out.theta = theta;
out.qbar  = qbar;
out.q     = q;
out.Mach  = Mach;
out.T     = T;
out.Iyy   = Iyy;
save('nominaltraj.mat','out');

