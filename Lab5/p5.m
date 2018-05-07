%% lab5: vehicle dynamics p5.m
% extended kalman filter using the linearized/nonlinear dynamics
clear; close all; clc;

%% parameters
% vehicle params
params.m = 2237; % kg
params.Jz = 5112; % kg*m^2
params.a = 1.46; % m
params.b = 1.55; % m

% tire params
params.A = -6.8357;
params.B = 0.0325;
params.C = 238.9874;
params.dt = .01;
% linearization point
params.q0 = [10; 0; 0; 0; 0; 0];

% probability information
params.m0 = [3; zeros(5,1)];
params.Sigma0 = diag([0.1*ones(6,1)]);

params.E = [eye(6), zeros(6,3)];

params.F = [zeros(3,6), eye(3)];
params.W = diag([1,1,0.1, 1, 1, 0.1, 1, 0.1, 0.05]);
%%
VehicleData = load('VehicleData.mat'); % Load linearized dynamics
linDyn = load('bicycle_lindyn.mat'); %%KF
N = length(VehicleData.time);
% initialize
xEst = zeros(6,N+1);
yEst = zeros(3,N);
xVar = zeros(6,6,N+1);
xEst(:,1) = params.m0;
xVar(:,:,1) = params.Sigma0;


for i = 1:N
    fprintf("iteration %d\n",i);
    uk = [VehicleData.Fx_in(i); VehicleData.delta_in(i)];
    yk = [VehicleData.vx_meas(i);...
        VehicleData.ay_meas(i)+...
        VehicleData.vx_meas(i)*VehicleData.yawRate_meas(i);...
        VehicleData.yawRate_meas(i)];
    xk = xEst(:,i);
    p = [params.m; params.Jz; params.a; params.b;...
        params.A ;params.B; params.C; params.dt]; [Ak,~,Ck,~] = getLinearDynamics_at_k(xk,uk,p);
    [xk1k,Sxk1k, ykk1] = myEKF231B(xk,xVar(:,:,i),...
        [],[],Ck,params.E,params.F,params.W,uk,yk,params);
    xEst(:,i+1) = xk1k;
    yEst(:,i) = ykk1;
    xVar(:,:,i+1) = Sxk1k;
end
betaEst = atan(xEst(2,:)./xEst(1,:));
ayEst = yEst(2,:)-xEst(1,1:end-1).*xEst(3,1:end-1);

%% Plotting
figure; hold on;
plot([0 VehicleData.time], betaEst,'r','linewidth',1);
plot([0 VehicleData.time],[0 VehicleData.beta_meas'],':b','linewidth',2);
grid on; grid minor;
title('\beta_{Est} vs \beta_{meas}');
legend('\beta_{Est}','\beta_{meas}');
xlabel('time (s)');
figure; hold on;
plot([VehicleData.time], ayEst,'r','linewidth',2);
plot([VehicleData.time],[VehicleData.ay_meas'],':b','linewidth',2);
grid on; grid minor;
title('ay_{Est} vs ay_{meas}');
legend('ay_{Est}','ay_{meas}');
xlabel('time (s)');

figure; hold on;
plot([VehicleData.time], yEst(3,:),'r','linewidth',1);
plot([VehicleData.time],[VehicleData.yawRate_meas'],':b','linewidth',2);
grid on; grid minor;
title('\omegaz_{Est} vs \omegaz_{meas}');
legend('\omegaz_{Est}','\omegaz_{meas}');
xlabel('time (s)');

figure; hold on;
plot([0 VehicleData.time], xEst(1,:),'r','linewidth',1);
plot([VehicleData.time],[VehicleData.vx_meas'],':b','linewidth',2);
grid on; grid minor;
title('vx_{Est} vs vx_{meas}');
legend('vx_{Est}','vx_{meas}');
xlabel('time (s)');
