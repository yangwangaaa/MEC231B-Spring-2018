function p2()
%% lab5: vehicle dynamics p2.m
% discrete dynamics
close all; clc;

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

% states
%x=[vx,vy,wz,X,Y,psi]
tspan = [0 20];

%% Lane Change
x0 = [20; 0; 0; 0;0;0];
Fx=0; %N
delta_Lane = @(t) deg2rad(5*(t>5 & t<=7)-5*(t>12 & t<=14)); % rads

odefun = @(t,x) bicycle_discrete_nonlin_dyn(t,x,[Fx; delta_Lane(t)],params);
odefun2 = @(t,x) bicycle_model_nonlin_dyn(t,x,[Fx; delta_Lane(t)],params);

% continuous
[t_cont, x_cont] = ode45(odefun2,tspan,x0);

deltaT = 0.01;
t_dist = 0:deltaT:20;
x_dist = x0';
for t = t_dist
xk = odefun(t,x0);
x_dist = [x_dist; xk'];
x0 = xk;
end
x_dist(end,:) = [];

drawPlots();

%%
function drawPlots()
figure;
subplot(1,2,1); hold on;
plot(x_cont(:,4),x_cont(:,5),'r','linewidth',1);
scatter(x0(4),x0(5),80,'x'); hold on;
plot(x_dist(:,4),x_dist(:,5),':b','linewidth',2);
grid on; title('X vs Y');xlabel('X'); ylabel('Y');
legend('xy path contin','initial state','xy path discrete','location','best');

subplot(1,2,2);hold on;
plot(t_cont,x_cont(:,1),'g',t_cont,x_cont(:,2),'m','linewidth',2); ...
hold on;
plot(t_dist,x_dist(:,1),':r',t_dist,x_dist(:,2),':b','linewidth',2);
grid on; legend('vx cont','vy cont','vx dist','vy dist','location','best');
title('vx and vy');xlabel('time (s)'); ylabel('m/s');
end
end