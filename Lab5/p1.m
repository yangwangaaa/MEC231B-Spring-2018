function p1()
%% lab5: vehicle dynamics p1.m
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
%% 1.a constant velocity skidpad 21
x0 = [20; 0; 0; 0; 0; 0];
Fx=2500; %N
delta = deg2rad(10); % rads
odefun = @(t,x) bicycle_model_nonlin_dyn(t,x,[Fx; delta],params);
[t, x] = ode45(odefun,tspan,x0);
drawPlots(t,x,Fx*ones(size(t)),delta*ones(size(t)));

%% 1.b constant acceleration skidpad
x0 = [5; 0; 0; 0; 0; 0];
Fx=8000; %N
delta = deg2rad(10); % rads
odefun = @(t,x) bicycle_model_nonlin_dyn(t,x,[Fx; delta],params);
[t, x] = ode45(odefun,tspan,x0);
drawPlots(t,x,Fx*ones(size(t)),delta*ones(size(t)));

%% 1.c lane change
x0 = [20; 0; 0; 0; 0; 0];
Fx=0; %N
delta_Lane = @(t) deg2rad(5*(t>5 & t<=7)-5*(t>12 & t<=14)); % rads
odefun = @(t,x) bicycle_model_nonlin_dyn(t,x,[Fx; delta_Lane(t)],params);
[t, x] = ode45(odefun,tspan,x0);
drawPlots(t,x,Fx*ones(size(t)),delta_Lane(t));

%% Plotting
function drawPlots(t,x,Fx,delta)
figure;
subplot(2,2,1);
plot(t,Fx,'r','linewidth',2);
grid on; title('Fx');
xlabel('time (s)'); ylabel('N');

subplot(2,2,2); hold on;
plot(x(:,4),x(:,5),'linewidth',2);
scatter(x(1,4),x(1,5),80,'x');
grid on; title('X vs Y');xlabel('X'); ylabel('Y');
legend('xy path','initial state','location','best');

subplot(2,2,3);
plot(t,delta,'r','linewidth',2);
grid on; title('\∆');
xlabel('time (s)'); ylabel('degrees');

subplot(2,2,4);
plot(t,x(:,1),'r',t,x(:,2),'b','linewidth',2);
grid on; legend('vx','vy','location','best');

end
end