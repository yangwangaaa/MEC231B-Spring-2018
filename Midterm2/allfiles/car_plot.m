function car_plot(u,z0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the motion of a car according to the model:
% {x}_{k+1} &=& x_k-R_k*sin(\theta_k)+R_k*sin(\theta_k+\beta_k) 
%{y}_{k+1} &=& y_k+R_k*cos(\theta_k)-R_k*cos(\theta_k+\beta_k) 
%{\theta}_{k+1} &=& \theta_k+\beta_k 
%---------------------------------------------------
% with inputs u=[R_k,\beta_k] and state z=[z,y,\theta]
% u is  2xN matrix where u(:,k) are the two inputs at time k
% z0 is a 3x1 vector of initial conditions.
%--------------------------------------------
% In the plot, red wheels are the  front  vehicle wheels

% Sampling parameter
rep=40;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of a simple kinematic vehicle
t=1;
xsim=[];
usim=[];
xsim(:,1)=z0;
N=size(u,2);
for k = 1:N        
    for j = 1:rep        
         t=t+1;
         xsim(1,t) =xsim(1,t-1) -u(1,k)*sin(xsim(3,t-1))+u(1,k)*sin(xsim(3,t-1)+1/rep*u(2,k));
         xsim(2,t) =xsim(2,t-1) +u(1,k)*cos(xsim(3,t-1))-u(1,k)*cos(xsim(3,t-1)+1/rep*u(2,k));
         xsim(3,t) = xsim(3,t-1) + 1/rep*u(2,k);
         usim(:,t-1)=[u(1,k);1/rep*u(2,k)];
    end
end

usim(:,t)=usim(:,t-1);
X=xsim';
U=usim';

% 1.4) Car Model 
auto.w = 2.0;                % car width [m]
auto.db = 1.2;               % rear axis position, from back [m]
auto.df = 1.0;               % front axis position, from front [m]
auto.l = 2.8+auto.df+auto.db % car length [m]
auto.tyr = 0.8;              % tyre diameter [m]
auto.dmax = 25*pi/180;       % maximum front wheel steering angle [rad]
auto.drat = 14.5;            % ratio between steering wheel angle and front wheel angle
auto.d = auto.l - auto.df - auto.db;  % axel distance [m] %L in the midterm

%Compute axis limits here
xmin=min(xsim(1,:));
xmax=max(xsim(1,:));
ymin=min(xsim(2,:));
ymax=max(xsim(2,:));
L=max(xmax-xmin,ymax-ymin);
limits=[(xmin+xmax)/2-L/2-2,(xmin+xmax)/2+L/2+2,(ymin+ymax)/2-L/2-2,(ymin+ymax)/2+L/2+2];

fig = figure()
hold on
axis(limits)
xlabel('x [m]')
ylabel('y [m]')
title('Navigation')

% Trajectory
plottraj.linestyle = '-';
plottraj.linewidth = 1;
plottraj.color = [0 0 1]; % RGB value
plottraj.marker = 'none';

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Vehicle trajectory (oversampled, to show the path)
for i = 1:N*rep+1
    p = plotcar(X(i,1),X(i,2),X(i,3),U(i,2),auto,fig,[0.3 0.3 0.3]);
    pause(0.002)
    if i<=N*rep & i>1 
        delete(p)
    end
    if i<N*rep
        plot([X(i,1),X(i+1,1)],[X(i,2),X(i+1,2)],plottraj);
    end
end
for i = 1:rep:N*rep+1
 plot(X(i,1),X(i,2),'go')
 [u,v] = pol2cart(X(i,3),1);
 quiver(X(i,1),X(i,2),u,v,'g','linewidth',2,'maxheadsize',2);
end

%%%%%%
% Plot model states


%% Attribution
% ME C231A and EECS C220B, UC Berkeley, Fall 2016