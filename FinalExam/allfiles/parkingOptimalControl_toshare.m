clear all
clear yalmip

%% Car Parameters
lr  = 1.7;
lf = 1.1;
L=lf+lr;

%% Ego car G and g, not used in this code
width=1; %car width 
G = [-1 0; 1 0;0 -1; 0 1]; g = [0; lr+lf; width; width]; % polyhedron Gy<=g as in the paper

%% Input constraints
umin(2)=-pi; 
delta_max=25*pi/180; % Max steering
umin(1)=(lf+lr)/tan(delta_max);
model.u.min=umin; 

%% State constraints 
model.z.min=-[20;10;2*pi]; 
model.z.max=[20;20;2*pi]; 

%Initial, terminal conditions and horizon
z0 = [-10;10;0];  
zT = [0;0;-pi/2];  
N=3; %navigation  horizon


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Obstacle list
i=1;
obs{i}.center=[-5;0];
obs{i}.LW=[2;8];
obs{i}.theta=0/180; %(in radiants)
i=i+1;
obs{i}.center=[5;0];
obs{i}.LW=[2;8];
obs{i}.theta=0/180; %(in radiants)
i=i+1;
obs{i}.center=[0;14];
obs{i}.LW=[1;8];
obs{i}.theta=0/180; %(in radiants)


%% Some obtacle postprocessing 
for j=1:length(obs)
    t=obs{j}.theta;
    % generate T matrix for each obstacle
    obs{j}.T=[cos(t), -sin(t);sin(t) cos(t)]*diag(obs{j}.LW/2);
    % polyehdral representaion
    obs{j}.poly=obs{j}.T*unitbox(2)+obs{j}.center;
    [AA{j},bb{j}]=double(obs{j}.poly);
end

%% Setup the Navigation Problem
%options = sdpsettings('solver','ipopt');
options = sdpsettings('solver','fmincon','verbose',1);

z = sdpvar(3,N+1);
u = sdpvar(2,N);
constr = [z(:,N+1)==zT, z(:,1) == z0];
cost = 0;
for k = 1:N
constr = constr+[z(1,k+1) == z(1,k)-u(1,k)*sin(z(3,k))+u(1,k)*sin(z(3,k)+u(2,k)),...    
                  z(2,k+1) == z(2,k)+u(1,k)*cos(z(3,k))-u(1,k)*cos(z(3,k)+u(2,k)),...
                  z(3,k+1) == z(3,k)+u(2,k),...
                  model.u.min(1) <= u(1,k),model.u.min(2) <= u(2,k),u(2,k) <= -model.u.min(2),...
                  model.z.min <= z(:,k+1),z(:,k+1)<=model.z.max];
                  cost = cost + (u(2,k)*u(1,k))^2;
end

%% Compute Navigation Solution
optimize(constr,cost,options)
zdata = double(z);
udata = double(u);

%% Plot Solution
figure
plot(zdata(1,:),zdata(2,:),'r')
hold on
car_plot(double(u),z0)
%% Plot obstacles
for j=1:length(obs)
plot(polytope(AA{j},bb{j}));
end
axis([model.z.min(1) model.z.max(1) model.z.min(2) model.z.max(2)])
