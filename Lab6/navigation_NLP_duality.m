clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Francesco Borrelli ME C231A 2015
% Kinematic Navigation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=20;
sampling=10;
%Var Defintions
z = sdpvar(2,N);

%Initial and terminal condition
z0 = [0;0];
zT = [10;10];
dzmin=-[2;2]/2;
dzmax=[2;2]/2;

%Obstacle list
i=1;
obs{i}.center=[4;4];
obs{i}.LW=[2;2];
obs{i}.theta=30*pi/180; %(in radiants)
i=i+1;
obs{i}.center=[5.5;7];
obs{i}.LW=[4;1];
obs{i}.theta=-45*pi/180; %(in radiants)
i=i+1;
obs{i}.center=[9;4];
obs{i}.LW=[2;2.5];
obs{i}.theta=-30*pi/180; %(in radiants)

% some obtacle postprocessing 
for j=1:length(obs)
    t=obs{j}.theta;
    % generate T matrix for each obstacle
    obs{j}.T=[cos(t), -sin(t);sin(t) cos(t)]*diag(obs{j}.LW/2);
    % polyehdral representaion
    obs{j}.poly=obs{j}.T*unitbox(2)+obs{j}.center;
    [AA{j},bb{j}]=double(obs{j}.poly);
    lambda{j} = sdpvar(size(AA{j},1),N,'full');
end


%try to remove/add this one


%Constraints
%Setup Optimization Problem
cost = 0;
constr = [z(:,1)==z0;z(:,N)==zT];
Q=eye(2);
%constr = [zmin<=z(:,N)<= zmax, z(:,1)==z0,z(:,N)==zT];
for t = 2:N
      cost=cost+(z(:,t)-z(:,t-1))'*Q*(z(:,t)-z(:,t-1));
      constr = constr +[dzmin<= z(:,t)-z(:,t-1)<=dzmax]; 
          for j=1:length(obs) 
              for k = 0:sampling-1
                xs=z(:,t-1)+k/sampling*(z(:,t)-z(:,t-1));
                %constr = constr +[(xs-obs{j}.center)'*inv(obs{j}.T)'*inv(obs{j}.T)*(xs-obs{j}.center)>=2]; 
                constr = constr +[(AA{j}*xs-bb{j})'*lambda{j}(:,t)>=0.3];% try 0.3, 0.1 restoration failed, try adding penetration
              end
              constr = constr +[lambda{j}(:,t)'*AA{j}*AA{j}'*lambda{j}(:,t)==1]; %here norm fails
              constr = constr +[lambda{j}(:,t)>=0]; 
          end
end
options = sdpsettings('solver','ipopt');
%options.ipopt=ipoptset('linear_solver','MUMPS');
solvesdp(constr,cost,options);
z_vec = double(z);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Functions % to add title and labels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
th = 0:pi/50:2*pi;
for j=1:length(obs)
    for l=1:length(th)
        z=[cos(th(l));sin(th(l))]*sqrt(2);
        y=obs{j}.T*z+obs{j}.center;
        xobs{j}(l) = y(1);
        yobs{j}(l) = y(2);
    end
end

%% plot routine
figure
plot(z_vec(1,:),z_vec(2,:),'o')
hold on
for j=1:length(obs)
%plot(obs{j}.T*unitbox(2)+obs{j}.center);
plot(polytope(AA{j},bb{j}));
%plot(xobs{j}, yobs{j},'b');
end
plot(z_vec(1,:),z_vec(2,:),'o')
