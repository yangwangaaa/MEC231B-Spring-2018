%% Final Exam - ME C231B Jun Zeng
% We also include all seperate files in the './allfiles/' folder for testing

addpath('./allfiles/');

%% Vehicle Navigation

%% 1.a
% Nothing to hand in

%% 1.b
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
N=4; %navigation  horizon


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
    lambda{j} = sdpvar(size(AA{j},1),N,'full');
end

%% Setup the Navigation Problem
%options = sdpsettings('solver','ipopt');
options = sdpsettings('solver','fmincon','verbose',1);

z = sdpvar(3,N+1);
u = sdpvar(2,N);
constr = [z(:,N+1)==zT, z(:,1) == z0];
cost = 0;
SampleNum = 10;
for k = 1:N
    constr = constr+...
         [z(1,k+1) == z(1,k)-u(1,k)*sin(z(3,k))+u(1,k)*sin(z(3,k)+u(2,k)),...
          z(2,k+1) == z(2,k)+u(1,k)*cos(z(3,k))-u(1,k)*cos(z(3,k)+u(2,k)),...
          z(3,k+1) == z(3,k)+u(2,k),...
          model.u.min(1) <= u(1,k),model.u.min(2) <= u(2,k),u(2,k) <= -model.u.min(2),...
          model.z.min <= z(:,k+1),z(:,k+1)<=model.z.max];
    cost = cost + (u(2,k)*u(1,k))^2;
    for p = 1:SampleNum-1
        for q = 1:size(obs,2)
            zs = z(:,k)+p/SampleNum*(z(:,k+1)-z(:,k));
            A = AA{q}; b = bb{q};
            constr = constr + [(AA{q}*zs(1:2)-bb{q})'*lambda{q}(:,k) > 0];
            constr = constr + [lambda{q}(:,k)'*AA{q}*AA{q}'*lambda{q}(:,k)<=1];
            constr = constr + [lambda{q}(:,k) >= 0];
        end
    end
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

%% 1.c
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
N=6; %navigation  horizon


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
    lambda{j} = sdpvar(size(AA{j},1),N,'full');
end

%% Setup the Navigation Problem
options = sdpsettings('solver','ipopt');
%options = sdpsettings('solver','fmincon','verbose',1);

z = sdpvar(3,N+1);
u = sdpvar(2,N);
constr = [z(:,N+1)==zT, z(:,1) == z0];
cost = 0;
SampleNum = 10;
for k = 1:N
    constr = constr+...
         [z(1,k+1) == z(1,k)-u(1,k)*sin(z(3,k))+u(1,k)*sin(z(3,k)+u(2,k)),...
          z(2,k+1) == z(2,k)+u(1,k)*cos(z(3,k))-u(1,k)*cos(z(3,k)+u(2,k)),...
          z(3,k+1) == z(3,k)+u(2,k),...
          model.u.min(1) <= u(1,k),model.u.min(2) <= u(2,k),u(2,k) <= -model.u.min(2),...
          model.z.min <= z(:,k+1),z(:,k+1)<=model.z.max];
    cost = cost + (u(2,k)*u(1,k))^2;
    for p = 1:SampleNum-1
        for q = 1:size(obs,2)
            zs = z(:,k)+p/SampleNum*(z(:,k+1)-z(:,k));
            A = AA{q}; b = bb{q};
            constr = constr + [(AA{q}*zs(1:2)-bb{q})'*lambda{q}(:,k) > 0];
            constr = constr + [lambda{q}(:,k)'*AA{q}*AA{q}'*lambda{q}(:,k)<=1];
            constr = constr + [lambda{q}(:,k) >= 0];
        end
    end
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

%% 1.d
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
N=6; %navigation  horizon


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
    lambda{j} = sdpvar(size(AA{j},1),N,'full');
    s{j} = sdpvar(size(bb{j},1),N,'full');
end

%% Setup the Navigation Problem
options = sdpsettings('solver','ipopt');
%options = sdpsettings('solver','fmincon','verbose',1);

z = sdpvar(3,N+1);
u = sdpvar(2,N);
constr = [z(:,N+1)==zT, z(:,1) == z0];
cost = 0;
SampleNum = 10;
for k = 1:N
    constr = constr+...
         [z(1,k+1) == z(1,k)-u(1,k)*sin(z(3,k))+u(1,k)*sin(z(3,k)+u(2,k)),...
          z(2,k+1) == z(2,k)+u(1,k)*cos(z(3,k))-u(1,k)*cos(z(3,k)+u(2,k)),...
          z(3,k+1) == z(3,k)+u(2,k),...
          model.u.min(1) <= u(1,k),model.u.min(2) <= u(2,k),u(2,k) <= -model.u.min(2),...
          model.z.min <= z(:,k+1),z(:,k+1)<=model.z.max];
    cost = cost + (u(2,k)*u(1,k))^2;
    for p = 1:SampleNum-1
        for q = 1:size(obs,2)
            zs = z(:,k)+p/SampleNum*(z(:,k+1)-z(:,k));
            A = AA{q}; b = bb{q};
            constr = constr + [(AA{q}*zs(1:2)-bb{q})'*lambda{q}(:,k) >= -s{q}(k)];
            % avoid unnecessary warning of inequality
            constr = constr + [lambda{q}(:,k)'*AA{q}*AA{q}'*lambda{q}(:,k)==1];
            constr = constr + [lambda{q}(:,k) >= 0];
            constr = constr + [s{q}(k) >= 0];
        end
    end
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

%% 1.e
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
N=4; %navigation  horizon


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
    lambda{j} = sdpvar(size(AA{j},1),N,'full');
end

%% Setup the Navigation Problem
options = sdpsettings('solver','ipopt');
%options = sdpsettings('solver','fmincon','verbose',1);

z = sdpvar(3,N+1);
u = sdpvar(2,N);
constr = [z(:,N+1)==zT, z(:,1) == z0];
cost = 0;
SampleNum = 10;
slack = 10000;
for k = 1:N
    constr = constr+...
         [z(1,k+1) == z(1,k)-u(1,k)*sin(z(3,k))+u(1,k)*sin(z(3,k)+u(2,k)),...
          z(2,k+1) == z(2,k)+u(1,k)*cos(z(3,k))-u(1,k)*cos(z(3,k)+u(2,k)),...
          z(3,k+1) == z(3,k)+u(2,k),...
          model.u.min(2) <= u(2,k),u(2,k) <= -model.u.min(2),...
          model.z.min <= z(:,k+1),z(:,k+1)<=model.z.max]; %model.u.min(1) <= u(1,k),... don't need anymore
    cost = cost + (u(2,k)*u(1,k))^2;
    % Here we introduce a slack variable for the radius of turning which we
    % have seen in ME231A to solve this problem
    cost = cost + slack*(1/(u(1,k)^2)-1/(model.u.min(1)^2));
    for p = 1:SampleNum-1
        for q = 1:size(obs,2)
            zs = z(:,k)+p/SampleNum*(z(:,k+1)-z(:,k));
            A = AA{q}; b = bb{q};
            constr = constr + [(AA{q}*zs(1:2)-bb{q})'*lambda{q}(:,k) > 0];
            constr = constr + [lambda{q}(:,k)'*AA{q}*AA{q}'*lambda{q}(:,k)<=1];
            constr = constr + [lambda{q}(:,k) >= 0];
        end
    end
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

%% Robust Control

%% 2
% At the line 3, we use norm calculator with a setting of accuracy to get
% the maximum frequency response and its corresponded frequency, thus in
% the line 4, we have A = abs(freqresp(G,B)). At the line 5, we use frd to
% compare the system G with a newly defined system, where at the frequency B,
% the frequency response is exactly the maximum frequency response of G.
% Finally, we see the superposition of bodegram of two these system at
% frequnecy B.

q2;

%% 3
% At the line 4, we use random complex number delta to generate D. We go to
% the file cnum2sys.m, as we have seen in the homework and the slide, the
% behavior of cnumsys.m is that the infinity norm of this system is the
% same as the complex number generated and the associated frequency is
% wBar.  At the line 5, we calculate the frequency 
% response of D and to compare it with delta (absolutely the same!). At the
% line 6, we verify that wBar is indeed the frequency. At the line 7, we
% generate the bode plot of D, we find out the gain is constant and equals
% to 20*log(delta).
q3;

%% q.4

%% 4.a
P = tf(1,[1,-1]);
C = tf([5.8 9],[0.04 1 0]);
isstable(feedback(P,C))

%% 4.b
G = -C/(1+P*C);

%% 4.c
normTol = 0.001;
fprintf('The smallest absolute value of Delta calculated by the small gain thm\n')
[InfNorm,freq] = norm(G,inf,normTol);
bound = 1/InfNorm
%verification: to get the smallest delta
M = freqresp(G,freq);
[U,S,V]=svd(M);
deltamin = -1/S(1,1)*V(:,1)*U(:,1)';
pole(feedback(P-deltamin,C))
%one pure imaginary pool appear at this critical value 

%% 4.d
% We see that it can tolerate up to 100% of the modeled uncertainty, the
% results found here is eactly as we have seen in 4.c
deltamin = ucomplex('delta',0,'Radius',0.1454);
%we use the smallest value found in 4.c
[stabmarg,destabunc,report] = robuststab(feedback(P+deltamin,C))

%% 4.e
deltamin = -1/S(1,1)*V(:,1)*U(:,1)';
deltanew = cnum2sys(deltamin,freq);
pole(feedback(P+deltanew, C))
%We can see clearly that there are two poles on the imaginary
% axis, thus the dynamic system generated in 4.e is unstable.

%% 4.f
deltanew = ultidyn('delta', [1 1], 'Bound', norm(deltamin));
[stabmarg,destabunc,report] = robuststab(feedback(P+deltanew,C))

%% 4.g
% G is the sensitivity function of the negative feedback system of P,C
P = tf(1,[1 -1]);
C = tf([5.8 9],[0.04 1 0]);
G_new = P*C/(1+P*C);

%% 4.h
[Inf_norm, freq] = norm(G_new, inf, normTol);
M = freqresp(G_new, freq);
[U,S,V] = svd(M);
deltamin = -1/S(1,1)*V(:,1)*U(:,1)';
fprintf('norm of delta\n')
disp(norm(deltamin))
fprintf('norm of delta calculated by small gain thm\n')
disp(1/Inf_norm)
P_tilde = P*(1+deltamin);
pole(feedback(P_tilde, C))

%% 4.i
deltanew = ultidyn('delta', [1 1], 'Bound', norm(deltamin));
uncertainSys = feedback((P*(1+deltanew)),C);
[stabmarg,destabunc,report] = robuststab(uncertainSys)
% By refering to the output in stabmarg and destabunc, the
% system is marginally unstable, which meets the result found in 4.h

%% q5

%% 5.a
P = tf(1,[1 -1]); 
C = tf([5.8 9],[0.04 1 0]);
S = feedback(1,P*C);
Wp = tf([0.667 3],[1 0.003]);
normTol = 0.001;
norm(Wp*S,inf,normTol)<=1 

%% 5.b
figure
bodemag(S,'r--',1/Wp,'bo')
legend('S','1/Wp')

%% 5.c
% Based on the results of 4.g, the system is stable, moreover we can verify
% it by using the function robstab.
delta = ultidyn('delta',[1 1],'bound',1);
Pu = P*(1+0.4*delta);
System1 = feedback(C,Pu);
[stabmarg,wcu,report] = robuststab(System1)

%% 5.d
[wcg,wcu] = wcgain(Wp/(1+Pu*C))
% specific stable linear system
wcu.delta

%% 5.e
% We can see the worst case gain from the intersection of 1/Wp and the
% sensitivity of the worse case Pu (refer to the small gain theorem), the
% intersection frequency is larger than 1, which confirm the worse case
% gain seen in 5.d. 
Pu_worst = P*(1+0.4*wcu.delta);
Sworst = 1/(1+Pu_worst*C);
figure
bodemag(S,'r-',Sworst,'go',1/Wp,'y.')
legend('Sensitivity for P','Sensitivity for the worst Pu','1/Wp')


%% 5.f
figure
step(S,Sworst,4)
legend('Step response of S for P','Step response of S for Worst case Pu')

%% 5.g Task 1
Wu = makeweight(0.4, 20, 400);
delta = ultidyn('delta',[1 1],'bound',1);
PuNew = P*(1+Wu*delta);
robuststab(feedback(PuNew,C))

%% 5.g Task 2
[wcg,wcu] = wcgain(Wp/(1+PuNew*C))

%% 5.g Task 3
figure
Pu_worst = P*(1+Wu*wcu.delta);
step(1/(1+P*C),1/(1+Pu_worst*C),4)

