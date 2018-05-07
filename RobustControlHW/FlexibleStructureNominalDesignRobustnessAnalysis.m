%% Flexible structure: Nominal Design
% Problem features
%
% * 4 mode flexible structure
% * 2 (uncollocated) sensor/actuators.
% * uncertain mass, damping and stiffness parameters
% * unmodeled dynamics in actuator models.
%
% Formulate disturbance rejection problem, and solve for nominal (no
% uncertainty) plant.  Assess robust stability of closed-loop system using 
% the actual uncertainty model.  

%% Parameter uncertainty in Flexible structure
% Model properties
%
% * 4 mode flexible structure
% * uncertain mass, damping and stiffness parameters
%
% Use |ureal| parameters (12 in all) to obtain |uss| model.  All parameters
% enter in _rank=1_ fashion, so each has only one copy in final LFT model

% UC Berkeley, ME C231B/EECS C220C, Spring 2018

%% Lumped parameter values
M = 1;
K = 1;
xi = 0.01;
C = 2*sqrt(M*K)*xi;

%% Distribute these over N elements
N = 4;
Kv = repmat(N*K,[1 N]);
Mv = repmat(M/N,[1 N]);
Cv = repmat(C,[1 N]);

%% Add uncertainty to each term
% For east access, store each uncertain element in an array, with
% appropriate name, for use later in building the state-space model.
percK = 5;
Kv = [ureal('k1',Kv(1),'Percentage',percK), ...
   ureal('k2',Kv(2),'Percentage',percK), ...
   ureal('k3',Kv(3),'Percentage',percK), ...
   ureal('k4',Kv(4),'Percentage',percK)];
perkM = 5;
Mv = [ureal('m1',Mv(1),'Percentage',perkM), ...
   ureal('m2',Mv(2),'Percentage',perkM), ...
   ureal('m3',Mv(3),'Percentage',perkM), ...
   ureal('m4',Mv(4),'Percentage',perkM)];
perkC = 5;
Cv = [ureal('c1',Cv(1),'Percentage',perkC), ...
   ureal('c2',Cv(2),'Percentage',perkC), ...
   ureal('c3',Cv(3),'Percentage',perkC), ...
   ureal('c4',Cv(4),'Percentage',perkC)];

%% Mass, Damping and Stiffness Matrices
% Standard $M \ddot z + C \dot z + K z = F$ model.  Define $x = [z ; \dot
% z]$.  Form $M$, $C$ and $K$ matrices, all of which are uncertain.
Mmat = diag(Mv)
%%
Cmat = -[-Cv(1)-Cv(2) Cv(2) 0 0;...
   Cv(2) -Cv(2)-Cv(3) Cv(3) 0;...
   0 Cv(3) -Cv(3)-Cv(4) Cv(4);...
   0  0 Cv(4) -Cv(4)]
%%
Kmat = -[-Kv(1)-Kv(2) Kv(2) 0 0;...
   Kv(2) -Kv(2)-Kv(3) Kv(3) 0;...
   0 Kv(3) -Kv(3)-Kv(4) Kv(4);...
   0  0  Kv(4) -Kv(4)]

%% Create A and B matrices of uncertain system.
% Model has independent forces at each mass as inputs.  All states are
% listed as outputs.  $A$ matrix depends on mass, damping and stiffness
% matrices, while $B$ only depends on mass matrix.  Since outputs consist
% of all states, $C$ and $D$ matrices are simple.
A = [zeros(N,N) eye(N); inv(Mmat)*[-Kmat -Cmat]];
B = [zeros(N,N); inv(Mmat)*eye(4)];
C = eye(2*N);
S = ss(A,B,C,0)

%% Check nominal system poles
pole(S.NominalValue)

%% Bode magnitude from $F_2$ to $z_4$
% Plot 25 samples of Bode Magnitude plot
nSamples = 25;
bodemag(usample(S(4,2),nSamples))
ylim([-50 50])
xlim([0.1 10])

%% Model disturbanace and actuator forces
% Disturbance acts directly on mass #2.  Control inputs act on mass #1, and
% between mass 2 and 3.
Bd = [0;1;0;0];
Bu = [[1;0;0;0] [0;-1;1;0]];

%% Model sensed and regulated vaiables
% Sensors meaure z4, and relative deflection at 1/2, so z2-z1.   Regulated
% variable is z4.
Ce = [0 0 0 1 0 0 0 0];
Cy = [0 0 0 1 0 0 0 0;-1 1 0 0 0 0 0 0];

%% Form system
G = [Ce;Cy]*S*[Bd Bu];

%% Add uncertainty (dynamic) in actuators
wDelta = makeweight(0.15, 100, 2);
G = G*blkdiag(1, eye(2)+wDelta*ultidyn('umD',[2 2]));

%% Performance weights
% Main performance objective is regulation of z4.  Weight wP is applied to
% z4, and has large value at low-frequency, aiming for nearly prefect
% steady-state regulation in the presence of constant disturbances.
% Identical noise weights are used to model spectrum of sensor noise at
% both measurements.
wP = makeweight(1000,7,0.75);
wN = makeweight(0.001,20,1.25)*eye(2);

%% Weighted open-loop interconnection
PlantModel = G;
systemnames = 'PlantModel wP wN';
inputvar = '[f2;n1;n2;u{2}]';
outputvar = '[wP;u;PlantModel(2)+wN(1);PlantModel(3)+wN(2)]';
input_to_wP = '[PlantModel(1)]';
input_to_wN = '[n1;n2]';
input_to_PlantModel = '[f2;u]';
P = sysic;
%%
% Capture nominal model for synthesis
Pnom = P.NominalValue;

%% Hinf synthesis (iterate to optimal)
nY = 2;
nU = 2;
[K,CLP,gAch] = hinfsyn(Pnom,nY,nU);
%%
% Note that "optimal" controller is by itself, unstable
isstable(K)
%%
% Of course, the closed-loop is stable...
isstable(lft(Pnom,K))

%% Unweighted open-loop interconnection: nominal model
% Set weights to 1, and use same interconnection code
wP = ss(1);
wN = ss(eye(2));
PlantModel = G;
systemnames = 'PlantModel wP wN';
inputvar = '[f2;n1;n2;u{2}]';
outputvar = '[wP;u;PlantModel(2)+wN(1);PlantModel(3)+wN(2)]';
input_to_wP = '[PlantModel(1)]';
input_to_wN = '[n1;n2]';
input_to_PlantModel = '[f2;u]';
Puw = sysic;
%%
% Capture nominal model as well
PuwNom = Puw.NominalValue;

%% Plot unweighted open and closed-loop responses, from [f2;n1;n2] to the
% errors, which are [z4;u1;u2].  In open-loop [n1;n2] have no effect, and
% [u1;u2] are identically zero, so the only open-loop response is from f2
% to z4.  Closed-loop involves all signals.
h = bodeplot(PuwNom(1:3,1:3),'r',lft(PuwNom,K),'g');
setoptions(h,'PhaseVisible','off');
allYlim = getoptions(h,'yLim');
for i=1:numel(allYlim)
    allYlim{i} = [-60 20];
end;
allXlim = getoptions(h,'xLim');
for i=1:numel(allXlim)
    allXlim{i} = [0.01 10];
end;
setoptions(h,'yLim',allYlim,'xLim',allXlim);

%% Samples of BodeMagnitude (unweighted closed-loop), d2-to-z4
clf
uwCLP = lft(Puw,K);
bodemag(uwCLP(1,1),'y',uwCLP(1,1).NominalValue,'k')
legend('Samples','Nominal')
ylim([-50 50]);
xlim([0.01 10]);
%%
% Clearly the robustness of the performance is quite poor.

%% "backoff" from optimal by 10%
Kb = hinfsyn(P,nY,nU,'gmax',1.1*gAch,'gmin',1.1*gAch);

%% Stability margin of "optimal" controller
[SM,DU] = robuststab(lft(G,K));
SM
%%
% Verify that |DU| destabilizes at reported frequency
pole(lft(usubs(G,DU),K))

%% Stability margin of "backed-off" controller
% Note that the stability margin improves, although it is still not
% robustly stable to modeled uncertainty.
SMb = robuststab(lft(G,Kb))

%% Closed-loop Time-domain simulation: both controllers
clf
TF = 8;
step(lft(Puw(:,[1 4 5]),K),'r',lft(Puw(:,[1 4 5]),Kb),'k--',TF)
legend('Original','10% backoff')
%%
% Again, now from time-domain, additional evidence that the robustness of
% the performance is quite poor.

%% File Information
disp(mfilename)
