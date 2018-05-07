%% Flexible structure: Nominal Design
% Problem features
%
% * 4 mode flexible structure
% * 2 (uncollocated) sensor/actuators.
%
% Formulate disturbance rejection problem, and solve for optimal H_infinity
% controller.  Examine effect of suboptimality on controller bandwidth.

% UC Berkeley, ME C231B/EECS C220C, Spring 2017

%% Settings
format compact

%% Lumped parameter values
M = 1;
K = 1;
xi = 0.005;
C = 2*sqrt(M*K)*xi;

%% Distribute these over N elements
N = 4;
Kv = repmat(N*K,[1 N]);
Mv = repmat(M/N,[1 N]);
Cv = repmat(C,[1 N]);
Cv = repmat(N*C,[1 N]);

%% Mass, Damping and Stiffness Matrices
% Standard $M \ddot z + C \dot z + K z = F$ model.  Define $x = [z ; \dot
% z]$.  Form $M$, $C$ and $K$ matrices.
Mmat = diag(Mv)
%%
Cmat = -[-Cv(1)-Cv(2) Cv(2) 0 0;...
   Cv(2) -Cv(2)-Cv(3) Cv(3) 0;...
   0 Cv(3) -Cv(3)-Cv(4) Cv(4);...
   0  0 Cv(4) -Cv(4)];
Kmat = -[-Kv(1)-Kv(2) Kv(2) 0 0;...
   Kv(2) -Kv(2)-Kv(3) Kv(3) 0;...
   0 Kv(3) -Kv(3)-Kv(4) Kv(4);...
   0  0  Kv(4) -Kv(4)];

%% Create A and B matrices of uncertain system.
% Model has independent forces at each mass as inputs.  All states are
% listed as outputs.  $A$ matrix depends on mass, damping and stiffness
% matrices, while $B$ only depends on mass matrix.  Since outputs consist
% of all states, $C$ and $D$ matrices are simple.
A = [zeros(N,N) eye(N); inv(Mmat)*[-Kmat -Cmat]];
B = [zeros(N,N); inv(Mmat)*eye(4)];
C = eye(2*N);
S = ss(A,B,C,0);

%% Check system poles
pole(S)

%% Bode magnitude from $F_2$ to $z_4$ (open-Loop)
bodemag(S(4,2))
ylim([-50 50])
xlim([0.1 10])

%% Model disturbance and actuator forces
% Disturbance acts directly on mass #2.  Control inputs act on mass #1, and
% between mass 2 and 3.
Bd = [0;1;0;0];
Bu = [[1;0;0;0] [0;-1;1;0]];

%% Model sensed and regulated vaiables
% Sensors measure z4, and relative deflection at 1/2, so z2-z1.   Regulated
% variable is z4.
Ce = [0 0 0 1 0 0 0 0];
Cy = [0 0 0 1 0 0 0 0;-1 1 0 0 0 0 0 0];

%% Form system
G = [Ce;Cy]*S*[Bd Bu];

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

%% Hinf synthesis (iterate to optimal)
nY = 2;
nU = 2;
[K,CLP,gAch] = hinfsyn(P,nY,nU);
%%
% Note that "optimal" controller is by itself, unstable
isstable(K)
%%
% Of course, the closed-loop is stable...
isstable(lft(P,K))
%%
% The closed-loop Hinf-norm agrees with |gAch|
[norm(lft(P,K),inf) gAch]

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

%% Plot unweighted open and closed-loop responses, from [f2;n1;n2] to the
% errors, which are [z4;u1;u2].  In open-loop [n1;n2] have no effect, and
% [u1;u2] are identically zero, so the only open-loop response is from f2
% to z4.  Closed-loop involves all signals.
h = bodeplot(Puw(1:3,1:3),'r',lft(Puw,K),'g');
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

%% BodeMagnitude (unweighted closed-loop), d2-to-z4
clf
uwCLP = lft(Puw,K);
bodemag(uwCLP(1,1),'k')
legend('d2-to-z4')
ylim([-50 10]);
xlim([0.01 10]);

%% "backoff" from optimal by 2%
boFac = 1.02;
[Kb,~,gAchb] = hinfsyn(P,nY,nU,'gmax',boFac*gAch,'gmin',boFac*gAch);
fprintf('Optimal Norm, %6.4g;    Achieved with 2%% backoff %6.4g\n',gAch,gAchb)
bodemag(K,'r',Kb,'b');
legend('Optimal','2% backoff', 'location', 'best')

%% Closed-loop Time-domain simulation: both controllers
clf
TF = 8;
step(lft(Puw(:,[1 4 5]),K), 'k', lft(Puw(:,[1 4 5]),Kb), 'b',TF)
legend('Optimal','2% backoff')

%% K: Check classical and disc loop-at-a-time margins at plant input
% |loopmargin(P,K)| uses negative feedback.  In this example, the
% generalized open-loop "plant" (called |P|, formed with |sysic|, and the
% argument to the |hinfsyn| call) does not include a negative-sign in the
% measured variables.  Hence, the call to |loopmargin| must include a
% negative-sign, so that when |loopmargin| creates the negative-feedback
% loop, it actually creates the closed-loop system that |hinfsyn| has
% designed.
Gyu = G(2:3,2:3);
[ci,di] = loopmargin(Gyu,-K,'ci,di');
%%
% The gain and phase margins in channel 1 are very small.  This
% design is probably not suitable for real-world implementation
ci(1)
%%
di(1)
%%
% Similarly, the gain and phase margins in 2nd channel are very small.  This
% design is probably not suitable for real-world implementation
ci(2)
%%
di(2)
%%
% You can verify with an additional calculation that the margins acheived
% by |Kb| are similar - fine in channel 1, but poor in channel 2.

%% Challenge
% Redo the design with different measurements.  For example, try measuring
% z4 and z3-z2 (relative displacement at 2/3 junction).  Verify if the same
% general conclusions hold as in the case we completed.

%% File Information
disp(mfilename)

