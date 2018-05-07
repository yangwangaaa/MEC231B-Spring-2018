%% FlightControlExWithSSV.m
% This example computes the robust stability margins and worst-case 
% disturbance rejection performance for a rigid body transport aircraft 
% with an output feedback control law.  The uncertainty model has 14 real 
% parametric uncertainties associated with the aerodynamic coefficients 
% of the aircraft. The uncertainty model also includes unmodeled actuator
% dynamics. This example is taken from the following textbook:
%
% * _A Practical Approach to Robustness Analysis with Aeronautical
%    Applications_ by G. Ferreres, Kluwer, 1999.
%
% UC Berkeley, ME C231B/EECS C220C, Spring 2017

%% Create Open-Loop Model
% The open-loop model includes the rigid-body aircraft dynamics. 
% These dynamics contain uncertainties in the aerodynamic coefficients. 
% The open-loop model also includes dynamics for the rudder and aileron
% actuators.  The actuator models also include dynamic uncertainty.

%%
% The model parameters are (Appendix A.1, p.187):
deg2rad = pi/180;
rad2deg = 1/deg2rad;
gV = 0.146418;         % g/V
tan_theta0 = 0.14;     % tan(theta0)
alpha0 = 8*deg2rad;    % (rad)

%%
% The uncertain aerodynamic coefficients are:
Ybeta = ureal('Ybeta',-0.082,'Percentage',10);
Yp = ureal('Yp',0.010827,'Percentage',10);
Yr = ureal('Yr',0.060268,'Percentage',10);
Ydeltap = ureal('Ydeltap',0.002,'Percentage',10);
Ydeltar = ureal('Ydeltar',0.0118,'Percentage',10);
Lbeta = ureal('Lbeta',-0.84,'Percentage',10);
Lp = ureal('Lp',-0.76,'Percentage',10);
Lr = ureal('Lr',0.74,'Percentage',10);
Ldeltap = ureal('Ldeltap',0.095,'Percentage',10);
Ldeltar = ureal('Ldeltar',0.06,'Percentage',10);
Nbeta = ureal('Nbeta',0.092,'Percentage',10);
Np = ureal('Np',-0.23,'Percentage',10);
Nr = ureal('Nr',-0.114,'Percentage',10);
Ndeltar = ureal('Ndeltar',-0.151,'Percentage',10);

%%
% The states, inputs, and outputs are given by:
%
% * States = [beta; p; r; phi] = [sideslip; roll rate; yaw rate; roll angle]
% * Inputs = [deltap; deltar] = [ aileron deflection; rudder deflection]
% * Outputs = [ny; p; r; phi] = [accel; roll rate; yaw rate; roll angle]
%
% The state equations (See Eq 2.1/2.2 on p.30 and p.188) are:

A = [Ybeta (Yp+sin(alpha0)) (Yr-cos(alpha0)) gV; ...
    Lbeta  Lp Lr 0; Nbeta Np Nr 0; 0 1 tan_theta0 0];
B = [Ydeltap Ydeltar; Ldeltap Ldeltar; 0 Ndeltar; 0 0];
C = [-1/gV*deg2rad*[Ybeta Yp Yr 0]; zeros(3,1) eye(3)];
D = [-1/gV*deg2rad*[Ydeltap Ydeltar]; zeros(3,2)];
AIRCRAFT = ss(A,B,C,D);

%%
% The nominal models for rudder and aileron actuators are:
N1 = [-1.77, 399];
D1 = [1 48.2 399];
deltap_act_nom = tf(N1,D1);

N2 = [2.6 -1185 27350];
D2 = [1 77.7 3331 27350];
deltar_act_nom = tf(N2,D2);

%%
% Include multiplicative uncertainty on the actuator dynamics.
% The uncertainty weight for each actuator specifies 20% uncertainty at 
% low frequencies and 200% uncertainty at high frequencies.  The 
% uncertainty weights for each actuator cross 100% uncertainty at roughly
% |2*wB| where |wB| denotes the -6dB bandwidth of the actuator.
Wup = makeweight(0.20,32,2);
deltap_act = deltap_act_nom*(1+Wup*ultidyn('Delp_act',[1 1]));

%%
% The rudder actuator has similar uncertainty.
Wur = makeweight(0.20,70,2);
deltar_act = deltar_act_nom*(1+Wur*ultidyn('Delr_act',[1 1]));

%% Compose open-loop model
% The uncertain aircraft model is
P = AIRCRAFT*blkdiag(deltap_act,deltar_act);

%% Control Law
% A constant gain output feedback law is used.  The gain below is taken
% from documentation in the SMT toolbox.
K = [-629.8858 11.5254 3.3110 9.4278; ...
  285.9496 0.3693 -2.6301 -0.5489];

%% Closed-Loop
% Use the FEEDBACK command to form the closed-loop system from an input
% disturbance to plant output.  The closed-loop is nominally stable. 
CLOOP = feedback(P,K);
isstable(CLOOP.Nominal)

%% Robust Stabilty
% The closed-loop is uncertain due to the real parametric and dynamic 
% uncertainty in the plant model. An uncertain frequency response for 
% the closed-loop system is computed with the UFRD command. Robust
% stability of the closed loop is assessed using the ROBUSTSTAB command.
% The closed loop is robustly stable and can tolerate up to 276% of the 
% modeled uncertainty.
w = logspace(-1,3,100);
CLOOPfr = ufrd(CLOOP,w);
[stabmarg,destabunc,report,info] = robuststab(CLOOP);

%% Robust stability calculation is a MUSSV analysis
% The ROBUSTSTAB performs a mu analysis in order to compute the 
% stability margins.  The plot below shows the mu upper and lower
% bounds which are stored in the info structure.   The stability margins
% returned by ROBUSTSTAB are inversely related to the peak of the mu
% lower/upper bound plots.
figure
semilogx(info.MussvBnds(1),'b',info.MussvBnds(2),'r')
xlim([w(1) w(end)])
xlabel('Frequency (rad/sec)');
ylabel('Mu Bounds');
%%
stabmarglb = 1/norm(info.MussvBnds(1),inf)
%%
stabmargub = 1/norm(info.MussvBnds(2),inf)
%%
stabmarg % from robuststab

%% Conclusions
% In this example, we relate the bounds from |robuststab| with
% frequency-dependent structured-singular plots, which are also returned by
% |robuststab| in the |info| structure.

%% Attribution
% Copyright 2016-17, Andy Packard.  This work is licensed under the Creative
% Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. To
% view a copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to
% Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.

%% File Information
disp(mfilename)