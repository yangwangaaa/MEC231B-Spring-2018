%% Robustness of Servo Controller for DC Motor
% This example shows how to use uncertain objects in Robust Control
% Toolbox(TM) to model uncertain systems. A follow-on to this demo
% will show how to automate robustness calculations using the robustness 
% analysis tools.

%   Copyright 1986-2012 The MathWorks, Inc.
% UC Berkeley, ME C231B/EECS C220C, Spring 2017

%% Data Structures for Uncertainty Modeling
% Robust Control Toolbox lets you create uncertain elements, such as
% physical parameters whose values are not known exactly, and combine these
% elements into uncertain models.  You can then easily analyze the impact of
% uncertainty on the control system performance.
%
% For example, consider a plant model
%
% $$P(s) = \frac{\gamma}{\tau s + 1}$$
%
% where |gamma| can range in the interval [3,5] and |tau| has average
% value 0.5 with 30% variability.  You can create an uncertain model of P(s)
% as in this example:

gamma = ureal('gamma',4,'range',[3 5]);
tau = ureal('tau',.5,'Percentage',30);
P = tf(gamma,[tau 1])

%%
% Suppose we have designed an integral controller |C| for the nominal 
% plant (|gamma|=4 and |tau|=0.5).  To find out how variations of 
% |gamma| and |tau| affect the plant and the closed-loop performance,
% we form the closed-loop system |CLP| from |C| and  |P|.

KI = 1/(2*tau.Nominal*gamma.Nominal);
C = tf(KI,[1 0]);
CLP = feedback(P*C,1)

%%
% We can now generate 20 random samples of the uncertain parameters 
% |gamma| and |tau| and plot the corresponding step responses of the
% plant and closed-loop models:

subplot(2,1,1); step(usample(P,20)), title('Plant response (20 samples)')
subplot(2,1,2); step(usample(CLP,20)), title('Closed-loop response (20 samples)')

%%
% *Figure 1:* Step responses of the plant and closed-loop models

%%
% The bottom plot shows that the closed-loop system is reasonably robust despite
% significant fluctuations in the plant DC gain.  This is a desirable, and
% common characteristic of a properly designed feedback system.

%% DC Motor Example with Parameter Uncertainty and Unmodeled Dynamics
% Now we'll build on the Control System Toolbox(TM) DC motor example by adding
% parameter uncertainty and unmodeled dynamics, and investigating the
% robustness of the servo controller to such uncertainty.
%
% The nominal model of the DC motor is defined by the resistance R, the
% inductance L, the emf constant Kb, armature constant Km, the linear
% approximation of viscous friction Kf and the inertial load J.  Each of
% these components varies within a specific range of values.  The resistance
% and inductance constants range within +/- 40% of their nominal values.
% We use the |ureal| function to construct these uncertain parameters:

R = ureal('R',2,'Percentage',40);
L = ureal('L',0.5,'Percentage',40);

%%
% For physical reasons, the values of Kf and Kb are the same, even if they
% are uncertain.  In this example, the nominal value is 0.015 with a range
% between 0.012 and 0.019.

K = ureal('K',0.015,'Range',[0.012 0.019]);
Km = K;
Kb = K;

%%
% Viscous friction, Kf, has a nominal value of 0.2 with a 50% variation in
% its value. 

Kf = ureal('Kf',0.2,'Percentage',50);

%% Electrical and Mechanical Equations
% The current in the electrical circuit, and the torque applied to the
% rotor can be expressed in terms of the applied voltage and the angular
% speed. Create the transfer function |H| relating these variables, and
% make |AngularSpeed| an output of |H| for later use:

H = [1;0;Km] * tf(1,[L R]) * [1 -Kb] + [0 0;0 1;0 -Kf];
H.InputName = {'AppliedVoltage';'AngularSpeed'};
H.OutputName = {'Current';'AngularSpeed';'RotorTorque'};

%%
% H is an multi-input, multi-output uncertain system as seen from its
% display.
H

%%
% The motor typically drives an inertia, whose dynamic characteristics
% relate the applied torque to the rate-of-change of the angular speed.
% For a rigid body, this is a constant.   A more realistic, but uncertain, model
% might contain unknown damped resonances.  Use the |ultidyn| object to model
% uncertain linear time-invariant dynamics.   The nominal value
% of the rigid body inertia is set to 0.02 and we add 15% dynamic uncertainty
% in multiplicative form:

J = 0.02*(1 + ultidyn('Jlti',[1 1],'Type','GainBounded','Bound',0.15,...
   'SampleStateDim',4));

%% Uncertain Model of DC Motor
% It is a simple matter to relate the |AngularSpeed| input to the |RotorTorque|
% output through the uncertain inertia, |J|, using the |lft| command.
% The |lft| command below connects the last output of |H|, RotorTorque,
% to the input of the uncertain inertia. It connects the output of the
% uncertain inertia to the last input of |H|, AngularSpeed.
%
% The AngularSpeed input equals RotorTorque/(J*s), hence "positive"
% feedback from the 3rd output to the 2nd input of |H| is used
% to make the connection.  This results in a system with 1 input
% (|AppliedVoltage|) and 2 outputs, (|Current| and |AngularSpeed|).

Pall = lft(H, tf(1,[1 0])/J);

%%
% Select only the AngularSpeed output for the remainder of the control
% analysis.

P = Pall(2,:)

%%
% P is a single-input, single-output uncertain model of the DC motor.


%% Open-Loop Analysis
% First, let's compare the step response of the nominal DC motor with 20 samples of the 
% uncertain model of the DC motor:

clf
step(P.NominalValue,'r-+',usample(P,20),'b',3)
legend('Nominal','Samples')
 
%%
% *Figure 2:* Open-loop step response analysis

%%
% Similarly, we can compare the Bode plot of the open-loop nominal (red) and 
% sampled (blue) uncertain models of the DC motor.

om = logspace(-1,2,80);
Pg = ufrd(P,om);
bode(usample(Pg,25),'b',Pg.NominalValue,'r-+');
legend('Samples','Nominal')

%%
% *Figure 3:* Open-loop Bode plot analysis

%% Attribution
% Copyright 2016-17, Andy Packard.  This work is licensed under the Creative
% Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. To
% view a copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to
% Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.

%% File Information
disp(mfilename)