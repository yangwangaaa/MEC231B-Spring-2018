%% Robust Stability, Mu Analysis, and Worst-Case Gain
% This demo shows how to use Robust Control Toolbox(TM) to analyze and
% quantify the robustness of feedback control systems, including modeling
% errors and parameter variations.  We'll look at how to test for robust
% stability with the |robuststab| functions and gain insight into the
% connection with mu analysis and the |mussv| function.   In addition the
% worst-case gain of the system is evaluated with the |wcgain| function.
%
% UC Berkeley, ME C231B/EECS C220C, Spring 2017

%   Copyright 1986-2011 The MathWorks, Inc.
%   $Revision: 1.1.6.7.2.1 $ $Date: 2011/07/01 16:40:30 $

%% System Description
% Figure 1 shows the block diagram of a closed-loop system. The plant model
% P is uncertain, the plant output yP must be regulated to remain small in
% the presence of disturbances d and measurement noise n.  Weighting functions
% Wd and Wn are used to define the performance objective as a closed-loop
% norm,
%
% $$\left\| P(1+KP)^{-1}W_d,  \ \ \ -PK(1+PK)^{-1}W_n \right\|_{\infty}$$
%
% In this example, Wd is large at low frequencies, and Wn is large at high
% frequencies.
%
% <<..\Images\rsrpmuIC.png>>
% 
% *Figure 1*: Block Diagram 

%% Creating an Uncertain Plant Model
% An uncertain plant model is a lightly-damped, second-order system with
% parametric uncertainty in the denominator coefficients and significant
% frequency-dependent unmodeled dynamics beyond 6 radians/second. The
% mathematical model looks like this: 
%
% $$\frac{16}{s^2+0.16 s+k} (1+W_u(s) \delta(s))$$
%
% The parameter |k| is assumed to be about 30% uncertain, with a nominal
% value of 16. The frequency-dependent uncertainty at the plant input is
% assumed to be about 30% at low frequency, rising to 100% at 10
% radians/second, and larger uncertainty beyond that. Construct the 
% uncertain plant model |P| by creating and combining the uncertain
% elements:

k = ureal('k',16,'Percentage',30);
delta = ultidyn('delta',[1 1],'SampleStateDim',4);
Wu = makeweight(0.3,10,20);
P = tf(16,[1 0.16 k])*(1+Wu*delta);

%% Designing a Controller
% We use the controller designed in the demo
% "Improving Stability While Preserving Open-Loop Characteristics".  The
% plant model used there happens to be the nominal value of the uncertain
% plant model created above.  For completeness, we repeat the commands used to
% generate the controller. The classical margins from |allmargin| indicate 
% good stability robustness to unstructured gain/phase variations within 
% the loop.  
K_PI = pid(1,0.8);
K_rolloff = tf(1,[1/20 1]);
Kprop = K_PI*K_rolloff;
[Gamma,K] = cprop2cact(P.NominalValue,-Kprop);

allmargin(P.NominalValue*K)

%% Stability Robustness with Respect to Modeled Uncertainty
% Assessing stability robustness to the explicitly modeled uncertainty in 
% |P| requires an alternate analysis using the command |robuststab|. Since
% this calculation does not involve other performance objectives (such as
% disturbance rejection and/or noise insensitivity), any representation of
% the closed-loop system is adequate for analysis.  Here we use |feedback|
% to form the closed-loop system from input disturbance to plant output,
% and then assess the stability robustness
% to the modeled uncertainty using |robuststab|.
D2Yclp = feedback(P,K);
stabmarg = robuststab(D2Yclp)
%%
% The margin is greater than 1, indicating that the closed-loop system is
% stable for all of the modeled uncertainty.


%% Expressing Performance as a Weighted Norm
% The robust performance analysis tools focus on how norms of closed-loop
% transfer functions are affected by model uncertainty.  Hence, we
% introduce frequency-dependent weighting functions to combine the multiple
% competing objectives (disturbance rejection and noise insensitivity) into
% one overall objective, 
%
% $$\left\| P(1+KP)^{-1}W_d,  \ \ \ -PK(1+PK)^{-1}W_n \right\|_{\infty}$$
%
clf
Wd = makeweight(100,.4,.15);
Wn = makeweight(0.5,20,100);
bodemag(Wd,'b--',Wn,'k--')
title('Performance Weighting Functions')
legend('Input Disturbance','Measurement Noise')

%% Creating a Closed-Loop System with CONNECT
% Use |connect| to build an uncertain model of the closed-loop system.
% Name the signals coming in and out of each block and let |connect|
% do the wiring:

C = K;
P.u = 'uP';  P.y = 'yP';
C.u = 'uC';  C.y = 'yC';
S1 = sumblk('uP = yC + D');
S2 = sumblk('uC = -yP - N');
Wn.u = 'n'; Wn.y = 'N';
Wd.u = 'd'; Wd.y = 'D';
ClosedLoop = connect(P,C,S1,S2,Wn,Wd,{'d','n'},'yP');

%%
% The variable |ClosedLoop| is an uncertain system with two inputs and
% one output. It depends on two uncertain elements: a real parameter |k|
% and an uncertain linear, time-invariant dynamic element |delta|.

ClosedLoop


%% Robust Stability Analysis
% Does the closed-loop system remain stable for all values of |k|, |delta| 
% in the ranges specified above?  Yes, as this assessment was performed earlier,
% using the unweighted closed-loop system.  Here, we repeat the analysis
% using |ClosedLoop|, which includes performance weights, to reiterate that
% the stability assessment is unaffected by performance objectives.  As
% expected, |stabmarg| is identical to it's earlier calculated value.
[stabmarg,destabunc,report,info] = robuststab(ClosedLoop);

stabmarg

%%
% The variable |stabmarg| gives upper and lower bounds on the *robust
% stability margin*, a measure of how much uncertainty on |k|, |delta| the
% feedback loop can  tolerate before becoming unstable. The margin is 1.5,
% which means that the closed loop will remain stable for up to 150% of the
% specified uncertainty. The variable |report| summarizes these results:   

report

%%
% The variable |destabunc| contains the combination of (|k|,|delta|),
% closest to their nominal values, that causes instability.  

destabunc

%%
% We can substitute these values into |ClosedLoop|, and verify that these
% values cause the closed-loop system to be unstable:
damp(usubs(ClosedLoop,destabunc))

%%
% Note that the natural frequency of the unstable closed-loop pole is
% given by |stabmarg.DestabilizingFrequency|:

stabmarg.DestabilizingFrequency

%% Connection with Mu Analysis
% The structured singular value, or mu, is the mathematical tool used 
% by |robuststab| to compute the robust stability margin. If you are 
% comfortable with structured singular value analysis, you can use 
% the |mussv| function directly to compute mu as a function of frequency
% and reproduce the results above. The function |mussv| is the underlying
% engine for all robustness analysis commands. 
%
% To use |mussv|, we first extract the |(M,Delta)| decomposition of the 
% uncertain closed-loop model |ClosedLoop|, where |Delta| is a 
% block-diagonal matrix of (normalized) uncertain elements.
% The 3rd output argument of |lftdata|, |BlkStruct|, describes the block-diagonal
% structure of |Delta| and can be used directly by |mussv|

[M,Delta,BlkStruct] = lftdata(ClosedLoop);

%%
% For a robust stability analysis, only the channels of |M| associated
% with the uncertainty channels are used. Based on the row/column size of
% |Delta|, select the proper columns and rows of |M|. Remember that the
% rows of |Delta| correspond to the columns of |M|, and vice versa.
% Consequently, the column dimension of |Delta| is used to specify the rows
% of |M|:    

szDelta = size(Delta);
M11 = M(1:szDelta(2),1:szDelta(1));

%%
% Mu-analysis is performed on a finite grid of frequencies. For comparison
% purposes, evaluate the frequency response of |M11| over the same
% frequency grid as used for the |robuststab| analysis.

omega = info.Frequency;
M11_g = frd(M11,omega);

%%
% Compute |mu(M11)| at these frequencies and plot the resulting lower and
% upper bounds: 

mubnds = mussv(M11_g,BlkStruct,'s');

LinMagopt = bodeoptions;
LinMagopt.PhaseVisible = 'off'; LinMagopt.XLim = [1e-1 1e2]; LinMagopt.MagUnits = 'abs';
bodeplot(mubnds(1,1),mubnds(1,2),LinMagopt);
xlabel('Frequency (rad/sec)');
ylabel('Mu upper/lower bounds');
title('Mu plot of robust stability margins (inverted scale)');

%%
% *Figure 3:* Mu plot of robust stability margins (inverted scale)

%%
% The robust stability margin is the reciprocal of the structured singular
% value. Therefore upper bounds from |mussv| become lower bounds
% on the stability margin. Make these conversions and find the
% destabilizing frequency where the mu upper bound peaks (that is, where
% the stability margin is smallest):  

[pkl,wPeakLow] = norm(mubnds(1,2),inf);
[pku] = norm(mubnds(1,1),inf);
SMfromMU.LowerBound = 1/pku;
SMfromMU.UpperBound = 1/pkl;
SMfromMU.DestabilizingFrequency = wPeakLow;

%%
% Compare |SMfromMU| to the robust stability margin bounds |stabmarg| computed with 
% |robuststab|; they are identical:

stabmarg
SMfromMU

%%
% Finally, note that the same mu bounds |mubnd| are available in the |info|
% output of |robuststab|: 

bodeplot(info.MussvBnds(1,1),info.MussvBnds(1,2),LinMagopt)
xlabel('Frequency (rad/sec)');
ylabel('Mu upper/lower bounds');
title('Mu plot of robust stability margins (inverted scale)');

%%
% *Figure 4:* Mu plot of robust stability margins (inverted scale)

%% Worst-Case Gain Analysis
% The closed-loop transfer function maps |[d;n]| to |e|. The worst-case
% gain of this transfer function indicates where disturbance rejection is
% worst. You can use |wcgain| to assess the worst (largest) value of this
% gain:

[maxgain,maxgainunc,info] = wcgain(ClosedLoop);
maxgain

%%
% The variable |maxgainunc| contains the values of the uncertain elements
% associated with maxgain.LowerBound:

maxgainunc

%%
% We can verify that this perturbation causes the closed-loop
% system to have gain at least equal to the value of |maxgain.LowerBound|:

maxgain.LowerBound
norm(usubs(ClosedLoop,maxgainunc),inf)

%%
% Note that there is a difference in the answers, and the gain is actually
% slightly larger than the lower bound. This is because
% |wcgain| uses a finite frequency grid to compute the worst-case gain,
% whereas |norm| uses estimates the peak gain more accurately. 
%
% Finally, compare the nominal and worst-case closed-loop gains:

bodemag(ClosedLoop.Nominal,'b-',usubs(ClosedLoop,maxgainunc),'r--',{.01,100})
legend('Nominal','Worst case')

%%
% *Figure 6:* Bode diagram comparing the nominal and worst-case closed-loop
% gains. 

%% 
% This analysis reveals that the nominal disturbance rejection and noise
% insensitivity properties of the controller K are not robust to the specified
% level of uncertainty.   The degree to which performance is not robust to
% the uncertainty level is quantified in the robust performance and
% worst-case gain calculations.

%% Attribution
% Copyright 2016-17, Andy Packard.  This work is licensed under the Creative
% Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. To
% view a copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to
% Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.

%% File Information
disp(mfilename)