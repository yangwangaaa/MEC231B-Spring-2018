%% (Solution) Constructing the smallest destabilizing perturbation (at Plant input)
% A nominal plant and controller are given.  There are several tasks to
% complete, culminating in construction of the smallest
% input-multiplicative plant perturbation that leads to instability.
%
% UC Berkeley, ME C231B/EECS C220C, Spring 2017

%% Plant, controller specification
% 2-input, 2-output plant and controller are specified in state-space form.
Ap = [ -0.2  10; -10  -0.2];
Bp = eye(2);
Cp = [1 8;-10 1];
P = ss(Ap,Bp,Cp,0);

Ac = blkdiag(zeros(2,2),-62.83,[-15.68 11;-11 -15.68],-14.33);
Bc = [0.5343   0.2178;...
     -0.2174   0.5255;...
      3.019   -3.903;...
      0.1649    2.265;...
     -1.308    -2.36;...
     -2.334  0.06362];
Cc = [ -4  0  -1.369  -0.02616  1.771  -3.681;...
        0 -4  -1.467  -3.616   -1.955  -0.6843];
C = ss(Ac,sqrt(3)*Bc,sqrt(3)*Cc,0);

%% Verify closed-loop system is stable
% Form all closed-loop maps with |LOOPSENS|
help loopsens
%%
H = loopsens(P,C)

%% Compute Hinf norm of Ti
% Ti is complementary sensitivity at plant input
[TiNorm,freq] = norm(H.Ti,inf,0.001)

%% Compute Frequency Response Matrix of Ti at frequency where peak occurs
M = freqresp(H.Ti,freq)

%%
% Find the smallest plant input-multiplicative uncertainty (a constant,
% complex matrix of dimension 2-by-2) which causes instability
[U,S,V] = svd(M);
Delta = -1/S(1,1)*V(:,1)*U(:,1)'; 
[norm(Delta)  1/TiNorm]

%% Perturbed closed-loop system
% Form perturbed plant, using input-multiplicative form and the
% perturbation matrix.
Ptilde = P*(eye(size(P,2))+Delta);

%% Verify Instability of perturbed closed-loop system
% Compute perturbed closed-loop system using |loopsens|, and check poles
pole(feedback(Ptilde,C))

%% Connection between unstable pole, and frequency where peak of Ti occurs
% Show connection between the two frequencies
freq

%% (Extra Credit!) Use |cnum2sys| to create real-dynamic uncertainty
% Rather than accepting a complex matrix as a destabilizing element, use
% the command |cnum2sys| to convert the complex matrix uncertainty into a
% real-valued, linear, dynamic system that leads to instability.  The Hinf
% norm of the dynamic system should match the maximum singular value of the
% original destabilizing complex matrix.
v1hat = cnum2sys(V(1,1),freq);
v2hat = cnum2sys(V(2,1),freq);
ubar1hat = cnum2sys(conj(U(1,1)),freq);
ubar2hat = cnum2sys(conj(U(2,1)),freq);
DeltaHat = -1/S(1,1)*[v1hat;v2hat]*[ubar1hat ubar2hat];
pole(DeltaHat)
[norm(DeltaHat,inf) norm(Delta)]
pole(feedback(P*(eye(2)+DeltaHat),C))

%% Frequency-dependent uncertainty
% Adopt a frequency-dependent uncertainty model, using a scalar weighting
% function, |wu|.  Take uncertainty model to be |P*(eye(2) + wu(s) Delta)|,
% using a first order |wu| with:
%
% * DC gain of 0.35 (35% uncertainty at low frequency)
% * gain of 1 at 40 rads/TimeUnit (100% uncertainty at 40)
% * High-Frequency gain of 10 (1000% uncertainty at high frequency)
% 
wu = makeweight(0.35,40,10);

%%
% Mimic/repeat procedure to find smallest destabilizing |Delta|, naming the
% perturbation |wDelta| in this case.
% Both *constant, complex* and *real, linear dynamic* perturbations can be
% constructed, following the same procedure, incorporating the weighting
% function |wu|.  Start by computing the norm of the weighted complementary
% sensitivity function.
[wTiNorm,wfreq] = norm(wu*H.Ti,inf,0.001);
M = freqresp(wu*H.Ti,wfreq);
[U,S,V] = svd(M);
%%
% Constant, complex perturbation is constructed with the same steps, based
% on SVD
wDelta = -1/S(1,1)*V(:,1)*U(:,1)'; 
[norm(wDelta)  1/wTiNorm]
%%
% Real, linear, dynamic perturbation is constructed in the same manner,
% using |cnum2sys| on the elements from the SVD construction.
v1hat = cnum2sys(V(1,1),wfreq);
v2hat = cnum2sys(V(2,1),wfreq);
ubar1hat = cnum2sys(conj(U(1,1)),wfreq);
ubar2hat = cnum2sys(conj(U(2,1)),wfreq);
wDeltaHat = -1/S(1,1)*[v1hat;v2hat]*[ubar1hat ubar2hat];
pole(wDeltaHat)
%%
% The constructed perturbations are of the same magnitude.
[norm(wDeltaHat,inf) norm(wDelta)]
%%
% Both result in closed-loop poles on the imaginary axis at the frequency
% associated with the peak value of the |wu*Ti|.
pole(feedback(P*(eye(2)+wu*wDelta),C))
%%
pole(feedback(P*(eye(2)+wu*wDeltaHat),C))

%% Conclusions
% Using the |svd| command, the smallest destabilizing uncertainty for an
% input-multiplicative uncertainty model is constructed.  It is initially
% constructed as a complex matrix, but can also be constructed as a real,
% linear, dynamic system, using |cnum2sys|.  If the input-multiplicative
% uncertainty model is frequency-dependent, the same process can be used,
% simply incorporating the weighting function

%% Attribution
% Copyright 2016-17, Andy Packard.  This work is licensed under the Creative
% Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. To
% view a copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to
% Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.

%% File Information
disp(mfilename)