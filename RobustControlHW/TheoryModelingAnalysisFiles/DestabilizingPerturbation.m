%% Constructing the smallest destabilizing perturbation (at Plant input).
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

%% Compute Hinf norm of Ti
% Ti is complementary sensitivity at plant input

%% Compute Frequency Response Matrix of Ti at frequency where peak occurs

%%
% Find the smallest plant input-multiplicative uncertainty (a constant,
% complex matrix of dimension 2-by-2) which causes instability

%% Perturbed closed-loop system
% Form perturbed plant, using input-multiplicative form and the
% perturbation matrix.

%% Verify Instability of perturbed closed-loop system
% Compute perturbed closed-loop system using |loopsens|, and check poles

%% Connection between unstable pole, and frequency where peak of Ti occurs
% Show connection between the two frequencies

%% (Extra Credit!) Use |cnum2sys| to create real-dynamic uncertainty
% Rather than accepting a complex matrix as a destabilizing element, use
% the command |cnum2sys| to convert the complex matrix uncertainty into a
% real-valued, linear, dynamic system that leads to instability.  The Hinf
% norm of the dynamic system should match the maximum singular value of the
% original destabilizing complex matrix.

%% Frequency-dependent uncertainty
% Adopt a frequency-dependent uncertainty model, using a scalar weighting
% function, |wu|.  Take uncertainty model to be |P*(eye(2) + wu(s) Delta)|,
% using a first order |wu| with:
%
% * DC gain of 0.35 (35% uncertainty at low frequency)
% * gain of 1 at 40 rads/TimeUnit (100% uncertainty at 40)
% * High-Frequency gain of 10 (1000% uncertainty at high frequency)
% 
% Use |makeweight| for the construction.

%%
% Mimic/repeat procedure to find smallest destabilizing |Delta|, naming the
% perturbation |wDelta| in this case.
% Both *constant, complex* and *real, linear dynamic* perturbations can be
% constructed, following the same procedure, incorporating the weighting
% function |wu|.  Start by computing the norm of the weighted complementary
% sensitivity function.

%%
% Constant, complex perturbation is constructed with the same steps, based
% on SVD

%%
% Real, linear, dynamic perturbation is constructed in the same manner,
% using |cnum2sys| on the elements from the SVD construction.

%%
% The constructed perturbations are of the same magnitude.

%%
% Both result in closed-loop poles on the imaginary axis at the frequency
% associated with the peak value of the |wu*Ti|.

%% Conclusions
% Using the |svd| command, the smallest destabilizing uncertainty for an
% input-multiplicative uncertainty model is constructed.  It is initially
% constructed as a complex matrix, but can also be constructed as a real,
% linear, dynamic system, using |cnum2sys|.  If the input-multiplicative
% uncertainty model is frequency-dependent, the same process can be used,
% simply incorporating the weighting function
