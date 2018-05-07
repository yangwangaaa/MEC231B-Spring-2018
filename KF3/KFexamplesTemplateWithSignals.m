%% Plant/Kalman filter/Controller simulation template
% This template can be used for closed-loop simulations involving three
% separate components: plant, Kalman filter, and controller.  The
% controller has its own states, as well as a reference input, and is
% implemented (in this script) as a dynamical system with function handles.
%
% It may be useful to extend this script to allow for a "true" plant, and a
% "modeled" plant.  The Kalman Filter (and controller) have access to (and 
% are based on) the "modeled" plant equations, but the actual evolution of
% the [disturbance;control] -> output would be based on the "true" plant.
% This type of arrangment would allow you to study the robustness of the
% overall control system (KF and controller) to unknown variations in the
% plant behavior.

%% Attribution
% ME C231B, UC Berkeley, Spring 2018


%% Notation
% This file aids in the simulation of the time-varying Kalman filter for N
% times steps, from i=0, to i=N-1.  Hence for each matrix of problem data,
% there needs to be N values (from i=0 to i=N-1).  We use 3-d arrays, where
% the 3rd dimension is N.  The variance of the initial condition is denoted
% by Sx0. The computed variances are notated as: Sxii1 means S^x_{i|i-1},
% Syii1 means S^y_{i|i-1} and Sxii means S^x_{i|i}.  The signal estimates
% are notated as follows: xii1 means xHat_{i|i-1}, xii means xHat_{i|i} and
% wii means wHat_{i|i}

%% Define sequence length
N = 

%% Create System
% The Kalman Filter code below assumes that the state-space arrays (A, E,
% C, F) represent time-varying dynamics, and hence should be defined as 3-d
% arrays.
nX = 
nU = 
nW = 
nY = 
arrayA =   % nX-by-nX-by-N
arrayB =   % nX-by-nU-by-N
arrayE =   % nX-by-nW-by-N
arrayC =   % nY-by-nX-by-N
arrayF =   % nY-by-nW-by-N

%% Create/Declare variance of disturbance and initial condition (mean too)
% The variance of the "w" sequence is allowed to be time-varying, and
% should be defined as a 3-d array.
arraySW =  % nW-by-nW-by-N
Sx0 =      % nX-by-nX
m0 =       % nX-by-1


%% Create representation of controller
% Controller has internal states (n_i), a reference input (r_i), and
% measures all aspects of the state estimates.   It is implemented using
% function handles as
%     OutputEquation:   ui = outCHan(ni,ri)
%     StateEquation:    ni1 = dynCHan(ni,ri,xi1i,Sxi1i,xii,Sxii)
%     InitialCondition: n0 = specified initial condition
outCHan = @(n,r) FILL-IN-FUNCTION
dynCHan = @(n,r,xhat1,Sx1,xhat,Sx) FILL-IN-FUNCTION
n0 = 

%% Initialize KF states with appropriate values
Sxii1 = Sx0;
xii1 = m0;

%% Initialize Controller states with appropriate values
ni = n0;

%% Create a specific initial condition and noise sequence
% Under ideal circumstances, this should be consistent with the statistical
% assumptions made in the previous code cell.  When studying robustness,
% namely how the filter performance degrades as assumptions are not met, it
% may be useful to create an initial condition and noise sequence which is
% not consistent with the assumptions
wSeq =     % nW-by-1-by-N
x0 =       % nX-by-1

%% Create a specific reference input
arrayR =   % nR-by-1-by-N

%% Simulate the system/KF one step at a time
ni = n0;
y = zeros(nY,N);
xSeq = zeros(nX,N);
xSeq(:,1) = x0;
xiiSeq = zeros(nX,N);
for i=0:T-1
   iMatlab = i+1;
   Ai = arrayA(:,:,iMatlab);
   Bi = arrayB(:,:,iMatlab);
   Ei = arrayE(:,:,iMatlab);
   Ci = arrayC(:,:,iMatlab);
   Fi = arrayF(:,:,iMatlab);
   Swi = arraySW(:,:,iMatlab);
   wi = wSeq(:,iMatlab);
   xi = xSeq(:,iMatlab);
   % Get r(i), value of reference input
   ri = arrayR(:,iMatlab);
   % Get u(i) from controller output map
   ui = outCHan(ni,ri);
   % Get y(i) from system model, using x(i) and w(i)
   y(:,iMatlab) = Ci*xi + Fi*wi;
   % Get estimates of x(i+1|i) using u(i) and y(i) from KF
   [xi1i,Sxi1i,xii,Sxii,Syii1,Lk,Hk,Gk,wkk] = ...
      KF231B(xii1,Sxii1,Ai,Bi,Ci,Ei,Fi,Swi,ui,y(:,iMatlab));
   % Save the state estimate xHat_{i|i}
   xiiSeq(:,iMatlab) = xii;
   % Get x(i+1) from system model, using x(i), u(i) and w(i)
   xSeq(:,iMatlab+1) = Ai*xi + Bi*ui + Ei*wi;
   % Get n(i+1) from controller state map
   ni1 = dynCHan(ni,ri,xi1i,Sxi1i,xii,Sxii);
   ni = ni1;
   % Save any signals or variances that will be used in later calculations.
   %
   % Shift the error-variance estimate so that when loop-index i advances,
   % the initial condition for this variance is correct.
   Sxii1 = Sxi1i;
   xii1 = xi1i;
end

