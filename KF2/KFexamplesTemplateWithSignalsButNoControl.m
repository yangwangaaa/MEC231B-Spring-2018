%% KFexamplesTemplateWithSignalsButNoControl
% This template can be used for open-loop simulations involving two
% separate components: plant and Kalman filter. 
%
% It may be useful to extend this script to allow for a "true" plant, and a
% "modeled" plant.  The Kalman Filter has access to (and is based on) the
% "modeled" plant equations, but the actual evolution of the disturbance
% -> output would be based on the "true" plant.  This type of arrangement
% allows you to study the robustness of the KF to unknown variations in the
% plant behavior.

%% Attribution
% ME C231B, UC Berkeley, Spring 2018

%% Create plant
T =  % sequence length, data runs from 0:(T-1)
nX = 
nW = 
nY = 
arrayA = 
arrayE = 
arrayC = 
arrayF = 

%% Create/Declare variance of disturbance and initial condition
arraySW =
Sx0 = 
m0 = 

%% Initialize KF states with appropriate values
Sxii1 = Sx0;
xii1 = m0;

%% Create a specific initial condition and noise sequence
% Under ideal circumstances, this should be consistent with the statistical
% assumptions made in the previous code cell.  When studying robustness,
% namely how the filter performance degrades as assumptions are not met, it
% may be useful to create an initial condition and noise sequence which is
% not consistent with the assumptions
wSeq = 
x0 = 

%% Simulate the system/KF one step at a time
emptyB = [];  % this template is for no control signal
emptyu = [];  % this template is for no control signal
y = zeros(nY,T);
xSeq = zeros(nX,T);
xSeq(:,1) = x0;
xEii1 = zeros(nX,T); xEii1(:,1) = m0;
xEii = zeros(nX,T);
for i=0:T-1
   iMatlab = i+1;
   Ai = arrayA(:,:,iMatlab);
   Ei = arrayE(:,:,iMatlab);
   Ci = arrayC(:,:,iMatlab);
   Fi = arrayF(:,:,iMatlab);
   Swi = arraySW(:,:,iMatlab);
   wi = wSeq(:,iMatlab);
   xi = xSeq(:,iMatlab);
   
   % Get y(i) from system model, using x(i) and w(i)
   y(:,iMatlab) = Ci*xi + Fi*wi;
   % Get estimates of x(i+1|i) using y(i) from KF
   [xi1i,Sxi1i,xii,Sxii,Syii1] = KF231B(xEii1(:,iMatlab),Sxii1,...
      Ai,emptyB,Ci,Ei,Fi,Swi,emptyu,y(:,iMatlab));
   xEii(:,iMatlab) = xii;
   xEii1(:,iMatlab+1) = xi1i;
   % Get x(i+1) from system model, using x(i) and w(i)
   xSeq(:,iMatlab+1) = Ai*xi + Ei*wi;
   % Shift the error-variance estimate so that when loop-index i advances,
   % the initial condition for this variance is correct.
   Sxii1 = Sxi1i;
end

