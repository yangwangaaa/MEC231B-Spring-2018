%% Plant/Kalman filter variance simulation template
% This template can be used for performance assessment of a Kalman filter,
% as applied to a known plant model.  The simulation loop only computes the
% state estimate error variances and various Kalman gains (but not a state
% estimate).  For this reason, a specific disturbance sequence and specific
% initial condition are not needed - only the variances of these quantities
% need to be specified.  The error variances and Kalman gains are stored in
% struct arrays.

%% Attribution
% ME C231B, UC Berkeley, Spring 2018

%% Notation
% This file aids in the simulation of the time-varying Kalman filter for N
% times steps, from i=0, to i=N-1.  So for each matrix in the problem data,
% there needs to be N values (from i=0 to i=N-1).  We use 3-d arrays, where
% the 3rd dimension is N.  The variance of the initial condition is denoted
% by Sx0. The computed variances are notated as: Sxii1 means S^x_{i|i-1},
% Syii1 means S^y_{i|i-1} and Sxii means S^x_{i|i}.

%% Define sequence length
N = 

%% Create System
% The Kalman Filter code below assumes that the state-space arrays (A, E,
% C, F) represent time-varying dynamics, and hence should be defined as 3-d
% arrays.
nX = 
nW = 
nY = 
arrayA =   % nX-by-nX-by-N
arrayE =   % nX-by-nW-by-N
arrayC =   % nY-by-nX-by-N
arrayF =   % nY-by-nW-by-N

%% Create/Declare variance of disturbance and initial condition
% The variance of the "w" sequence is allowed to be time-varying, and
% should be defined as a 3-d array.
arraySW =  % nW-by-nW-by-N
Sx0 =      % nX-by-nX

%% Initialize KF states with appropriate values
Sxii1 = Sx0;

%% Compute performance of filter without real-time implementation
% The performance, in terms of estimation error variance can be computed
% directly form the system model, disturbance model, and initial variance
% in state.
emptyx = [];
emptyB = [];
emptyu = [];
emptyy = [];
cArray = cell(N,1);  % used to allocate STRUCT arrays
arraySigma = struct('Sxi1i',cArray,'Sxii',cArray,'Syii1',cArray);
arrayGains = struct('Li',cArray,'Hi',cArray,'Gi',cArray);
for i=0:N-1
   iMatlab = i+1;
   Ai = arrayA(:,:,iMatlab);
   Ei = arrayE(:,:,iMatlab);
   Ci = arrayC(:,:,iMatlab);
   Fi = arrayF(:,:,iMatlab);
   Swi = arraySW(:,:,iMatlab);
   [~,Sxi1i,~,Sxii,Syii1,Li,Hi,Gi,~] = ...
      KF231B(emptyx,Sxii1,Ai,emptyB,Ci,Ei,Fi,Swi,emptyu,emptyy);
   arraySigma(iMatlab).Sxi1i = Sxi1i;
   arraySigma(iMatlab).Sxii = Sxii;
   arraySigma(iMatlab).Syii1 = Syii1;
   arrayGains(iMatlab).Li = Li;
   arrayGains(iMatlab).Hi = Hi;
   arrayGains(iMatlab).Gi = Gi;
   % Shift the error-variance estimate so that when loop-index i advances,
   % the initial condition for this variance is correct.
   Sxii1 = Sxi1i;
end
