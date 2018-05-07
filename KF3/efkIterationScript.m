%% Plant/Extended Kalman filter simulation template
% This template can be used for open-loop simulations involving two
% separate components: plant and Extended Kalman filter. 

%% Attribution
% ME C231B, UC Berkeley, Spring 2018
N = 1000;   % should be sufficient number of steps
Delta = 0.025;  % "sample" time for 1st-order Euler approximation

%% Create plant as a dynamical system
nX = 3;
nW = 2;
nY = 1;
rho0 = 1.125; gamma = 7200; g = 9.8; 
fH = @(x) FallingObjectDynEq(x,Delta,rho0,gamma,g);
arrayE = repmat([0 0;Delta 0;0 0],[1 1 N]);
hH = @(x) FallingObjectOutEq(x,Delta,rho0,gamma,g);
arrayF = repmat([0 1],[1 1 N]);
partialFX = @(x) FallingObjectDynEqJac(x,Delta,rho0,gamma,g);
partialHX = @(x) FallingObjectOutEqJac(x,Delta,rho0,gamma,g);

%% Create/Declare variance of disturbance and initial condition
SW = [(0.1*g)^2 0;0 35^2];
arraySW = repmat(SW,[1 1 N]);
Sx0 = diag([1000 1e4 8e9]);
m0 = [3e4;-1900;3.6e5];

%% Initialize KF states with appropriate values
Sxii1 = Sx0;
xii1 = m0;

%% Create a specific initial condition and noise sequence
w = sqrtm(SW)*randn(nW,N);
x0 = m0 + 2*(rand(nX,1)-0.5).*diag(sqrt(Sx0));

%% Simulate the system/KF one step at a time
y = zeros(nY,N);
x = zeros(nX,N);
x(:,1) = x0;
xEii1 = zeros(nX,N); xEii1(:,1) = m0;
xEii = zeros(nX,N);
arraySxii = zeros(nX,nX,N);
for i=0:N-1
   iMatlab = i+1;
   E = arrayE(:,:,iMatlab);
   F = arrayF(:,:,iMatlab);
   Swi = arraySW(:,:,iMatlab);
   % Get y(i) from system model, using x(i) and w(i)
   y(:,iMatlab) = hH(x(:,iMatlab)) + F*w(:,iMatlab);
   % Get estimates of x(i+1|i) using y(i) from KF (EKF231BToFinish)
   [xi1i,Sxi1i,xii,Sxii,Syii1] = EKF231B(xEii1(:,iMatlab),Sxii1,...
      partialFX,partialHX,E,F,fH,hH,Swi,y(:,iMatlab));
   xEii(:,iMatlab) = xii;
   xEii1(:,iMatlab+1) = xi1i;
   arraySxii(:,:,iMatlab) = Sxii;
   % Get x(i+1) from system model, using x(i) and w(i)
   x(:,iMatlab+1) = fH(x(:,iMatlab)) + E*w(:,iMatlab);
   % Shift the error-variance estimate so that when loop-index i advances,
   % the initial condition for this variance is correct.
   Sxii1 = Sxi1i;
   if x(1,iMatlab+1)<0
      % Stop simulation when Altitude drops below 0
      break
   end
end
% Set N to last computed iteration (recall that loop terminates when
% altitude reaches 0).
N = iMatlab;

%% Plot estimation-error and SQRT(ErrVariance)
% For each state, plot the estimation error (since this is just a
% simulation, we have the actual state, as well as the estimate).  Also
% plot the sqrt of the diagonal entries of the Sigma^x_{k|k} matrix, as
% that is an "approximate" confidence interval of the size of the error
% that we could compute in a real implementation.
tVec = Delta*(0:N-1);
for i=1:3
   subplot(3,2,2*(i-1)+1);
   errVar = arraySxii(i,i,1:N);
   plot(tVec,x(i,1:N)-xEii(i,1:N),...
      tVec,sqrt(errVar(:)),'r',tVec,-sqrt(errVar(:)),'r')
   % legend('EstError','UpBndEst','LowBndEst','location','SouthEast');
   ylabel(['Error, x_' int2str(i)]);
   xlabel('Time, s')
   subplot(3,2,2*i);
   plot(tVec,x(i,1:N),tVec,xEii(i,1:N))
   ylabel(['x_' int2str(i)]);
   xlabel('Time, s')
end

