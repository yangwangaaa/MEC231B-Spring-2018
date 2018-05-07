function [K,L,H,G,nIter,estErrVar] = formSteadyStateKF(A,E,C,F,W,TS)
% TS is specified sample time.  Use -1 to indicate a discrete-time
% system with unspecified sample time.
% Get dimensions
nX = size(A,1);
nY = size(C,1);
nW = size(W,1);
% Note: Your code should check for consistency of dimensions across
% A, E, C, F and W, but this is not necessary in this asignment.
% Of course, in general, please use good programming practices.

% Initialize value for Sxkk1.  Convergence will occur for
% any positive-definite initialization.
Sxkk1 = eye(nX);

% Initialize ``previous'' values for L, H, G
prevL = zeros(nX,nY);
prevH = zeros(nX,nY);
prevG = zeros(nW,nY);
iterConverged = false;
nIter = 0;
while ~iterConverged
   nIter = nIter + 1;
   [~,Sxk1k,~,Sxkk,~,Lk,Hk,Gk,~] = ...
      KF231B([],Sxkk1,A,[],C,E,F,W,[],[]);
   if norm(Lk-prevL)<1e-8 && norm(Hk-prevH)<1e-8 && ...
         norm(Gk-prevG)<1e-8 && norm(Sxk1k-Sxkk1)<1e-8
      iterConverged = true;
      L = Lk;
      H = Hk;
      G = Gk;
      K = ss(A-L*C,L,[eye(nX);eye(nX)-H*C;-G*C], [zeros(nX,nY);H;G], TS);
      K.OutputGroup.xkk1 = 1:nX;
      K.OutputGroup.xkk = nX+1:2*nX;
      K.OutputGroup.wkk = 2*nX+1:2*nX+nW;
      % The steady-state variance of {hat x}_{k|k} - x_k is Sxkk.
      % This is a useful output to understand.  Using the
      % Chebychev inequality, you can (probabilistically)
      % bound the potential error in the state-estimate.
      estErrVar = Sxkk;
   else
      Sxkk1 = Sxk1k;
      prevL = Lk;
      prevH = Hk;
      prevG = Gk;
   end
end