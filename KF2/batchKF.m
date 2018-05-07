function [LkBatch,VkBatch,ekVar] = batchKF(Amat,Emat,Cmat,Fmat,Sigma0,SigmaW,k)
% Batch Kalman filter formulation from ME C231B, UC Berkeley, Spring 2018

nx = size(Amat,1);
% Initialize the 6 batch matrices
[calA,Psi,N,Gam,R,S] = LOCALbuildBatch(nx);

% Fill in batch matrices in order to estimate x0 to xk, as linear combinations
% of y from 0,1,...,k-1.  This needs the state-space matrices defined on
% i = [0,1,...,k-1]
for i=0:k-1
   Ai = Amat(:,:,i+1);  % i+1 is needed because Matlab indices start at 1
   Ei = Emat(:,:,i+1);
   Ci = Cmat(:,:,i+1);
   Fi = Fmat(:,:,i+1);
   [calA,Psi,N,Gam,R,S] = LOCALbuildBatch(calA,Psi,N,Gam,R,S,Ai,Ei,Ci,Fi);
end
% The next line assumes that the variance of w_k is constant across k,
% given by the single matrix SigmaW.  The KRON command simply repeats this
% matrix in a block diagonal form, as in the lecture slides.
SigmaBarW = kron(eye(k),SigmaW);
SigmaZ = Psi*Sigma0*Psi' + Gam*SigmaBarW*Gam';
SigmaP = R*Sigma0*R' + S*SigmaBarW*S';
SigmaZP = Psi*Sigma0*R' + Gam*SigmaBarW*S';
LkBatch = SigmaZP/SigmaP;
VkBatch = Psi - LkBatch*R;
ekVar = SigmaZ - LkBatch*SigmaZP';

function [calA,Psi,N,Gam,R,S] = LOCALbuildBatch(calA,Psi,N,Gam,R,S,A,E,C,F)
if nargin==1
   nx = calA;  % first (and only) argument is NX
   calA = eye(nx); Psi = eye(nx); Gam = zeros(nx,0);
   N = zeros(nx,0); S = zeros(0,0); R = zeros(0,nx);
else
   % Determine k, ny, nw, nx from dimensions of A, E, C, F and (eg) S
   [ny,nx] = size(C);
   nw = size(E,2);
   k = size(S,1)/ny;
   % Use recursions to build new, larger matrices.  Use temporary
   % variable for the updated calA and N, since the original values are
   % needed to define a few other of the updated matrices.
   newCalA = A*calA;
   newN = [A*N E];
   Psi = [Psi;newCalA];
   Gam = [[Gam zeros((k+1)*nx,nw)];newN];
   R = [R;C*calA];
   S = [S zeros(k*ny,nw);C*N F];
   calA = newCalA;
   N = newN;
end