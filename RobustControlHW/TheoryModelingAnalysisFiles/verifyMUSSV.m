%% Verifying the bounds from |MUSSV|
% Through elementary, but detailed calculations, the lower and upper bounds
% from |mussv| can be verified.
%
% UC Berkeley, ME C231B/EECS C220C, Spring 2017

%% Create block structure and compute |mussv| of a random matrix
blk = [-2 0;2 0];  % 1-by-1 real, repeated twice; 1-by-1 complex, repeated twice
M = complex(randn(4),randn(4));
[bnds,muinfo] = mussv(M,blk);
%%
% Check that upper bound is less than or equal to norm(M)
bnds(1)<=norm(M)
%%
% Check that lower bound is greater than all real eigenvalues 
eVal = eig(M);
realeVal = eVal(imag(eVal)==0);
max([realeVal;0])<=bnds(2)
%%
% Extract verification parameters
[VDelta,VSigma,VLmi] = mussvextract(muinfo);

%% Check perturbation and lower bound
% Verify that Delta has correct dimension
size(VDelta)
%%
% Examine Delta, and confirm it has the proper
% structure.
VDelta
%%
% M*VDelta should have an eigenvalue at 1 (so det(I-M*VDelta)=0)
eig(M*VDelta)
%%
% Lower bound should be equal to 1/norm(VDelta)
[bnds(2) 1/norm(VDelta)]

%% Verify upper bound, using certificates from |VSigma|
% Fields |GLeft|, |GMiddle| and |GRight| should be diagonal, with nonzero
% real entries associated with *real* uncertainties.  These matrices should
% be "equal" although they may be different dimensions when |M| is
% nonsquare.
VSigma.GLeft
%%
VSigma.GMiddle
%%
VSigma.GRight
%%
% Fields |DLeft| and |DRight| should be block-diagonal, and structured
% consistent with the structure defined by |blk|.  Specifically, for any
% matrix |Delta| with the structure defined by |blk|, it should be that
% |Delta*DLeft = DRight*Delta|.   Both |DLeft| and |DRight| should be
% square and invertible, and "essentially" equal, although dimensions can
% be different.  If |blk| contains any nonsquare full blocks, then the
% dimensions of |DLeft| and |DRight| will be different.
VSigma.DLeft
%%
VSigma.DRight
%%
% The upper bound is verified by a norm test on a matrix involving the
% matrix |M|, the upper bound (|bnds(1)|, and the scalings.  Form the
% matrix in 3 steps.
TL = (eye(4)+VSigma.GLeft^2)^-0.25;
TM = (1/bnds(1))*VSigma.DLeft*M/VSigma.DRight - sqrt(-1)*VSigma.GMiddle;
TR = (eye(4)+VSigma.GRight^2)^-0.25;
%%
% The upper bound is certified if NORM(TL*TM*TR)<=1.  Check this
norm(TL*TM*TR)

%% Verify upper bound, using certificates from |VLmi|
% Fields |Grc|, |Gcr| should be block-diagonal, with nonzero square matrix
% component entries associated with *real* uncertainties.  These square
% block entries should be Hermitian, and the entries in |Grc| should equal
% those in |Gcr|.  
VLmi.Gcr
%%
VLmi.Grc
%%
% Fields |Dc| and |Dr| should be block-diagonal, hermtian, positive-definite
% and structured consistently with the structure defined by |blk|. 
% For any matrix |Delta| with the structure defined by |blk|, it should be
% that |Delta*Dr = Dc*Delta|.   The matrices |Dc| and |Dc| should be
% "essentially" equal, although dimensions will be different if |blk|
% contains any nonsquare full blocks.
VLmi.Dc
%%
eig(VLmi.Dc)
%%
VLmi.Dr
%%
% The upper bound is verified by a semidefinite test on a Hermitian matrix
% formed from the matrix |M|, the upper bound (|bnds(1)|, and the scalings.
% Construct the matrix, and verify that it is negative semidefinite.
T = M'*VLmi.Dr*M - bnds(1)^2*VLmi.Dc + sqrt(-1)*(VLmi.Gcr*M-M'*VLmi.Grc);
eig(T)

%% Verfication for only square, complex blocks is easier
% Repeat the steps for a problem involving only square complex blocks.
blk = [3 0;2 2;1 1];
M = complex(randn(6),randn(6));
[bnds,muinfo] = mussv(M,blk);
[VDelta,VSigma,VLmi] = mussvextract(muinfo);

%% Check perturbation and lower bound
% Verify that Delta has correct dimension
size(VDelta)
%%
% Examine Delta, and confirm it has the proper structure.
VDelta
%%
% M*VDelta should have an eigenvalue at 1 (so det(I-M*VDelta)=0)
eig(M*VDelta)
%%
% Lower bound should be equal to 1/norm(VDelta)
[bnds(2) 1/norm(VDelta)]

%% Verify upper bound, using certificates from |VSigma|
% Fields |GLeft|, |GMiddle| and |GRight| will all be zero
VSigma.GLeft
%%
isequal(VSigma.GLeft,VSigma.GMiddle)
%%
isequal(VSigma.GRight,VSigma.GMiddle)
%%
% Fields |DLeft| and |DRight| will be block-diagonal, and structured
% consistent with the structure defined by |blk|, and equal.
VSigma.DLeft
%%
isequal(VSigma.DRight,VSigma.DLeft)
%%
% With no |G|-scalings, the upper bound is simply the "DMDinv" norm.
[norm(VSigma.DLeft*M/VSigma.DRight) bnds(1)]

%% Verify upper bound, using certificates from |VLmi|
% Fields |Grc|, |Gcr| are zero
VLmi.Gcr
%%
isequal(VLmi.Grc,VLmi.Gcr)
%%
% Fields |Dc| and |Dr| should be block-diagonal, hermitian, positive-definite
% and structured consistently with the structure defined by |blk|, and equal. 
VLmi.Dc
%%
eig(VLmi.Dc)
%%
isequal(VLmi.Dr,VLmi.Dc)
%%
% With no |G|-scalings, the upper bound is simply the "DMDinv" norm.  In
% the LMI certificates, the matrix square roots must be taken of the |D|
% matrices. 
[norm(sqrtm(VLmi.Dr)*M/sqrtm(VLmi.Dc)) bnds(1)]

%% Upper-Bound quality
% By making approximations, the convex optimization for the upper-bound can
% be solved quickly, at the expense of bound quality.  The |'f'| option
% forces a faster, but less complete minimization, usually resulting in 
% lower quality (ie., larger) upper bounds.  Test this option.  Run this
% section several times to see the effect.
Delta = [-3 0;-2 0;1 0;1 0];
M = crandn(7,7);
tic; bndDefault = mussv(M,Delta); tDefault = toc;
tic; bndFaster = mussv(M,Delta,'f'); tFaster = toc;
disp(['Default upper bound: ' num2str(bndDefault(1)) ', Time: ' num2str(tDefault)]);
disp(['Faster upper bound: ' num2str(bndFaster(1)) ', Time: ' num2str(tFaster)]);

%% Ordering of |MU| based on order of |DELTA| set
% Create 3 different block structures, each of superset of the preceeding
DeltaA = [-1 0;-1 0;1 0];
DeltaB = [-1 0;1 0;1 0];   % contains DeltaA
DeltaC = [-1 0;2 2];       % contains DeltaB
M = crandn(3,3);
% Compute bounds for MU
bndA = mussv(M,DeltaA);
bndB = mussv(M,DeltaB);
bndC = mussv(M,DeltaC);
% Based on Delta, muA <= muB <= muC.  The upper bounds should satisfy this
% ordering as well.
[bndA(1) bndB(1) bndC(1)]

%% Conclusions
% The lower and upper bounds for |MUSSV| can be certified using all of the
% information returned by |mussv| in the 2nd output argument.  The lower
% bound is certified by simply exhibiting a perturbation matrix with the
% corect structure that causes singularity.  The lower bound is then the
% reciprocal of that matrix.  The upper bound is certified using the
% S-procedure, which amounts to the existence of "D" and "G" matrices with
% the appropriate structure (consistent with the uncertainty structure)
% that together with |M|, atisfy an inequality.   While the formulae in
% notes and literature are simple, in practice they are slighly more
% complicated due to the possibility of non-square full blocks, which make
% the dimension bookkeeping a bit clumsy.

%% File Information
disp(mfilename)