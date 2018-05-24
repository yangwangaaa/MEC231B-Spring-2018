%% Midterm (Jun Zeng)
close all
options = optimoptions(@fminunc,'Display','none');

%% 4.g
load('midtermModel1.mat')
NList = [1, 2, 3, 4, 5, 8, 12, 15];
cost1 = zeros(1,length(NList));
for i = 1:length(NList)
    N = NList(i);
    f = @(Lmat) evaluateCost(A,E,C,F,Lmat);
    LmatInit = zeros(nX,(N+1)*nY);
    [Lopt, optCost] = fminunc(f, LmatInit, options);
    [AK,BK,CK,DK] = createKfromLmat(Lopt,nX,nY,N);
    G1 = combineProcessEstimator(A,E,C,F,AK,BK,CK,DK);
    cost1(i) = norm(G1,2);
end
plot(NList,cost1,'ro');

W = eye(size(E, 2));
TS = -1;
[K,L,H,G,nIter,estErrVar] = formSteadyStateKF(A,E,C,F,W,TS);
fprintf('Estimation xkk')
K.OutputGroup.xkk

[K,L,H,G,nIter,estErrVar] = formSteadyStateKF(A,E,C,F,W,TS);
GsysK = combineProcessEstimator(A,E,C,F,A-L*C,L,eye(nX)-H*C, H);
costKF = norm(GsysK, 2);

% calculate covar
covarKF = covar(GsysK, eye(nW))
Gsys15 = combineProcessEstimator(A,E,C,F,AK,BK,CK,DK);
covarL15 = covar(Gsys15, eye(nW))

cost2 = zeros(length(NList),1);
for i = 1:length(NList)
    cost2(i) = costKF; 
end
figure
plot(NList, cost1, 'c*')
hold on 
plot(NList, cost2,'r')
title('2-norm versus N, for several estimator designs')


%% 4.h
% variables that should exist:
% process model (A, E, C, F) consistent dimensions
% N, positive integer, dictates complexity of estimator
% nX and nY, consistent with (A,C)
% Create the function handle for the cost function. The decision variable
% is the matix Lmat = @(Lmat) evaluateCost(A,E,C,F,Lmat);
load('testCaseSingleProcessModel.mat')
NList = [1, 2, 3, 4, 5, 8, 12, 15];
cost1 = zeros(1,length(NList));
costList2 = zeros(1,length(NList));
CovarList = zeros(1,length(NList));
for i = 1:length(NList)
    N = NList(i);
    f = @(Lmat) evaluateCost(A,E,C,F,Lmat);
    LmatInit = zeros(nX,(N+1)*nY);
    [Lopt, optCost] = fminunc(f, LmatInit, options);
    [AK,BK,CK,DK] = createKfromLmat(Lopt,nX,nY,N);
    G1 = combineProcessEstimator(A,E,C,F,AK,BK,CK,DK);
    cost1(i) = norm(G1,2);
end

% Steady State Testing
N = 15;
TS = -1;
W = eye(size(E,2));
[K,L,H,G,nIter,estErrVar] = formSteadyStateKF(A,E,C,F,W,TS);
disp(estErrVar)

%% 4.j

%% Robust Design
load('midtermMultipleModels.mat')
NList = [1,2,3,4,5,8,12];
cost1 = zeros(1,length(NList));
for i = 1:length(NList)
   N = NList(i);
   f = @(Lmat) evaluateWorstCostManyModels(A(:,:,i),E(:,:,i),C(:,:,i),F(:,:,i),Lmat);
   LmatInit = zeros(nX,(N+1)*nY);
   [Lopt, optCost] = fminunc(f, LmatInit, options);
   [AK,BK,CK,DK] = createKfromLmat(Lopt,nX,nY,N);
   G = combineProcessEstimator(A(:,:,i),E(:,:,i),C(:,:,i),F(:,:,i),AK,BK,CK,DK);
   cost1(i) = norm(G,2);
end
figure
plot(NList,cost1);
xlabel('N')
ylabel('Estimator Design for different N')

%% Robustness of the estimator

figure
for i = [3,6,9]
    Ai = A(:,:,i);
    Ei = E(:,:,i);
    Ci = C(:,:,i);
    Fi = F(:,:,i);
    W = eye(2);
    [K,L,H,G,nIter,estErrVar] = formSteadyStateKF(Ai,Ei,Ci,Fi,W,-1);
    cost1 = [];
    for j = 1:11
        Aj = A(:,:,j);
        Ej = E(:,:,j);
        Cj = C(:,:,j);
        Fj = F(:,:,j);
        M = [zeros(nX),eye(nX),zeros(nX,nW)];
        ICA = [Aj zeros(size(Aj,1),size(K.A,2));K.B*Cj K.A];
        ICB = [Ej;K.B*Fj];
        ICC = [eye(size(M,1))-M*K.D*Cj -M*K.C];
        ICD = -M*K.D*Fj;
        IC = ss(ICA,ICB,ICC,ICD,-1);
        cost1 = [cost1 norm(IC,2)];
    end
    plot(cost1);
    hold on 
end

N = 12;
cost2 = zeros(11,1);
for i = 1:11
   f = @(Lmat) evaluateWorstCostManyModels(A(:,:,i),E(:,:,i),C(:,:,i),F(:,:,i),Lmat);
   LmatInit = zeros(nX,(N+1)*nY);
   [Lopt, optCost] = fminunc(f, LmatInit,options);
   [AK,BK,CK,DK] = createKfromLmat(Lopt,nX,nY,N);
   G = combineProcessEstimator(A(:,:,i),E(:,:,i),C(:,:,i),F(:,:,i),AK,BK,CK,DK);
   cost2(i) = norm(G,2);
end
plot(1:11,cost2);

N = 12;
costgamma = zeros(11,1);
for i = 1:size(A,3)
    [AK,BK,CK,DK] = createKfromLmat(Lopt,nX,nY,N);
    Gsys = combineProcessEstimator(A(:,:,i),E(:,:,i),C(:,:,i),F(:,:,i),AK,BK,CK,DK);
    cost = norm(Gsys, 2);
    costgamma(i) = cost;
end
plot(1:size(A,3), costgamma, '-*')
legend('\alpha_{ij} (j = 3)','\alpha_{ij} (j = 6)','\alpha_{ij} (j = 9)','\alpha_{ii}','\gamma_i')

%% Comments
% At first, we can see that the time variant kalman filter is better than the fixed
% time kalman filter and alternate one. Let's look at \alpha_{ij}, for a fixed value of j, the
% covariance of error in the
% conbination of jth kalman filter with i process model will reach the
% minimal value when i = j. We also plot \alpha_{ii} to reveal this comparaison.

