function cost = evaluateCost(A,E,C,F,Lmat)
% compute nX from A
nX = size(A,1);
% compute nY from C
nY = size(C,1);
% compute N from size of Lmat
N = size(Lmat,2)/nY-1;
% use createKfromLmat to obtain matrices AK,BK,CK,DK
[AK,BK,CK,DK] = createKfromLmat(Lmat,nX,nY,N);
% use combineProcessEstimator to create G
G = combineProcessEstimator(A,E,C,F,AK,BK,CK,DK);
% use norm to obtain cost.
cost = nomr(G,2);
end