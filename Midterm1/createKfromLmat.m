function [AK,BK,CK,DK] = createKfromLmat(Lmat,nX,nY,N)
if size(Lmat) == [nX,nY*(N+1)]
    fprintf('');
else
    fprintf('ERROR!');
end
CK = Lmat(:,1:nY*N);
DK = Lmat(:,nY*N+1:nY*(N+1));
AK = [zeros(nY*(N-1),nY),eye(nY*(N-1));zeros(nY,nY),zeros(nY,nY*(N-1))];
BK = [zeros(nY*(N-1),nY);eye(nY,nY)];
end