function Gsys = combineProcessEstimator(A,E,C,F,AK,BK,CK,DK)
Lx = size(A,1);
Lz = size(AK,1);
ANew = [A,zeros(Lx,Lz);BK*C,AK];
BNew = [E;BK*F];
CNew = [-DK*C+eye(size(DK*C,1),size(DK*C,2)),-CK];
DNew = -DK*F;
Gsys = ss(ANew, BNew, CNew, DNew, -1);
end