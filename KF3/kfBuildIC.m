function IC = kfBuildIC(G,P,K,M)
% this function generate the interconnected state space of the augmented
% system based on the noise processor state space, plant state space, and
% the kalman filter state space.

% sampling time:
tsG = G.Ts;
tsP = P.Ts;
tsK = K.Ts;
ts = 0;

if tsG == tsP == tsK
    ts = tsG;
end

% extract the input matrices
AG = G.A;
BG = G.B;
CG = G.C;
DG = G.D;
AP = P.A;
BP = P.B;
CP = P.C;
DP = P.D;
AK = K.A;
BK = K.B;
CK = K.C;
DK = K.D;

% concatenate to form the new system
AIC = [AG zeros(size(AG,1),size(AP,2)) zeros(size(AG,1),size(AK,2));
    BP*CG AP zeros(size(AP,1),size(AK,2));
    BK*DP*CG BK*CP AK];
BIC = [BG; BP*DG; BK*DP*DG];
CIC = [-M*DK*DP*CG eye(size(M,1), size(CP,2))-M*DK*CP -M*CK];
DIC = -M*DK*DP*DG;

IC = ss(AIC, BIC, CIC, DIC, ts);

end