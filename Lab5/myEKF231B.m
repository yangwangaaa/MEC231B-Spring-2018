function [xk1k,Sxk1k,ykk1] = myEKF231B(xkk1,Sxkk1,A,B,C,E,F,Swk,uk,yk,params)
%%
Sykk1 = C*Sxkk1*C' + F*Swk*F';
tmpK = Sxkk1*C'/Sykk1;
Sxkk = Sxkk1-tmpK*C*Sxkk1;
Gk = Swk*F'/Sykk1;
ykk1 = bicycle_measurements(0,xkk1,uk,params);
ek=yk-ykk1;
xkk = xkk1 + tmpK*ek;
p = [params.m; params.Jz; params.a; params.b; params.A ;params.B; ...
      params.C; params.dt];
[Akk,~,~,~] = getLinearDynamics_at_k(xkk,uk,p);
fk = (xkk+params.dt*bicycle_model_nonlin_dyn(0,xkk,uk,params) )+ E*Gk*ek;
xk1k = fk + E*Gk*ek;
tmp = Akk*tmpK*F*Swk*E';
Sxk1k = Akk*Sxkk*Akk' + E*(Swk-Gk*F*Swk)*E'-tmp-tmp';
end