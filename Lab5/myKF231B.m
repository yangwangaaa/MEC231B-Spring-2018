function [xk1k,Sxk1k, ykk1] = myKF231B(xkk1,Sxkk1,A,B,C,D,E,F,Swk,uk,yk)

tmp1 = (C*Sxkk1*C' + F*Swk*F');
tmp2 = tmp1\(yk- C*xkk1 -D*uk);
% state update
% ------------
% hat{x}_{k+1|k} = xk1k
xk1k = A*xkk1 + B*uk + A*Sxkk1*C'*tmp2;
% variance update
% ---------------
tmp3 = tmp1\(C*Sxkk1);
Sxk1k = A*(Sxkk1 - Sxkk1*C'*tmp3)*A' + E*Swk*E';
ykk1 = C*xkk1 + D*uk;
end