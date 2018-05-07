function [xk1k,Sxk1k,xkk,Sxkk,Sykk1] = ...
   EKF231B(xkk1,Sxkk1,dfdx,dhdx,E,F,fhandle,hhandle,Swk,yk)
% Implements one step of the Extended Kalman filter, using the notation in
% slides

C = dhdx(xkk1);
Sykk1 = C*Sxkk1*C' + F*Swk*F';
ek = yk - hhandle(xkk1);
tmpK = Sxkk1*C'/Sykk1;
xkk = xkk1 + tmpK*ek;

Sxkk = Sxkk1-tmpK*C*Sxkk1;
Gk = Swk*F'/Sykk1;
A = dfdx(xkk);

xk1k = fhandle(xkk) + E*Gk*ek;
tmp = A*tmpK*F*Swk*E';
Sxk1k = A*Sxkk*A' + E*(Swk-Gk*F*Swk)*E'-tmp-tmp';

% UC Berkeley, Spring 2018, ME C231B

