function [xk1k,Sxk1k,xkk,Sxkk,Sykk1] = ...
   KF231BtoBuildOn(xkk1,Sxkk1,A,B,C,E,F,Swk,uk,yk)
% Implements one step of the Kalman filter, using the notation in slides at
% the end of Estimation3.  The input arguments are:
%    xkk1 = xhat_{k|k-1}
%    Sxkk1 = Sigma^x_{k|k-1}
%    A, B, C, E, F: state matrices at time k
%    Swk = variance in zero-mean disturbance w_k
%    uk = deterministic control input u_k
%    yk = measurement y_k
% The output arguments are:
%    xk1k = xhat_{k+1|k}
%    Sxk1k = Sigma^x_{k+1|k}
%    xkk = xhat_{k|k}
%    Sxkk = Sigma^x_{k|k}
%    Sykk1 = Sigma^y_{k|k-1} 
% The input signals, xkk1, uk, yk may all be empty matrices, implying that
% the function will only updated the error variances, and will not provide
% any state estimates (so xk1k and xkk will be returned empty as well).

% ME C231B, UC Berkeley, Spring 2018

% Group the variance calculations together.  Recall that the evolution of
% these quantities does not depend on measurements.
Sykk1 = C*Sxkk1*C' + F*Swk*F';
tmpK = Sxkk1*C'/Sykk1;
Sxkk = Sxkk1-tmpK*C*Sxkk1;
Gk = Swk*F'/Sykk1;
tmp = A*tmpK*F*Swk*E';
Sxk1k = A*Sxkk*A' + E*(Swk-Gk*F*Swk)*E'-tmp-tmp';

% Calculate state estimations, if measurements and previous estimates are
% supplied. 
if ~isempty(xkk1) && ~isempty(yk)
   ek = yk - C*xkk1;
   xkk = xkk1 + tmpK*ek;
   xk1k = A*xkk + E*Gk*ek;
   % If there is a control input, apply that to the evolution of the state
   % estimate.
   if ~isempty(B)
      xk1k = xk1k + B*uk;
   end
else
   % Copy empty inputs over to empty outputs
   xkk = xkk1;
   xk1k = xkk1;
end
