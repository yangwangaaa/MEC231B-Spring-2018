function [y] = bicycle_measurements(t,x,u,params)

% function defining the nonlinear measurement model for vehicle dynamics
% y = h(x,u)
%% unpacking parameters
% vehicle parameters
m = params.m;
Jz = params.Jz;
a = params.a;
b = params.b;

% tire model parameters
A = params.A;
B = params.B;
C = params.C;

%% extracting states and inputs
[vx, vy, wz, X, Y, psi] = deal(x(1), x(2), x(3), x(4), x(5), x(6));

Fx = u(1);
delta = (u(2));

%% measurement
% slip angles
if vx~=0
    af = atan((vy+a*wz)/vx) - delta;
    ar = atan((vy-b*wz)/vx);
else
    af = 0;
    ar = 0;
end
y = [vx;2*(Ff(af)*cos(delta)+Fr(ar))/m;wz];
%% tire forces
    function Ff_val = Ff(alpha_f)
        Ff_val = m*(a/(a+b))*A*sin(C*atan(B*alpha_f));
    end
    function Fr_val = Fr(alpha_r)
        Fr_val = m*(b/(a+b))*A*sin(C*atan(B*alpha_r));
    end
end