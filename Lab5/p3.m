%% lab5: vehicle dynamics p3.m
% symbolic linearization discrete nonlinear equation
clear; close all; clc;
%%
syms vx vy wz X Y psi Fx delta real
syms m Jz a b A B C dt real

q=[vx;vy;wz;X;Y;psi];
u=[Fx;delta];

af = atan((vy+a*wz)/vx) - delta;
ar = atan((vy-b*wz)/vx);
Ff = @(alpha_f) m*(a/(a+b))*A*sin(C*atan(B*alpha_f));
Fr = @(alpha_r) m*(b/(a+b))*A*sin(C*atan(B*alpha_r));
f = [vy*wz - (2*Ff(af)*sin(delta))/m + Fx/m;
     -vx*wz + (2/m)*(Ff(af)*cos(delta) + Fr(ar)); 
     (2/Jz)*(a*Ff(af)*cos(delta)-b*Fr(ar)); 
     vx*cos(psi)-vy*sin(psi); 
     vx*sin(psi)+vy*cos(psi);
     wz];
y = [vx; 2*(Ff(af)*cos(delta)+Fr(ar))/m; wz];
qk1 = q+dt*f;
Amat = jacobian(qk1,q);
Bmat = jacobian(qk1,u);
Cmat = jacobian(y,q);
Dmat = jacobian(y,u);

q0=[10;0;0;0;0;0];
u0=[0;0];
params_sym = [m; Jz; a; b; A; B; C; dt];
params_val = [2237; 5112; 1.46; 1.55; -6.8357; 0.0325; 238.9874; 0.01];
Aval = double(subs(Amat,[q; u; params_sym], [q0; u0; params_val]));
Bval = double(subs(Bmat,[q; u; params_sym], [q0; u0; params_val]));
Cval = double(subs(Cmat,[q; u; params_sym], [q0; u0; params_val]));
Dval = double(subs(Dmat,[q; u; params_sym], [q0; u0; params_val]));

syms vxk vyk wzk Xk Yk psik Fxk ∆k real
qk = [vxk; vyk; wzk; Xk; Yk; psik];
uk = [Fxk; deltak];

matlabFunction(Amat, Bmat, Cmat, Dmat, 'File',...
'getLinearDynamics_at_k', 'Vars', {q, u, params_sym});

save('bicycle_lindyn.mat','Aval','Bval','Cval','Dval',...
'Amat','Bmat','Cmat','Dmat','q','u');