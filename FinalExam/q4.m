%% q.4

%% 4.a
P = tf(1,[1,-1]);
C = tf([5.8 9],[0.04 1 0]);
isstable(feedback(P,C))

%% 4.b
G = P/(1+P*C);

%% 4.c
normTol = 0.001;
%smallest value of delta
bound = 1/norm(G,inf,normTol)
error = 0.1;

% Verification
% unstable
delta = complex(bound-error/sqrt(2),bound-error/sqrt(2));
istable(feedback(P+delta,C))
pole(feedback(P+delta,C))



