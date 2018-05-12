%% q.4

%% 4.a
P = tf(1,[1,-1]);
C = tf([5.8 9],[0.04 1 0]);
isstable(feedback(P,C))

%% 4.b
G = -C/(1+P*C);

%% 4.c
normTol = 0.001;
fprintf('The smallest absolute value of Delta calculated by the small gain thm\n')
[InfNorm,freq] = norm(G,inf,normTol);
bound = 1/InfNorm
%verification: to get the smallest delta
M = freqresp(G,freq);
[U,S,V]=svd(M);
deltamin = -1/S(1,1)*V(:,1)*U(:,1)';
pole(feedback(P-deltamin,C))
%one pure imaginary pool appear at this critical value 

%% 4.d
% We see that it can tolerate up to 100% of the modeled uncertainty, the
% results found here is eactly as we have seen in 4.c
deltamin = ucomplex('delta',0,'Radius',0.1454);
%we use the smallest value found in 4.c
[stabmarg,destabunc,report] = robuststab(feedback(P+deltamin,C))

%% 4.e
deltamin = -1/S(1,1)*V(:,1)*U(:,1)';
deltanew = cnum2sys(deltamin,freq);
pole(feedback(P+deltanew, C))
%We can see clearly that there are two poles on the imaginary
% axis, thus the dynamic system generated in 4.e is unstable.

%% 4.f
deltanew = ultidyn('delta', [1 1], 'Bound', norm(deltamin));
[stabmarg,destabunc,report] = robuststab(feedback(P+deltanew,C))

%% 4.g
% G is the sensitivity function of the negative feedback system of P,C
P = tf(1,[1 -1]);
C = tf([5.8 9],[0.04 1 0]);
G_new = P*C/(1+P*C);

%% 4.h
[Inf_norm, freq] = norm(G_new, inf, normTol);
M = freqresp(G_new, freq);
[U,S,V] = svd(M);
deltamin = -1/S(1,1)*V(:,1)*U(:,1)';
fprintf('norm of delta\n')
disp(norm(deltamin))
fprintf('norm of delta calculated by small gain thm\n')
disp(1/Inf_norm)
P_tilde = P*(1+deltamin);
pole(feedback(P_tilde, C))

%% 4.i
deltanew = ultidyn('delta', [1 1], 'Bound', norm(deltamin));
uncertainSys = feedback((P*(1+deltanew)),C);
[stabmarg,destabunc,report] = robuststab(uncertainSys)
% By refering to the output in stabmarg and destabunc, the
% system is marginally unstable, which meets the result found in 4.h