%% q.4

%% 4.a
P = tf(1,[1,-1]);
C = tf([5.8 9],[0.04 1 0]);
isstable(feedback(P,C))

%% 4.b
G = -C/(1+P*C);

%% 4.c
normTol = 0.001;
fprintf('The smallest absolute value of Delta\n')
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
delta = ucomplex('delta',0,'Radius',0.1454);
%we use the smallest value found in 4.c
[stabmarg,destabunc,report] = robuststab(feedback(P+delta,C))

%% 4.e

