%% q2
G = tf([-5.8 -3.2 9],[0.02 0.98 4.8 9]);
normTol = 1e-5;
[A,B] = norm(G,inf,normTol);
[A abs(freqresp(G,B))]
bodemag(G, frd(freqresp(G,B),B), 'ro')