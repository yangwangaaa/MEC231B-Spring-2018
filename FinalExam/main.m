%% Final Exam - ME C231B Jun Zeng

addpath('./allfiles/');

%% Vehicle Navigation

%% 1.a
% Nothing to hand in

%% 1.b
q1b;

%% 1.c
q1c;

%% 1.d
q1d;

%% 1.e
q1e;

%% Robust Control

%% 2
% In the line 3, we use norm calculator with a setting of accuracy to get
% the maximum frequency response and its corresponded frequency, thus in
% the line 4, we have A = abs(freqresp(G,B)). In the line 5, we use frd to
% compare the system G with a newly defined system, where at the frequency B,
% the frequency response is exactly the maximum frequency response of G.
% Finally, we see the superposition of bodegram of two these system at
% frequnecy B.

q2;

%% 3
% 

q3;
