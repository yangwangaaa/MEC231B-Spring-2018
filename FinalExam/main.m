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
% At the line 3, we use norm calculator with a setting of accuracy to get
% the maximum frequency response and its corresponded frequency, thus in
% the line 4, we have A = abs(freqresp(G,B)). At the line 5, we use frd to
% compare the system G with a newly defined system, where at the frequency B,
% the frequency response is exactly the maximum frequency response of G.
% Finally, we see the superposition of bodegram of two these system at
% frequnecy B.

q2;

%% 3
% At the line 4, we use random complex number delta to generate D. We go to
% the file cnum2sys.m, as we have seen in the homework and the slide, the
% behavior of cnumsys.m is that the infinity norm of this system is the
% same as the complex number generated and the associated frequency is
% wBar.  At the line 5, we calculate the frequency 
% response of D and to compare it with delta (absolutely the same!). At the
% line 6, we verify that wBar is indeed the frequency. At the line 7, we
% generate the bode plot of D, we find out the gain is constant and equals
% to 20*log(delta).
q3;
