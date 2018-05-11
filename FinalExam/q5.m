%% q5

%% 5.a
P = tf(1,[1 -1]); 
C = tf([5.8 9],[0.04 1 0]);
S = 1/(1+P*C);
Wp = tf([0.667 3],[1 0.003]);
normTol = 0.001;
norm(Wp*S,inf,normTol)<=1 

%% 5.b
figure
bodemag(S,'r--',1/Wp,'bo')
legend('S','1/Wp')

%% 5.c
% Based on the results of 4.g, the system is stable, moreover we can verify
% it by using the function robstab.
delta = ultidyn('delta',[1 1],'bound',1);
Pu = P*(1+0.4*delta);
Snew = 1/(1+Pu*C);
[stabmarg,wcu] = robstab(Snew);

%% 5.d
[wcg,wcu] = wcgain(Wp/(1+Pu*C))
% specific stable linear system
wcu.delta

%% 5.e
% We can see the worst case gain from the intersection of 1/Wp and the
% sensitivity of the worse case Pu (refer to the small gain theorem), the
% intersection frequency is larger than 1, which confirm the stability and
% the worse case gain seen in 5.c and 5.d.  
Pu_worst = P*(1+0.4*wcu.delta);
Sworst = 1/(1+Pu_worst*C);
figure
bodemag(S,'r-',Sworst,'go',1/Wp,'y.')
legend('Sensitivity for P','Sensitivity for the worst Pu','1/Wp')


%% 5.f
figure
step(S,P,Pu_worst,4)
legend('S','P','Worst case Pu')

