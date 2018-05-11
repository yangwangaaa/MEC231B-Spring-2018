normTol = 1e-5;
delta = complex(randn,randn);
wBar = exp(2*randn);
D = cnum2sys(delta,wBar);
[freqresp(D,wBar) delta]
[norm(D,inf,normTol) abs(delta)]
bode(D)