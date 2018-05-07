%% MU explanation: large, coupled sensitivities in MIMO problems
% See also: |mimoMotivate.m| and |mimoSolve.m|.
%
% UC Berkeley, ME C231B/EECS C220C, Spring 2017

%% Nominal Model
% Constant N-by-N real matrix.  The plant is *not* diagonal.  It represents
% an idealized manufacturing process, such as paper pulp spraying, where
% many nearly identical actuators contribute in a spatially invariant
% manner to a final result.
N = 5;
G = toeplitz([1;.75;.40;zeros(N-3,1)],[1 .75 .40 zeros(1,N-3)]);
G = ss(G)

%% Performance/Uncertainty interaction
% Since the plant |G| is nondynamic, any/all frequency-dependence will
% arise in the interaction between the uncertainty weightings (|Wunc|) and
% the performance weightings (|Wp|).  A parameter |BWsep| will be used to
% quantify the separation between:
%
% * the desired bandwidth of sensitivity reduction, [0 1/BWsep], and
% * the frequency range for which plant uncertainty exceeds 100%, [BWsep inf).
%
BWsep = 4.5;
wp = makeweight(100,1/BWsep,0.33);
wu = makeweight(0.33,BWsep,20);

%% Uncertainty in each InputChannel
% The plant model will have independent uncertainty in all input channels.
% This would be typical in almost every situation.  
Wunc = ss(zeros(N,N));
for i=1:N
   Wunc(i,i) = wu;
end

%% Output Sensitivity Performance Objective
Wp = wp*eye(N);

%% Inverse-Based Controller
Kinv = ss(0,1,1,0)*inv(G.d);

%% RobustPerformance MUSSV calculation
% Create closed-loop system with uncertainty and performance input/outputs
systemnames = 'G Kinv Wunc Wp';
iN = int2str(N);
inputvar =  ['[w{' iN '};d{' iN '}]'];
outputvar = '[Kinv;Wp]';
input_to_G = '[Wunc+Kinv]';
input_to_Kinv = '[-d-G]';
input_to_Wunc = '[w]';
input_to_Wp = '[d+G]';
M = sysic;
%%
% Compute frequency response, since |mussv| only works on a
% frequency-by-frequency basis.
Mg = frd(M,logspace(-3,3,200));
%%
% Create |mussv| block structure definition, for 
%
% * N, 1-by-1 actuator uncertainties
% * N-by-N "performance" block
%
blk = [ones(N,2);N N];
%%
% Compute structured singular value with |mussv|, and plot bounds.
bnds = mussv(Mg,blk,'as');
semilogx(bnds)
xlabel('Frequency')
ylabel('RobustPerformance MU bounds')

%% Conclusions
% The |mussv| analysis confirms what was observed in simulation.  The
% excellent nominal performance, along with excellent robust stability
% properties *do not* guarantee robustness of performance.  In this example,
% with the robust performance mu value equal to approximately 2.8, it is
% proven that *merely* 37% (1/2.8) of the modeled uncertainty can lead
% to performance degradation of 2.8 times that what was desired.

%% Attribution
% Copyright 2016-17, Andy Packard.  This work is licensed under the Creative
% Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License. To
% view a copy of this license, visit
% http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to
% Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.

%% File Information
disp(mfilename)