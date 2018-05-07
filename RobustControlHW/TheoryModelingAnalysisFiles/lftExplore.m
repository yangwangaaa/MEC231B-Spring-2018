%% lftExplore.m
% Investigate how LFTs are used to represent uncertain objects with many
% uncertain atoms.
%
% UC Berkeley, ME C231B/EECS C220C, Spring 2017

%% Create a simple linear-fractional transformation
% The transformation 3 + 4x is just linear, so there will be no fractional
% term.
x = ureal('x',0);
G = 3 + 4*x;
[M,D] = lftdata(G);
M
D

%% Equivalent expression
% Same value, expressed differently, yields different M
H = 3 + 2*x*2;
[M,D] = lftdata(H);
M

%% Linear in two variables
y = ureal('y',0);
H = 5 + 2*x*3 - 4*y;
[M,D] = lftdata(H);
M
D

%% Simplest quadratic
% Quadratic (in 1 variable) requires 2 "copies" of the UREAL object
G = x*x;
[M,D] = lftdata(G);
M
D

%% Unexpected
% Expressions are created sequentially, and can have "extra" unnecessary
% copies of x
G = x + x*x;
[M,D] = lftdata(G);
M
D
%%
% Redefining, in an equivalent expression results in the expected number of
% copies
G = x*(1+x);
[M,D] = lftdata(G);
M
D

%% Meaning of |AutoSimplify|
% The |AutoSimplify| flag controls what reduction methods for simplification
% are attempted when LFTs are created.
x.AutoSimplify
%%
% Setting it to 'full' directs numerical model-reduction-like schemes to
% be used in order to simplify LFTs
x.AutoSimplify = 'full';
G = x + x*x;
[M,D] = lftdata(G);
M
D
%%
% Note the fewer copies of the UREAL (now 2, as expected), but numerical
% methods have clearly been used, resulting in non-integer (in this simple
% example) entries in M

%% Conclusions
% LFTs are used to represent dependencies on uncertain elements.  It is
% (unfortunately) easy to create LFTs with many "extra" copies of uncertain
% elements, as automated reduction schemes are difficult to perfect.
% Research in efficient representations continues on (for example, at NASA
% and throughout Europe aerospace community) and often exploits
% problem-dependent structure to obtain "simple" LFT models.

%% File Information
disp(mfilename)