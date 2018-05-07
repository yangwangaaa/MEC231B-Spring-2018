%% Robust Control HW
addpath('./TheoryModelingAnalysisFiles')

%% cnum2sys.m
deltaNum = complex(randn,randn);
w = exp(3*randn);
deltaSys = cnum2sys(deltaNum,w);
zpk(deltaSys)
[norm(deltaSys,inf) abs(deltaNum)]
freqresp(deltaSys,w)

%% 1. Uncertain Analysis: theory

%% Summary for verifyMUSSV.m
% In order to use the high-level command mussv, the input argument needs to
% follow the specific rule stated in the official documentation to create
% the block structure. The output of the mussv, the muinfo gives us a
% comprehensive result about the upper bound information (the newlin/young
% method consisting finding the scailing beta, D, and G, the combination of
% which is smaller than 1 would certify the upperbound; the semidefinite
% verification of the upperb ouds is connected to the lecture slides, which
% means that the G term is included to shift the complex disk to cover the 
% real term in the delta) and lower bound information (VDelta, which is for the
% lower bound, the product of VDelta and M would make the det(I-M*VDelat =
% 0),
% and the system will be unstable, thus we could find the structured delta 
% for the system to be unstable)? If the block is purely complex, then the 
% general formulations of the upperbounds and lowerbounds can be simplified,
% which is consistent with the lecture notes. Additionally, for three blocks,
% where one contains another, the largest block would give the smallest minimum singular value,
% and the reciprocal of which is the uppper bound, and thus the biggest upperbound among the three.

%% Summary for DestabilizingPerturbation.m
% The DestabilizingPerturbationSolution.m show us the process. we firstly
% find out the upper bound of uncertainty block, assume it as beta and prove the infinity norm
% of the plant is supperior than 1 divided beta.

%% Summary for UncertainModeling.m
% The uncertain objects |ureal| and |ultidyn| can be used to form uncertain
% matrices and uncertain state-space models.  Many of the usual operations
% and manipulations for numerical matrices and state-space models are
% implemented, and behave as expected.

%% Summary for lftExplore.m
% Matlab provide lftdata.m function to help us decompose the uncertainty
% LFT problems. LFTs are used to represent dependencies on uncertain elements.  It is
% (unfortunately) easy to create LFTs with many "extra" copies of uncertain
% elements, as automated reduction schemes are difficult to perfect.

%% 2. Uncertainty Modeling and Analysis

%% Summary for WorstCaseAnalysisIntroduction.m
% Illustration of the use of some tools for modeling uncertain systems, and analyzing their worst-case behavior.

%% Summary for stateSpaceParameterUncertainty.m
% This file use uncertainty to create LFT problems. Model has independent forces at each mass as inputs.  All states are
% listed as outputs.  $A$ matrix depends on mass, damping and stiffness
% matrices, while $B$ only depends on mass matrix.  Since outputs consist
% of all states, $C$ and $D$ matrices are simple.

%% Summary for dcmotor demo1.m
% This example shows how to use uncertain objects in Robust Control
% Toolbox(TM) to model uncertain systems. A follow-on to this demo
% will show how to automate robustness calculations using the robustness 
% analysis tools. The closed-loop system is reasonably robust despite
% significant fluctuations in the plant DC gain.  This is a desirable, and
% common characteristic of a properly designed feedback system.

%% Summary for rsrpmu.m
% An uncertain plant model is a lightly-damped, second-order system with
% parametric uncertainty in the denominator coefficients and significant
% frequency-dependent unmodeled dynamics beyond 6 radians/seconds. Then the
% Stability Robustness with Respect to Modeled Uncertainty is done with a
% sensitivity analysis. Finally a robust analysis for all values of |k|,
% |delta| is carried out.

%% Summary for FlightControlExWithSSV.m
% The open-loop model includes the rigid-body aircraft dynamics. 
% These dynamics contain uncertainties in the aerodynamic coefficients. 
% The open-loop model also includes dynamics for the rudder and aileron
% actuators.  The actuator models also include dynamic uncertainty.
% The closed-loop is uncertain due to the real parametric and dynamic 
% uncertainty in the plant model. An uncertain frequency response for 
% the closed-loop system is computed with the UFRD command. Robust
% stability of the closed loop is assessed using the ROBUSTSTAB command.

%% Summary for mimoMotivate.m
% This is an example showing coupled sensitivities and the final part
% clearly illustrates how the inverse-based controller completely 
% decouples the closed-loop response of the plant to output disturbances.

%% Summary for mimoMotivateResolveMU.m
% By implementing sensitivity in each chanel and use the |mussv| analysis,
% the code confirms what was observed in simulation.  The 
% excellent nominal performance, along with excellent robust stability
% properties do not guarantee robustness of performance.


