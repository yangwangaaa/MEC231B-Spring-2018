function [xnew] = bicycle_discrete_nonlin_dyn(t,x,u,params)
% discrete dynamics through euler discretization
deltaT = 0.01;
xnew = x + deltaT*bicycle_model_nonlin_dyn(t,x,u,params); 
end