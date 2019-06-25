% Function to compute Gaussian lineshape and associated local field
% strength, w_loc.
% 
% [g,w_loc] = gauss_lineshape(T2s,f)
% 
% INPUTS:       T2s = T2 of semisolid compartment, seconds
%               f   = frequency at which function is to be evaluated (can
%               be vector)
% OUTPUTS:      g   = lineshape value, s
%               w_loc = local field term (rad/s)
%               
% Shaihan Malik (c), King's College London, April 2019

function [g,w_loc] = gauss_lineshape(T2s,f)

g = (T2s/(sqrt(2*pi)))*exp(-0.5*(2*pi*f*T2s).^2);

%%% w_loc is the local fields term, related to second moment. Morrison 1995
w_loc = 1/(sqrt(3)*T2s);

end