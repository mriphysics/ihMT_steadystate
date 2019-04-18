function [g,D] = gauss_lineshape(T2b,f)

g = (T2b/(sqrt(2*pi)))*exp(-0.5*(2*pi*f*T2b).^2);

%%% D is the local fields term, related to second moment. Morrison 1995
D = 1/(sqrt(3)*T2b);

end