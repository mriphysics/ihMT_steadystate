%%% function to generate Ernst signal curve for given flip (deg), TR, R1
function sig = ernst(flip,TR,R1)
e1=exp(-TR*R1);
sig = sind(flip).*(1-e1)./(1-e1.*cosd(flip));
end