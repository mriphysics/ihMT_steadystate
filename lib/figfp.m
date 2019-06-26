function [h] = figfp(a)

% figfp does exactly as figure, but also clears the screen
% SJM: also set paper position mode

h = figure(a);
set(h,'paperpositionmode','auto');

clf(a);
if nargout==0
    h=[];
end