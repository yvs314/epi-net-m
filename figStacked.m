%% aux: Draw Stacked I+S+R for 1 city
% take a node's s_i(t), i_i(t), and r_i(t) series and output a stacked plot
% with i_i(t) as the *first* series: thus its shape is not distorted
 
%YS: need to get used to working with figure as an output object
function [outFig] = figStacked(s_in, i_in, r_in, t_f)
    outFig = area([i_in; s_in; r_in]');
    outFig(1).FaceColor=[0.85 0.325 0.098];%Infected are reddish
    outFig(2).FaceColor=[0 0.447 0.741];%Susceptible are blueish
    outFig(3).FaceColor=[0.9 0.9 0.9];%Removed are greyish/whatever
 
end