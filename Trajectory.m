%% Scirpt by Shuang Gao March 27, 2020
% Function: Illustrating the figures. 
% Change 1, YS, 2020-04-13: re-wrote with switch expression 
% Change 2, YS, 2020-04-15: added line color parameter
% Change 3, YS, 2020-04-15: moved the invariants outside the switch 
% Change 4, YS, 2020-05-26: rm figure export; can well export from the main
% Change 5, YS, 2020-05-27: add proper node labels, the city names
% Change 6, YS, 2020-05-27: fix figure tiling

function [outFig] = Trajectory(zApp,t,lColor,plotTitle,nodeLabels)
% Assume zApp rows are state, columns are time steps, zApp has the same
% dimension at t
Y   = diag(1:size(zApp,1))* ones(size(zApp,1),size(t,2));

outFig = plot3(t,Y,zApp(:,:),'-','MarkerSize',1.5,'LineWidth',0.1, 'Color',lColor);

%yticks: space the nodes' curves
yticks(1:size(zApp,1)); %set ticks to integers; effectively, city IDs
%yticklabels:a row-vec of strings, e.g.['ATL';'BOS';'CLT';'DEN']
yticklabels(nodeLabels);

%ylabel not necessary?
xlabel('Time (days)'),zlabel('Compartment Population');
view(19,31);
grid on;
title(sprintf('%s', plotTitle));

end
