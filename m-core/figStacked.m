%% aux: Draw Stacked I+S+R for 1 city
% take a node's s_i(t), z_i(t), and r_i(t) series and output a stacked plot
% i_i(t) is the *first* series, to keep its shape intact
% OPT: nodeTag is optional, omit to suppress the title 
function [outFig] = figStacked(X_ts, s_in, z_in, r_in,nodeTag)   
    
    outFig = area(X_ts,[z_in; s_in; r_in]');
 
    %set the colors
    outFig(1).FaceColor=[0.85 0.325 0.098];%Infected are reddish
    outFig(2).FaceColor=[0 0.447 0.741];%Susceptible are blueish
    outFig(3).FaceColor=[0.941 0.863 0.51];%Removed are buff 
    outFig(3).FaceAlpha=0.7; % make the removed even lighter
    
    %title and axis labels
    if(nargin>=5)
        title(nodeTag);
    end
    xlabel('Time (days)');
    ylabel('Compartment Population');
    
    % legend
    lgd = legend('Infected','Susceptible','Recovered');
    title(lgd,'Compartment Types');
    legend('boxoff');
end