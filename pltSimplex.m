%% plot: Draw a trajectory of one node in the (z,s) simplex
% given two series of susceptible `s` and infected `z` fractions
% OPT: nodeTag is optional, omit to suppress the title 
function [outFig] = pltSimplex(s_in, z_in, nodeTag)   
        
outFig = newplot;
fplot(@(x) 1-x, [0,max(z_in)]);
hold on;
scatter(z_in,s_in);
hold off;
    
%title and axis labels
if(nargin>=3)
    title(nodeTag);
end
xlabel('z');
ylabel('s');
    
    % legend
  %  lgd = legend('s','Susceptible','Recovered');
   % legend('boxoff');
end