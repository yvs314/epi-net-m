%% Scirpt by Shuang Gao March 27, 2020
% Function: Illustrating the figures. 
% Change 1, YS, 2020-04-13: re-wrote with switch expression 

function Trajectory(zApp,t,titleName,figName)
% Assume zApp rows are state, columns are time steps, zApp has the same
% dimesion at t

hfig3DTra = figure('name','3-D trajectory');
nodes   = diag(1:size(zApp,1))* ones(size(zApp,1),size(t,2));
%plot3(t,nodes,zApp(:,:),'.'); 
plot3(t,nodes,zApp(:,:),'-','MarkerSize',1.5,'LineWidth',0.1, 'Color','b');
xlabel('time (s)'),ylabel('node'), zlabel('value')
view(19,31)
grid on

switch nargin
    
    case 3
        title(sprintf('%s', titleName));
    
    case 4
      title(sprintf('%s', titleName));  
      set(hfig3DTra, 'Units','inch','Position',[2 8 4 3],'PaperUnits', 'inch','PaperPosition',[0 8 4 3] ); 
      print(sprintf('%s', figName),'-depsc');

    %case 5   
end

end
