%% Infected, absolute, in time-vectorized form
    %Z(x(t))= z(t)\bullet N 

function LZt = LZt(xArr,N)
n = size(xArr,1)/2; 
ZArr = xArr(n+1:end,:).*N; %per-node infected
%sum each column to get total infected at each time
LZt = sum(ZArr); 
end

