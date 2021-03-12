%% Objective function's running cost, in time-vectorized form
    %L(\cdot,t)= e^{r_1t}cz^(t)\bullet N 
    %+ e^{r_2t}\frac{l}{2}u^{2}(t)\bullet N

function Ltxu = Ltxu(uArr,xArr,tArr,r1,r2,c,l,N)
n = size(uArr,1); zArr = xArr(n+1:end,:);
%time exponents with coeffs c or l/2 for each t \in tArr
myexp1 = arrayfun(@(t) exp(r1*t) * c,tArr); 
myexp2 = arrayfun(@(t) exp(r2*t) * (l/2),tArr);
%time exponents multiplied by infecteds `z`, resp. controls `u`
rcz = repmat(myexp1,n,1) .* zArr;
rcu = repmat(myexp2,n,1) .* uArr;
%transform the costs into absolute forms
rc = (rcz + rcu) .* N; %these are per-node costs
Ltxu = sum(rc); %sum each column to get overall running cost at each time
end

