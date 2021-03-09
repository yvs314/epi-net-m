%% Optimal Control, in the time-vectorized form 
% u=\frac{\beta}{e^{r_2t}l }(\lambda_z - \lambda_s )
% \mathbf{D}(s\oslash N)\mathbf{A}z

% tArr is a vector of time points, e.g. 0:T; 
%[x,lax]Arr are [2n \times |tArr|], values of x and lax at tArr
function u2 = utxla(tArr,xtArr,laxtArr,umin,umax, beta,l,A,r2,iN)
n = size(xtArr,1)/2; %x is [1 \times 2n], per-node susceptible then per-node infs
s = xtArr(1:n,:);  z = xtArr(n+1 : end,:); %carve up the state
las = laxtArr(1:n,:); laz = laxtArr(n+1:end,:); %likewise for costate

myexp = arrayfun( @(t) exp(-r2*t), tArr); %make the exponents for t\in tArr
time_dep_term = repmat((beta/l) * myexp, n,1);%[n \times |tArr|]

u1 = time_dep_term .* (laz - las) .* s .* iN .* (A * z);
u2 = min(max(u1,umin),umax); %apply the bounds
end

