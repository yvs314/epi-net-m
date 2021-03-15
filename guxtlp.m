%% Costate Dynamics for Net-SIR, with control
    % \dot{\lambda}_s}=\beta(\lambda_s - \lambda_z)
    %   \odot\mathbf{D}(\bar{u})\mathbf{A}z \\
    % \dot{\lambda_{z}} =-e^{r_{1}t}cN + 
    %   \beta{\mathbf{A}}^\top (\lambda_s - \lambda_z)
    %   \odot\bar{u}\odot s 
    %   + \gamma\lambda_z}
%input x as ode45's solution at t; lax=[la_s;la_z] 
%output as dlax=[las;laz]


function dlax = guxtlp(u,x,t,lax,beta,gamma,A,r1,c,NN) 
n = size(u,1); %u is [1 \times n], 1 scalar control per node    
s = x(1:n); z = x(n+1 : end); %get susceptible & infected components
las = lax(1:n); laz = lax(n+1:end); %likewise for costate
szd = las - laz; %precompute \lambda_s - \lambda_z, or let JIT do it?

dlas= beta*(szd .* (1-u).* (A * z) );
dlaz= -exp(r1*t)*c*NN ... 
        + beta* ((A') * szd) .* (1-u) .* s ...
        + gamma * laz; 
dlax = [dlas;dlaz];
end

