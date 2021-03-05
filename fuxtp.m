%% State Dynamics for Net-SIR, with control
    %\dot{s}=-\beta \mathbf{D}(\bar{u}\odot s)\mathbf{A}z \\
    %\dot{z}=\beta \mathbf{D}(\bar{u}\odot s)\mathbf{A}z - \gamma z
%input x=[s;z] output as dx=[s;z]

%tilde ~ instead of t since state eqs are autonomous (costate ain't)
function dx = fuxtp(u,~,x,beta,gamma,A) 
n = size(u,1); %u is [1 \times n], 1 scalar control per node    
s = x(1:n); 
z = x(n+1 : end);
ds= (-beta)*( (1-u).* s .* (A * z) );
dz= -ds - gamma * z;
dx = [ds;dz];
end

%shorthand for getting the parameters from the main workspace
%fxt1 = @(x) fuxt(zeros(nodeNum,1),0,x,beta,gamma,A);