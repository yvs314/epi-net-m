%% Forward-Backward Sweep for Net-SIR with Social Distancing Control

%Author: Yaroslav Salii, 2021+
% Directories: 
    %./data/by-tract is for problem instances (IVs + travel matrix)
    %./fig is for the figures (not tracked by git)
    %./out is for tabular output (not tracked by git)    
% Naming conventions: fields separated by "-" subfields by "_" 
    % Instances:    (with string $NAME and integer $SIZE)
        % $NAME_$SIZE-init.csv holds initial conditions; 
        % $NAME_$SIZE-trav.dat holds daily travelers matrix (air passengers +
        % commuters)
        %sample: a~NW~cty_75[-init.csv,-trav.dat], 75 counties in OR and WS
    % Output tables: $NAME_$SIZE-$CPG-[abs | frac]-out.csv
% 2021-03-04 v.0.1 up to reading IVs and travel matrix
% 2021-03-08 v.0.2 numeric solutions for state & co-state eqs
% 2021-03-08 v.0.3 control computed from state & costate values
% 2021-03-09 v.0.4 control spline-fitted, monotone piecewise cubic H.(pchip)

%% TODO
% switch to tArr instead of “magical numbers” 0:T / T:0 / [0 T]
%% Clear the workspace
clear; close all; %chuck all variables, close all figures etc.
%% Naming coventions setup
%./data/by-tract is for problem instances (IVs + travel matrix)
instDir="data/by-tract";
oDirFig = "fig"; %write the output figures and tables here
if(~exist(oDirFig,"dir"))
    mkdir(oDirFig); %make sure it exists
end
oDirTab = "out"; %write the output figures and tables here
if(~exist(oDirTab,"dir"))
    mkdir(oDirTab); %make sure it exists
end

fnSep="-"; %use - to separate file name fields
fnSubSep="_"; %use _ to subdivide file name fields
IV_suff="init.csv"; %all instances IVs end like this 
trav_suff="trav.dat"; %all instances' travel data end like this

%INPUT PATHS & FILENAMES
pathIV= @(iname) fullfile(instDir,iname+fnSep+IV_suff); %path to the IVs
pathTrav = @(iname) fullfile(instDir,iname+fnSep+trav_suff); %path to the travel data, if any

%OUTPUT PATHS & FILENAMES
pathOutTabAbs = @(iname) fullfile(outDir,join([iname,"abs"],fnSep));
pathOutTabFrac = @(iname) fullfile(outDir,join([iname,"frac"],fnSep));


%% Model Parameters
% these follow (El Ouardighi, Khmelnitsky, Sethi, 2020)

%INFECTION and RECOVERY RATES
%baseline R_0=2.74; beta/gamma 63%
%set to harmonic mean of \alpha and \beta from idem
beta = 0.1196; % 1/(1/ 0.2977 + 1/0.2); %infection rate; 
%1/beta = 8.36 mean time to be (symptomatic) infected
gamma  = 0.0437; % 1/ (beta/R_0); R_0 = beta / gamma 
%1/gamma = 22.904 mean time to recovery

%RUNNING COSTS      l for lockdown (control)
c = 200; %running cost of infections; mean(c_1=100, c_2=300, c_3=200)
l = 450; %running cost of control; (q3=450, securing social interactions)
%TERMINAL COST
k = 2000; %terminal cost of infections; mean(k_1=1K, k_2=3K, k_3=1K)
%FATIGUE RATES      (PI/PD-SF) 
r1 = 0.002; %infection rates fatigue rate
r2 = 0.002; %lockdown control fatigue rate

%%TIME
T = 180; %time is [0,T], in days
%tArr = 0:T that's the de fact usage

%CONTROL bounds (for each component, at each time)
umin = 0; umax = 1; % u_i \in [0,1] \forall i \in nodes

%% Problem Instance (initial values, populations, and travel matrix)
%inst="a~NW~tra_2072"; %by-tract OR + WS, with flights & commute
%inst="a~NW~cty_75"; %by-county OR + WS, with flights & commute
inst="a~NW~ste_2"; %by-state OR + WS, with flights & commute

% $IV_Path is a .CSV {id,AP_code,N_i,S_i,I_i,R_i,Name,LAT,LNG},
tIVs = readtable(pathIV(inst));

nodeNum = size(tIVs,1); %as many nodes as there are rows
N = table2array(tIVs(:,3)); %the population vector
iN = arrayfun(@(x) 1/x,N); %inverse pops, for Hadamard division by N etc.

%STATE: initial conditions
s0=table2array(tIVs(:,4)) .* iN; %susceptibles at t=0, frac
z0=table2array(tIVs(:,5)) .* iN; %infecteds at t=0, frac

%COSTATE: terminal values from transversality conditions (column vectors!)
lasT = zeros(nodeNum,1); %\lambda_s(T) = 0_n
lazT = exp(r1*T)*k*N; %\lambda_z(T) = e^{r_1T}kN

% $pathTrav is just a matrix (floating point vals)
Araw = load(pathTrav(inst)); 
% set its diagonal to pops N, then divide row-wise by N (traveling fracs)
A = diag(iN)* (Araw - diag(Araw) + diag(N));


%% Sweep setup

%0-INIT: state, co-state, and control; caveat: wrong discretization
s = zeros(nodeNum,T+1);
z = zeros(nodeNum,T+1);

las = zeros(nodeNum,T+1);
laz = zeros(nodeNum,T+1);

%CONTROL
u = zeros(nodeNum,T+1); %initial guess: constant zero 
ppu = pchip(0:T,u); %fit with monotone Fritsch--Carlson splines

s(:,1) = s0; z(:,1) = z0; %set initial conditions for state
las(:,T+1)=lasT; laz(:,T+1)=lazT; %set terminal conditions for costate

%combined state & costate vectors, [2n \times tSpan]
x = [s; z];
lax = [las; laz];
%% State equation


%compute state eq. dx = [ds; dz], compatible with ode45, spline-fitted control
ftx1 = @(t,x) futxp(ppval(ppu,t),t,x,beta,gamma,A);

%SOLVE with zero control (the zeros(n,1) in fxt2)
disp('Run ode45 on IVP for state x---forwards from 0 to T')
tic
    x_sln = ode45(ftx1, [0 T], x(:,1)); %consider setting 'NonNegative' flag
    x = deval(x_sln, 0:T );
toc 

%set the parameters for costate's RHS column, compat with ode45
%costate is supplied through deval(x_sln, t)
gtx1 = @(t,lax) guxtlp(ppval(ppu,t), deval(x_sln, t) ...
    , t, lax ...
    , beta, gamma, A, r1, c, N);

                    
disp('Run ode45 on IVP for costate lax---backwards from T to 0');
tic
    lax_sln = ode45(gtx1, [T 0], lax(:,end));
    lax = deval(lax_sln, 0:T);
toc

%% Compute optimal control---test

%set the parameters for the vector control computation
utxla1 = @(tArr,xtArr,laxtArr) utxla(tArr,xtArr,laxtArr,umin,umax,beta,l,A,r2,iN);
disp('Compute the optimal control at points 0..T');
tic
    unew = utxla1(0:T,x,lax);
toc

