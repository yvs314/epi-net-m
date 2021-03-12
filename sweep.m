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
    % Output tables: $NAME_$SIZE-$CPG-[abs | abs0 | frac | frac0].csv
% 2021-03-04 v.0.1 up to reading IVs and travel matrix
% 2021-03-08 v.0.2 numeric solutions for state & co-state eqs
% 2021-03-08 v.0.3 control computed from state & costate values
% 2021-03-09 v.0.4 control spline-fitted, monotone piecewise cubic H.(pchip)
% 2021-03-10 v.0.5 done forward-backward sweep
% 2021-03-12 v.0.6 done obj function and results export

%% TODO
% 1 debug output (time, errors, J, &c) to a log file
% 2 switch to tArr instead of “magical numbers” 0:T / T:0 / [0 T]
%% Clear the workspace
clear; close all; %chuck all variables, close all figures etc.
%% Naming coventions setup
%./data/by-tract is for problem instances (IVs + travel matrix)
instDir="data/by-tract";
ofigDir = "fig"; %write the output figures and tables here
if(~exist(ofigDir,"dir"))
    mkdir(ofigDir); %make sure it exists
end
otabDir = "out"; %write the output figures and tables here
if(~exist(otabDir,"dir"))
    mkdir(otabDir); %make sure it exists
end

fnSep="-"; %use - to separate file name fields
fnSubSep="_"; %use _ to subdivide file name fields
IV_suff="init.csv"; %all instances IVs end like this 
trav_suff="trav.dat"; %all instances' travel data end like this

%INPUT PATHS & FILENAMES
pathIV= @(iname) fullfile(instDir,iname+fnSep+IV_suff); %path to the IVs
pathTrav = @(iname) fullfile(instDir,iname+fnSep+trav_suff); %path to the travel data, if any

%OUTPUT PATHS & FILENAMES
pathotabs = @(iname) fullfile(otabDir,iname+"-abs.csv");
pathotabs0 = @(iname) fullfile(otabDir,iname+"-abs0.csv");
pathotfrac = @(iname) fullfile(otabDir,iname+"-frac.csv");
pathotfrac0 = @(iname) fullfile(otabDir,iname+"-frac0.csv"); %for the NULL control


%% Model Parameters
% these follow (El Ouardighi, Khmelnitsky, Sethi, 2020)

%INFECTION and RECOVERY RATES
%baseline R_0=2.74; beta/gamma 63%
%set to harmonic mean of \alpha and \beta from idem
beta = 0.1196; % 1/(1/ 0.2977 + 1/0.2); %infection rate; (idem)
%1/beta = 8.36 mean time to be (symptomatic) infected
%beta = 1/2.5; %compat with older
gamma  = 0.0437; % 1/ (beta/R_0); R_0 = beta / gamma (idem)
%1/gamma = 22.904 mean time to recovery
%gamma = 1/8.3; %compat with older


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
%tArr = 0:T that's the de facto usage

%CONTROL bounds (for each component, at each time)
umin = 0; umax = 1; % u_i \in [0,1] \forall i \in nodes

%% Problem Instance (initial values, populations, and travel matrix)
%inst="a~NW~tra_2072"; %by-tract OR + WS, with flights & commute
%inst="a~NW~cty_75"; %by-county OR + WS, with flights & commute
inst="a~NW~ste_2"; %by-state OR + WS, with flights & commute

% $IV_Path is a .CSV {id,AP_code,N_i,S_i,I_i,R_i,Name,LAT,LNG},
tIVs = readtable(pathIV(inst));

n = size(tIVs,1); %as many nodes as there are rows
N = table2array(tIVs(:,3)); %the population vector
iN = arrayfun(@(x) 1/x,N); %inverse pops, for Hadamard division by N etc.

%STATE: initial conditions
s0 = table2array(tIVs(:,4)) .* iN; %susceptibles at t=0, frac
z0 = table2array(tIVs(:,5)) .* iN; %infecteds at t=0, frac
x0 = [s0;z0]; %overall state is [s;z]

%COSTATE: terminal values from transversality conditions (column vectors!)
lasT = zeros(n,1); %\lambda_s(T) = 0_n
lazT = exp(r1*T)*k*N; %\lambda_z(T) = e^{r_1T}kN
laxT = [lasT;lazT]; %costate is [\lambda_s; \lambda_z]

% $pathTrav is just a matrix (floating point vals)
Araw = load(pathTrav(inst)); 
% set its diagonal to pops N, then divide row-wise by N (traveling fracs)
A = diag(iN)* (Araw - diag(Araw) + diag(N));


%% Sweep setup

%STATE: 0-init [2n \times |tArr|] + initial conditions
x = zeros(n*2,T+1); x(:,1) = x0;
%COSTATE: 0-init [2n \times |tArr|] + terminal conditions
lax = zeros(n*2,T+1); lax(:,end) = laxT;

%CONTROL
u = zeros(n,T+1); %initial guess: constant zero 
ppu = pchip(0:T,u); %fit with monotone Fritsch--Carlson splines

%MISC
delta = 0.001; %min. relative error for norms of u,x,lax in stopping conditions
test = -1.0; %ensure the loop is entered
stop_u = false; stop_x = false; stop_lax = false;%ensure the loop is entered
ct = 0; %set the loop counter


%% Forward-Backward Sweep Loop
while( ~stop_u || ~stop_x || ~stop_lax ) %while at least one rerr is > delta
    fprintf('Loop no. %d\n',ct); ct = ct+1;
    
    oldu = u; ppu = pchip(0:T,u); %fit with monotone Fritsch--Carlson splines
    oldx = x;
    oldlax = lax;
    
    %STATE
    ftx1 = @(t,x) futxp(ppval(ppu,t),t,x,beta,gamma,A);
    disp('Run ode45 on IVP for state x---forwards from 0 to T')
    tic
        x_sln = ode45(ftx1, [0 T], x(:,1));
        x = deval(x_sln, 0:T );
    toc 

    %COSTATE
    gtx1 = @(t,lax) guxtlp(ppval(ppu,t), deval(x_sln, t) ...
        , t, lax, beta, gamma, A, r1, c, N);
    disp('Run ode45 on IVP for costate lax---backwards from T to 0');
    tic
        lax_sln = ode45(gtx1, [T 0], lax(:,end));
        lax = deval(lax_sln, 0:T);
    toc

    %OBJECTIVE FUNCTION
    Lt1 = @(tArr) Ltxu(ppval(ppu,tArr),deval(x_sln,tArr),tArr,r1,r2,c,l,N); %running cost
    disp('Compute the objective function J(u,x,T)');
    tic 
        PsiT1 = PsiT(x(:,end),T,r1,N,k); %terminal cost
        J = PsiT1 + integral(Lt1,0,T); %the objective function
    toc 
    
    %CONTROL
    utxla1 = @(tArr,xtArr,laxtArr) utxla(tArr,xtArr,laxtArr,umin,umax,beta,l,A,r2,iN);
    disp('Compute the optimal control at points 0..T');
    tic; u1 = utxla1(0:T,x,lax); toc
    u = 0.5*(u1 + oldu); %gentle update of u (convex combination)
    
    fprintf('\nJ = %E\n',J); %print the objective function
    %STOPPING CONDITIONS (rel. err. \delta||_|| - ||old_ - _|| > 0)
    rerr_u = delta*norm(u,1) - norm(oldu - u,1); stop_u = rerr_u > 0;
    rerr_x = delta*norm(x,1) - norm(oldx - x,1); stop_x = rerr_x > 0;
    rerr_lax = delta*norm(lax,1) - norm(oldlax - lax,1); stop_lax = rerr_lax > 0;
    fprintf('rerr_u = %4.4f   rerr_x = %4.4f   rerr_lax = %4.4f\n',rerr_u,rerr_x,rerr_lax);
    
    %HUMAN-READABLE STOPPING CONDITIONS (relative error ||old_ - _|| / ||_|| < delta)
    hrerr_u = norm(oldu - u,1)/ norm(u,1);
    hrerr_x = norm(oldx - x,1) / norm(x,1);
    hrerr_lax = norm(oldlax - lax,1) / norm(lax,1); 
    fprintf('hrerr_u = %4.4f   hrerr_x = %4.4f   hrerr_lax = %4.4f\n\n' ... 
        ,hrerr_u,hrerr_x,hrerr_lax);
    
    %KEEP NULL-CONTROL solution
    if(ct == 1)
        xNull = x; laxNull = lax;
    end
end %next sweep iteration

%% Tabular output
%slice the state into (s,z,r) compartments
s = x(1:n,:); z = x(n+1:end,:); r = (1 - s - z); %[s z r] for output
sNull = xNull(1:n,:); zNull = xNull(n+1:end,:); rNull = (1 - sNull - zNull);

cns = [arrayfun( @(n) 's'+string(n),0:T) ...
     arrayfun( @(n) 'z'+string(n),0:T) ...
     arrayfun( @(n) 'r'+string(n),0:T)];
 
otfrac = horzcat(tIVs,array2table([s z r],'VariableNames',cns));
otfrac0 = horzcat(tIVs,array2table([sNull zNull rNull],'VariableNames',cns));
otabs = horzcat(tIVs,array2table([s z r] .* N,'VariableNames',upper(cns)));
otabs0 = horzcat(tIVs,array2table([sNull zNull rNull] .* N,'VariableNames',upper(cns)));

writetable(otfrac,pathotfrac(inst));
writetable(otfrac0,pathotfrac0(inst));
writetable(otabs,pathotabs(inst));
writetable(otabs0,pathotabs0(inst));

%% AUXILIARY FUNCTIONS

%Hi, I'm the terminal cost in the objective J and I'm too small for a separate file
function PsiT = PsiT(xT,T,r1,N,k)
n = size(xT,1) / 2; zT = xT(n+1:end,:);
PsiT = exp(r1*T) * k * dot(zT,N);
end

%% BIT BUCKET
