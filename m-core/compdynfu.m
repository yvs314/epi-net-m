%% Compute the Dynamics of Net-SIR with Fixed Control.
% Better IVs by running null-control pandemic

%Run the Net-SIR model with null control 
%to determine when does infection level reaches 1 in 10^5 everywhere
%and use it as initial conditions for the optimal control problem

%Author: Yaroslav Salii, 2022+
% Directories: 
    %./data/by-tract is for problem instances (IVs + travel matrix)
    %./fig is for the figures (not tracked by git)
    %./out is for tabular output (not tracked by git)    
% Naming conventions: fields separated by "-" subfields by "_" 
    % Instances:    (with string $NAME and integer $SIZE)
        % $NAME_$SIZE-init.csv holds initial conditions; 
        % $NAME_$SIZE-trav.dat holds daily travelers matrix (air passengers +
        % commuters)
        %sample: NWcty_75[-init.csv,-trav.dat], 75 counties in OR and WS
    % Output tables: as input, to /data/inst-100K
% 2022-02-23 v.0.0 chiseling away the unnecessary

%% TODO
%% Clear the workspace
clear; close all; %chuck all variables, close all figures etc.
%% Naming coventions setup
%./data/by-tract is for problem instances (IVs + travel matrix)
instDir="../data/by-tract";
ofigDir = "../fig"; %write the output figures and tables here
if(~exist(ofigDir,"dir"))
    mkdir(ofigDir); %make sure it exists
end
otabDir = "../data/inst-100K"; %write the output figures and tables here
if(~exist(otabDir,"dir"))
    mkdir(otabDir); %make sure it exists
end

fnSep="-"; %use - to separate file name fields
fnSubSep="_"; %use _ to subdivide file name fields


%YS: using regexp/extract as opposed to split by fnSep="-" because
%it's too complicated to extract the first element of split's return array
%in an anonymous function

%instance name as regexp pattern; 
rpInst = regexpPattern('[A-Z]+[0-9]*(tra|cty|ap|ste)_[0-9]+');

IV_suff="init.csv"; %all instances IVs end like this 
trav_suff="trav.dat"; %all instances' travel data end like this

%INPUT PATHS & FILENAMES
pathIV= @(iname) fullfile(instDir,iname+fnSep+IV_suff); %path to the IVs
pathTrav = @(iname) fullfile(instDir,iname+fnSep+trav_suff); %path to the travel data, if any

%OUTPUT PATHS & FILENAMES
pathotabs = @(iname) fullfile(otabDir,iname+"-abs.csv");
pathotabs0 = @(iname) fullfile(otabDir,iname+"-abs0.csv"); %for the NULL control
pathotfrac = @(iname) fullfile(otabDir,iname+"-frac.csv");
pathotfrac0 = @(iname) fullfile(otabDir,iname+"-frac0.csv"); %for the NULL control
pathotavg = @(iname) fullfile(otabDir,iname+"-avg.csv"); %all averaged
%LOG FILE, a CSV with a line for each iteration
pathotlog= @(iname) fullfile(otabDir,iname+"-log.csv"); 

%% Model Parameters
% these follow (El Ouardighi, Khmelnitsky, Sethi, 2020)

%INFECTION and RECOVERY RATES
%(Rader, Scarpino, Nande, &al., 2020): alpha = 0.15, gamma = 0.1
%baseline R_0=2.74; beta/gamma 63%
%set to harmonic mean of \alpha and \beta from idem
beta = 0.1196; % 1/(1/ 0.2977 + 1/0.2); %infection rate; (idem)
%1/beta = 8.36 mean time to be (symptomatic) infected
%beta = 1/2.5; %compat with older, R_0 = 3.3
gamma  = 0.0437; % 1/ (beta/R_0); R_0 = beta / gamma (idem)
%1/gamma = 22.904 mean time to recovery
%gamma = 1/8.3; %compat with older


%RUNNING COSTS      l for lockdown (control)
c = 0.00200; %running cost of infections; mean(c_1=100, c_2=300, c_3=200)
l = 0.000450; %running cost of control; (q3=450, securing social interactions)
%TERMINAL COST
k = 0.02000; %terminal cost of infections; mean(k_1=1K, k_2=3K, k_3=1K)
%FATIGUE RATES      (PI/PD-SF) 
r1 = 0.002; %infection rates fatigue rate
%r1 = 0; %infection rates fatigue rate
r2 = 0.002; %lockdown control fatigue rate
%r2 = 0;

%%TIME
T = 500; %time is [0,T], in days; default to 500 to hopefully reach endemic
%tArr = 0:T that's the de facto usage

%CONTROL bounds (for each component, at each time)
umin = 0; umax = 1; % u_i \in [0,1] \forall i \in nodes

%% Set Problem Instance Name

% read all instances IV files names' into array-of-structs
sAllIVs = dir(fullfile(instDir,'*-init.csv')); %array o' structs
cellAllIVs = extractfield(sAllIVs,'name'); %all IV file names (cell array)
cellAllInst = cellfun(@(s) extract(s,rpInst),cellAllIVs); %all instance names

%inst="NWtra_2072"; %by-tract OR + WS, with flights & commute
%inst="NWcty_75"; %by-county OR + WS, with flights & commute
%inst="NWap_23"; %by-airport OR + WS, with filghts & commute
inst="NWste_2"; %by-state OR + WS, with flights & commute

%inst="WCTtra_9110"; %WCT is CA + OR + WS
%inst="WCTcty_133";
%inst="WCTap_52";
%inst="WCTste_3";

%inst="ALLcty_3109";
%inst="ALLap_417";
%inst="ALLste_49";

%inst = "NW1tra_2834";
%inst = "NW1cty_136";
%inst = "NW1ap_37";
%inst = "NW1ste_4";

%inst = "NW2ste_8";
%inst = "NW2ap_74";
%inst = "NW2cty_259";
%inst = "NW2tra_4830";

%inst = "CAste_1"; %implementation doesn't support 1-node instances :(
%inst = "CAap_32";
%inst = "CActy_58";
%inst = "CAtra_7038";

% for inst = cellAllInst %it works!
%     disp(pathIV(inst));


%% Read Problem Instance (initial values, populations, and travel matrix)
% $IV_Path is a .CSV {id,AP_code,N_i,S_i,I_i,R_i,Name,LAT,LNG},
tIVs = readtable(pathIV(inst));

n = size(tIVs,1); %as many nodes as there are rows
N = table2array(tIVs(:,3)); %the population vector
iN = arrayfun(@(x) 1/x,N); %inverse pops, for Hadamard division by N (--> A)

NN = N / max(N); %normalized to max pop---for J, etc.
iNN = arrayfun(@(x) 1/x,NN); %inverse norm-pops, for Hadamard division (--> lax,u,J)
%NN = ones(n,1); iNN = ones(n,1); %FAT CAVEAT: "egalitarian" pricing

%STATE: initial conditions
s0 = table2array(tIVs(:,4)) .* iN; %susceptibles at t=0, frac
z0 = table2array(tIVs(:,5)) .* iN; %infecteds at t=0, frac
x0 = [s0;z0]; %overall state is [s;z]

%COSTATE: terminal values from transversality conditions (column vectors!)
%unused but I may find a use for these later
lasT = zeros(n,1); %\lambda_s(T) = 0_n
lazT = exp(r1*T)*k*NN; %\lambda_z(T) = e^{r_1T}k\tilde{N}
laxT = [lasT;lazT]; %costate is [\lambda_s; \lambda_z] 

% $pathTrav is just a matrix (floating point vals)
Araw = load(pathTrav(inst)); 
A = Araw'; %easier to add one transpose here than in every equation
% set its diagonal to pops N, then divide row-wise by N (traveling fracs)
A = diag(iN)* (A - diag(A) + diag(N));
%A = eye(n); %dumb debug: isolated nodes

%% Run Initialization

%STATE: 0-init [2n \times |tArr|] + initial conditions
x = zeros(n*2,T+1); x(:,1) = x0;
%COSTATE: 0-init [2n \times |tArr|] + terminal conditions
lax = zeros(n*2,T+1); 
lax(:,end) = laxT; 

%CONTROL
u = zeros(n,T+1); %constant zero control
ppu = pchip(0:T,u); %fit with monotone Fritsch--Carlson splines

%MISC
% delta = 1E-3; %min. relative error for norms of u,x,lax in stopping conditions
% ct = 0; %set the loop counter
% opts = odeset('InitialStep',1); %ensure the 1st step is at most 1 day long

%Run Flags
flg.computeCostate = 0;
flg.computeObjval = 0;
%% Sweep Log Setup
% a `table`, to be written out as csv
tablogNames = ["Iter","tcm","tdt",...
    "J","cZR","RR","ZZ"];

tablogTypes = ["uint32","duration","duration",...
    "double","double","double","double"];
    
%consider tcm::double to prevent truncation of ms

tablogSize=[1 size(tablogNames,2)];
tablog = table('Size',tablogSize,'VariableType',tablogTypes,'VariableNames',tablogNames);
tablog.tcm(1) = seconds(0);
%% Compute the Dynamics
tNumdynStart = tic;

    tIterStart = tic; %start the iteration stopwatch
    
    %store the previous loop's stuff
    oldu = u; ppu = pchip(0:T,u); %fit with monotone Fritsch--Carlson splines
    oldx = x; oldlax = lax; 
    
    %STATE
    ftx1 = @(t,x) futxp(ppval(ppu,t),t,x,beta,gamma,A);
   % disp('Run ode45 on IVP for state x---forwards from 0 to T')
    x_sln = ode45(ftx1, [0 T], x(:,1),opts);
    x = deval(x_sln, 0:T );

    %COSTATE
    if (flg.computeCostate)
        gtx1 = @(t,lax) guxtlp(ppval(ppu,t), deval(x_sln, t) ...
            , t, lax, beta, gamma, A, r1, c, NN);
    %   disp('Run ode45 on IVP for costate lax---backwards from T to 0');
        lax_sln = ode45(gtx1, [T 0], lax(:,end), opts);
        lax = deval(lax_sln, 0:T);
    end
    
    %OBJECTIVE FUNCTION
    if (flg.computeObjval)
        Lt1 = @(tArr) Ltxu(ppval(ppu,tArr),deval(x_sln,tArr),tArr,r1,r2,c,l,NN); %running cost
     %   disp('Compute the objective function J(u,x,T)');
        PsiT1 = PsiT(x(:,end),T,r1,NN,k); %terminal cost
        J = PsiT1 + integral(Lt1,0,T); %the objective function
    end
    %STATISTICS
    
    %recovered at day T
    RR_T = sum(N - (x(1:n,end) + x(n+1:end,end) ) .* N); 
    %infected at day T + recovered at day T: eff. all but susceptible at T
    cZR_T = sum( (ones(n,1) - x(1:n,end)) .* N);
    % infected at day T
    ZZ_T = sum(x(n+1:end,end) .* N);
    
    %cumulative infected, including recovered, daily
    %repmat(N, [1,numel([0:T])] ) - x(1:n,:).* N
    
    %ITERATION DONE
    %record the iteration's duration
    tIterEnd = toc(tIterStart); tdtsec = seconds(tIterEnd);
    %fill out the log, to be added as a row 
    %tablogNames = ["Iter","tcm","tdt","J","cZR","RR","ZZ"];
    currLog = {ct,tablog.tcm(1)+tdtsec,tdtsec,...
    J,round(cZR_T,4),round(RR_T,4),round(ZZ_T,4)};
    %KEEP NULL-CONTROL solution
    tablog(1,:) = currLog; 
    xNull = x; laxNull = lax; JNull = J; ZZNull = ZZ_T; cZRNull = cZR_T;
    sNull = xNull(1:n,:); zNull = xNull(n+1:end,:); rNull = (1 - sNull - zNull);
    
    disp(tablog(1,:)); %show the current log
    
tNumdynEnd = toc(tNumdynStart);

%% Evaluating the numerical solution

%slice the state into (s,z,r) compartments
s = x(1:n,:); z = x(n+1:end,:); r = (1 - s - z); %[s z r] for output

%normalize to 1-node model: sum all absolutes, and divide by pop
A1 = @(c) sum( c .* N ); 
a1 = @(c) A1(c) / sum(N);
%make the absolutes
%S = s .* N; Z = z.*N; R = r.*N; %multiply to get the absolute state
%S01 = sum(sNull .* N); Z01 = sum(zNull .* N); R0 = sum(rNull .* N);
s11 = a1(s); z11 = a1(z); r11 = a1(r); 
s01 = a1(sNull); z01 = a1(zNull); r01 = a1(rNull);
avgOut = [z01; z11; sum(u) / n]'; %[total zNull; total z; avg u]


%% Tabular output
%column names for frac and abs tables (will be set to uppercase for abs)
cns = [arrayfun( @(n) 's'+string(n),0:T) ...
     arrayfun( @(n) 'z'+string(n),0:T) ...
     arrayfun( @(n) 'r'+string(n),0:T)];

%cns_u = [ arrayfun( @(n) 'u'+string(n),0:T) ]; %col names for controls
%cns_frac_u = horzcat(cns,cns_u); %combined frac with control

%col. names for per-node average infected-null, infected-opt, control effort
cns_avgc = ["zNull_avg","z_avg","u_avg"];
 
%construct the solution output tables
otabs = horzcat(tIVs,array2table([s z r] .* N,'VariableNames',upper(cns)));

%construct the average output table
otabavgc = array2table(avgOut, 'VariableNames',cns_avgc);
%note that the log table is already constructed

%write the solution output tables
writetable(otabs,pathotabs(inst));
%write the FBsweep iteration log table
writetable(tablog,pathotlog(inst));

%end %end this dumb instance name loop

%% AUXILIARY FUNCTIONS
heatmaplog=@(x) heatmap(x,'GridVisible','off','Colormap',flip(autumn),'ColorScaling','log');
heatmap1=@(x) heatmap(x,'GridVisible','off','Colormap',flip(autumn));
heatmap3 = @(u) heatmap(u,'GridVisible','off','Colormap',cool);

%Hi, I'm the terminal cost in the objective J and I'm too small for a separate file
function PsiT = PsiT(xT,T,r1,NN,k)
n = size(xT,1) / 2; zT = xT(n+1:end,:);
PsiT = exp(r1*T) * k * dot(zT,NN);
end




