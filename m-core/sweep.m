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
        %sample: NWcty_75[-init.csv,-trav.dat], 75 counties in OR and WS
    % Output tables: $NAME_$SIZE-[abs | abs0 | frac | frac0 | avg | log].csv
% 2021-03-04 v.0.1 up to reading IVs and travel matrix
% 2021-03-08 v.0.2 numeric solutions for state & co-state eqs
% 2021-03-08 v.0.3 control computed from state & costate values
% 2021-03-09 v.0.4 control spline-fitted, monotone piecewise cubic H.(pchip)
% 2021-03-10 v.0.5 done forward-backward sweep
% 2021-03-12 v.0.6 done obj function and results export
% 2021-03-15 v.0.6.1 normalized pop sizes in objective function
%      "     v.0.6.2 flags to (a) bound lambda (b) set term. cost to 0
% 2021-03-16 v.1.0 done adding control update strategies
%      "     v.1.1 fix the criterion for zero controls
% 2021-03-18 v.1.1.1 add the transpose, made control 10 times cheaper
% 2021-03-25 v.1.2 add ZZ: abs. infected at T and cZR: abs. Z + R (T)
%            v.1.2.1 add norm-to-1-node output avgOut [z01; z11; avg_u]
% 2021-10-04 v.1.2.2 rehaul avgOut to output as a CSV with header
% 2021-10-24 v.1.2.3 add control effort output to _-frac.csv
% 2022-02-17 v.1.3   add iteration logging via _-log.csv (rough-in)
% 2022-03-02 v.1.4   add multy-instance logging via "sweep-summary.csv"

%% TODO
% 0 combine all output into a single file
% 1 debug output (time, errors, J, &c) to a log file
% 2 switch to tArr instead of “magical numbers” 0:T / T:0 / [0 T]
% 3 set asserts(condition,"error message") instead of comment caveats
%% Clear the workspace
clear; close all; %chuck all variables, close all figures etc.
%% Naming coventions setup
%./data/by-tract is default for problem instances (IVs + travel matrix)
instDir="../data/by-tract";%the default
%instDir="../data/inst-100K";
travDir="../data/by-tract";
ofigDir = "../fig"; %write the output figures and tables here
if(~exist(ofigDir,"dir"))
    mkdir(ofigDir); %make sure it exists
end
otabDir = "../out"; %write the output tables here
%otabDir = "../out-pz-rerun"; %write the output tables here
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
pathTrav = @(iname) fullfile(travDir,iname+fnSep+trav_suff); %path to the travel data, if any

%OUTPUT PATHS & FILENAMES
pathotabs = @(iname) fullfile(otabDir,iname+"-abs.csv");
pathotabs0 = @(iname) fullfile(otabDir,iname+"-abs0.csv"); %for the NULL control
pathotfrac = @(iname) fullfile(otabDir,iname+"-frac.csv");
pathotfrac0 = @(iname) fullfile(otabDir,iname+"-frac0.csv"); %for the NULL control
pathotavg = @(iname) fullfile(otabDir,iname+"-avg.csv"); %all averaged
%LOG FILE, a CSV with a line for each iteration
pathotlog= @(iname) fullfile(otabDir,iname+"-log.csv"); 
%SUMMARY, with all 
pathotsummary = fullfile(otabDir,"sweep-summary.csv");

tabsummary = table; %empty table, one row per instance

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
T = 180; %time is [0,T], in days; default to 180
%tArr = 0:T that's the de facto usage

%CONTROL bounds (for each component, at each time)
umin = 0; umax = 1; % u_i \in [0,1] \forall i \in nodes

%TEST FLAGS
killPsi = false; %(if true, set the terminal cost to 0)
laxmax=1000; laxmin = -laxmax;
boundlax = false;
%boundlax = true; %enforce lax is within bounds
%% Set Problem Instance Name

% read all instances IV files names' into array-of-structs
sAllIVs = dir(fullfile(instDir,'*-init.csv')); %array o' structs
cellAllIVs = extractfield(sAllIVs,'name'); %all IV file names (cell array)
cellAllInst = cellfun(@(s) extract(s,rpInst),cellAllIVs); %all instance names

%inst="NWtra_2072"; %by-tract OR + WS, with flights & commute
%inst="NWcty_75"; %by-county OR + WS, with flights & commute
%inst="NWap_23"; %by-airport OR + WS, with filghts & commute
%inst="NWste_2"; %by-state OR + WS, with flights & commute

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

%for inst = cellAllInst(16:16) %default to no. 16, "NWste_2" 
for inst = cellAllInst %run all instances found in instDir
    disp(pathIV(inst));


%% Read Problem Instance (initial values, populations, and travel matrix)
% $IV_Path is a .CSV {id,AP_code,N_i,S_i,I_i,R_i,Name,LAT,LNG},
tabIVs = readtable(pathIV(inst));

n = size(tabIVs,1); %as many nodes as there are rows
N = table2array(tabIVs(:,3)); %the population vector
iN = arrayfun(@(x) 1/x,N); %inverse pops, for Hadamard division by N (--> A)

NN = N / max(N); %normalized to max pop---for J, etc.
iNN = arrayfun(@(x) 1/x,NN); %inverse norm-pops, for Hadamard division (--> lax,u,J)
%NN = ones(n,1); iNN = ones(n,1); %FAT CAVEAT: "egalitarian" pricing

%STATE: initial conditions
s0 = table2array(tabIVs(:,4)) .* iN; %susceptibles at t=0, frac
z0 = table2array(tabIVs(:,5)) .* iN; %infecteds at t=0, frac
x0 = [s0;z0]; %overall state is [s;z]

%COSTATE: terminal values from transversality conditions (column vectors!)
lasT = zeros(n,1); %\lambda_s(T) = 0_n
lazT = exp(r1*T)*k*NN; %\lambda_z(T) = e^{r_1T}k\tilde{N}
laxT = [lasT;lazT]; %costate is [\lambda_s; \lambda_z] 

% $pathTrav is just a matrix (floating point vals)
Araw = load(pathTrav(inst)); 
A = Araw'; %easier to add one transpose here than in every equation
% set its diagonal to pops N, then divide row-wise by N (traveling fracs)
A = diag(iN)* (A - diag(A) + diag(N));
%A = eye(n); %dumb debug: isolated nodes

%% Sweep setup

%STATE: 0-init [2n \times |tArr|] + initial conditions
x = zeros(n*2,T+1); x(:,1) = x0;
%COSTATE: 0-init [2n \times |tArr|] + terminal conditions
lax = zeros(n*2,T+1); 
if(~killPsi) 
    lax(:,end) = laxT; end

%CONTROL
u = zeros(n,T+1); %initial guess: constant zero 
%u = ones(n,T+1); %other knee-jerk initial guess: constant 1
ppu = pchip(0:T,u); %fit with monotone Fritsch--Carlson splines

%MISC
delta = 1E-3; %min. relative error for norms of u,x,lax in stopping conditions
stop_u = false; stop_x = false; stop_lax = false;%ensure the loop is entered
ct = 0; %set the loop counter
opts = odeset('InitialStep',1,'AbsTol',1e-9,'RelTol',1e-6); %ensure the 1st step is at most 1 day long

%STATUS FLAGS
flg.ct_within_limit = true;
flg.objective_improved = true;
%flg.within_time_limit = true;
%% Sweep Log Setup
% a `table`, to be written out as csv
tablogNames = ["Iter","tcm","tdt",...
    "J","cZR","RR","ZZ",...
    "hrerr_u","hrerr_x","hrerr_lax",...
    "rerr_u","rerr_x","rerr_lax"];

tablogTypes = ["uint32","duration","duration",...
    "double","double","double","double",...
    "double","double","double",...
    "double","double","double"];
%consider tcm::double to prevent truncation of ms

tablogSize=[1 size(tablogNames,2)];
tablog = table('Size',tablogSize,'VariableType',tablogTypes,'VariableNames',tablogNames);
tablog.tcm(1) = seconds(0);
%% Forward-Backward Sweep Loop
tSweepStart = tic;
while( ~stop_u || ~stop_x || ~stop_lax ) %while at least one rerr is > delta
   % fprintf('Loop no. %d\n',ct); 
    ct = ct+1; %%inc the iteration counter
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
    gtx1 = @(t,lax) guxtlp(ppval(ppu,t), deval(x_sln, t) ...
        , t, lax, beta, gamma, A, r1, c, NN);
%   disp('Run ode45 on IVP for costate lax---backwards from T to 0');
    lax_sln = ode45(gtx1, [T 0], lax(:,end), opts);
    lax = deval(lax_sln, 0:T);
    if(boundlax)
        lax = min(laxmax,max(lax,laxmin));
    %    pplax = spline(0:T,lax);
    end


    %OBJECTIVE FUNCTION
    Lt1 = @(tArr) Ltxu(ppval(ppu,tArr),deval(x_sln,tArr),tArr,r1,r2,c,l,NN); %running cost
 %   disp('Compute the objective function J(u,x,T)');
    if (~killPsi)
        PsiT1 = PsiT(x(:,end),T,r1,NN,k); %terminal cost
    else
        PsiT1 = 0; %no terminal cost, no transversality
    end
    J = PsiT1 + integral(Lt1,0,T); %the objective function

    %STATISTICS
    
    %recovered at day T
    RR = sum(N - (x(1:n,end) + x(n+1:end,end) ) .* N); 
    %infected at day T + recovered at day T: eff. all but susceptible at T
    cZR = sum( (ones(n,1) - x(1:n,end)) .* N);
    % infected at day T
    ZZ = sum(x(n+1:end,end) .* N);
    
    %cumulative infected, including recovered, daily
    %repmat(N, [1,numel([0:T])] ) - x(1:n,:).* N
    
    %CONTROL
    utxla1 = @(tArr,xtArr,laxtArr) utxla(tArr,xtArr,laxtArr,umin,umax,beta,l,A,r2,iNN);
    u1 = utxla1(0:T,x,lax); 
    
    %UPDATE THE CONTROL
    %pick the update that better improves J??? 'd require to keep 2 last
    %u = 0.5*(u1 + oldu); %gentle update of u (convex combination)
    u = 0.9*oldu + 0.1*u1; %a more gentle update of u (convex combination)
    %u = u1; %just forget the old stuff: feels bad
    %bbupd1 = @(ct,u1,oldu) bbupd(0.99,umin,umax,ct,u1,oldu);
    %u = bbupd1(ct,u1,oldu);
    
    
    %STOPPING CONDITIONS (rel. err. \delta||_|| - ||old_ - _|| > 0)
    rerr_u = delta*norm(u,1) - norm(oldu - u,1); stop_u = rerr_u >= 0;
    rerr_x = delta*norm(x,1) - norm(oldx - x,1); stop_x = rerr_x >= 0;
    rerr_lax = delta*norm(lax,1) - norm(oldlax - lax,1); stop_lax = rerr_lax >= 0;
    
    %HUMAN-READABLE STOPPING CONDITIONS (relative error ||old_ - _|| / ||_|| < delta)
    hrerr_u = norm(oldu - u,1)/ norm(u,1);
    hrerr_x = norm(oldx - x,1) / norm(x,1);
    hrerr_lax = norm(oldlax - lax,1) / norm(lax,1); 
    
    
    %ITERATION DONE
    %record the iteration's duration
    tIterEnd = toc(tIterStart); tdtsec = seconds(tIterEnd);
    %fill out the log, to be added as a row 
    %tablogNames = ["Iter","tcm","tdt","J","cZR","RR","ZZ","hrerr_u","hrerr_x","hrerr_lax","rerr_u","rerr_x","rerr_lax"];
    currLog = {ct,tablog.tcm(1)+tdtsec,tdtsec,...
        J,round(cZR,4),round(RR,4),round(ZZ,4),hrerr_u,hrerr_x,hrerr_lax,rerr_u,rerr_x,rerr_lax};
    %KEEP NULL-CONTROL solution
    if(ct == 1) %first iteration-special
        tablog(1,:) = currLog; 
        xNull = x; laxNull = lax; JNull = J; ZZNull = ZZ; cZRNull = cZR;
        sNull = xNull(1:n,:); zNull = xNull(n+1:end,:); rNull = (1 - sNull - zNull);
    else
        tablog = [currLog; tablog]; %new iter. log on top
        %I'd preallocate if I knew the number of iterations in advance
    end
    
    disp(tablog(1,:)); %show the current log
    
    if(ct > 500)
        %error("That probably wouldn't converge. Terminating.");
        disp("That probably wouldn't converge. Stopping the sweep.");
        flg.ct_within_limit = false;
        break;
    end
end %next sweep iteration
tSweepEnd = toc(tSweepStart);

%% Evaluating the sweep results
fprintf("Sweep done in "); disp(seconds(tSweepEnd));
% fprintf('\n J / JNull = %4.4f   improved by %4.4f\n',J / JNull, 1 - J/JNull);
% fprintf('ZZNull = %d   ZZ = %d  cZRNull = %d cZR = %d\n', ... 
%     ceil(ZZNull), ceil(ZZ), ceil(cZRNull), ceil(cZR) );

if(J / JNull > 1)
   % error("Didn't improve over initial guess. Terminating.");
   disp("Didn't improve over initial guess. Setting failure flag");
   flg.objective_improved = false;
end

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

rowsummary = table;
rowsummary.Inst = inst;
rowsummary.Time = sum(tablog.tdt);
rowsummary.JNull = round(tablog.J(end),4);
rowsummary.JOpt = round(tablog.J(1),4);
rowsummary.JImprPercent = round( 1 - rowsummary.JOpt / rowsummary.JNull ,4)*100;
rowsummary.Iter = tablog.Iter(1);
rowsummary.cZRNull = ceil(tablog.cZR(end));
rowsummary.cZROpt = ceil(tablog.cZR(1));
rowsummary.zPeakOpt = round(max(z11),4);
rowsummary.zPeakNull = round(max(z01),4);
rowsummary.success = flg.ct_within_limit && flg.objective_improved;

disp(rowsummary);

%add to the summary if it wasn't empty
if(isempty(tabsummary))
    tabsummary = rowsummary;
else
    tabsummary = [rowsummary; tabsummary];
end

%% Tabular output
%column names for frac and abs tables (will be set to uppercase for abs)
cns = [arrayfun( @(n) 's'+string(n),0:T) ...
     arrayfun( @(n) 'z'+string(n),0:T) ...
     arrayfun( @(n) 'r'+string(n),0:T)];

cns_u = [ arrayfun( @(n) 'u'+string(n),0:T) ]; %col names for controls
cns_frac_u = horzcat(cns,cns_u); %combined frac with control

%col. names for per-node average infected-null, infected-opt, control effort
cns_avgc = ["zNull_avg","z_avg","u_avg"];
 
%construct the solution output tables
otfrac = horzcat(tabIVs,array2table([s z r u],'VariableNames', cns_frac_u) ); %adding the controls here
otfrac0 = horzcat(tabIVs,array2table([sNull zNull rNull],'VariableNames',cns));
otabs = horzcat(tabIVs,array2table([s z r] .* N,'VariableNames',upper(cns)));
otabs0 = horzcat(tabIVs,array2table([sNull zNull rNull] .* N,'VariableNames',upper(cns)));
%construct the average output table
otabavgc = array2table(avgOut, 'VariableNames',cns_avgc);
%note that the log table is already constructed

%write the solution output tables
writetable(otfrac,pathotfrac(inst));
writetable(otfrac0,pathotfrac0(inst));
writetable(otabs,pathotabs(inst));
writetable(otabs0,pathotabs0(inst));
%write the average output table
writetable(otabavgc,pathotavg(inst));
%write the FBsweep iteration log table
writetable(tablog,pathotlog(inst));

end %end this dumb instance name loop

% write the summary
writetable(tabsummary,pathotsummary);

%% AUXILIARY FUNCTIONS
heatmaplog=@(x) heatmap(x,'GridVisible','off','Colormap',flip(autumn),'ColorScaling','log');
heatmap1=@(x) heatmap(x,'GridVisible','off','Colormap',flip(autumn));
heatmap3 = @(u) heatmap(u,'GridVisible','off','Colormap',cool);

%Hi, I'm the terminal cost in the objective J and I'm too small for a separate file
function PsiT = PsiT(xT,T,r1,NN,k)
n = size(xT,1) / 2; zT = xT(n+1:end,:);
PsiT = exp(r1*T) * k * dot(zT,NN);
end

%update control with exponential backoff to umin and additive to umax
function u = bbupd(a,umin,umax,ct,u1,oldu)
u = zeros(size(u1)); %preallocate the control
b = a^ct; %precompute the backoff value
for t = 1:size(u,2)
    for nd = 1:size(u,1)
        if(u1(nd,t) > oldu(nd,t)) %increased since last loop
           u(nd,t) = umax*(1-a) + a*oldu(nd,t); %additive to umax
           %u(nd,t) = umax*(1-b) + b*oldu(nd,t); %exponential to umax
           %u(nd,t) = umax*(b) + (1-b)*oldu(nd,t); 
           %u(nd,t) = oldu(nd,t);
           %u(nd,t) = 0.999*oldu(nd,t) + 0.001*u1(nd,t);
        else
            u(nd,t) = umin*(1-b) + b*oldu(nd,t); %exponential to umin
        end
    end
end
end
%would look swell in zip-like statement. ohwait, it's Matl~



%% BIT BUCKET
