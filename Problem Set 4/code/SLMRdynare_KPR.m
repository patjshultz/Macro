%% ================================================================= 
% SLMRdynare.m
%
% matlab script that generates dynare output
%===================================================================
close all
clear
format shortg

%% =================================================================
% Preliminaries 
%===================================================================

% Utility function is U(C_t, 1- h_t)= [C_t^theta](1-h_t)^(1-theta)^(1-eta)/(1-eta)
% Variables are in logs
%  (the steady state is first computed in levels)
% ================================================================
% Set parameters
%==================================================================

delta=  0.025;                % depreciation rate of capital
gamma = 1.008;                % steady-state growth rate
alpha = 2/3;                 % Labor share of income

bbeta=  0.99;                % household discout factor

rho_a =  0.6;                % persistence of productivity shock
sigma_a =  0.0066;             % volatility of productivity shock. 
%abar= 0;                     % log mean of productivity shock

%{
rho_g =  0.95;                % persistence of government shock
sigma_g =  0.007;             % volatility of government shock. 
gy= 1/6;                      % mean government to GDP
%}

% steady state targets
hss = 0.24;                      % target steady-state hours

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% =======================================================================
% Solve for Deterministic Steady-State Values and Implied Parameters
%=========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

steady = modelsteady_tbd(delta, alpha, bbeta, hss, gamma)

% Steady-state quantities:
rss = steady(1,1);                              
iss = steady(2,1); 
yss = steady(3,1);                              
css = steady(4,1); 
kss = steady(5,1);                              
wss= steady(6,1); 
abar = steady(7,1);                                  
psi = steady(8,1); 
lambdass = steady(9,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ===========================================================
% Write Dynare Source File: 'SLMR.mod' 
%=============================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create file
fid = fopen('SLMR.mod','w+');
fprintf(fid, 'var a y i k c h r w lambda; \n');
fprintf(fid, 'varexo eps ;\n');
fprintf(fid, '\n');
fprintf(fid, 'parameters bbeta psi alpha delta gamma abar rho_a sigma_a; \n');
fprintf(fid, '\n');
fprintf(fid, 'bbeta = %3.6f;\n',bbeta);
fprintf(fid, 'psi = %3.6f;\n',psi);
fprintf(fid, 'alpha = %3.6f;\n',alpha);
fprintf(fid, 'delta = %3.6f;\n',delta);
fprintf(fid, 'gamma = %3.6f;\n',gamma);
fprintf(fid, 'rho_a= %3.6f;\n',rho_a);
fprintf(fid, 'sigma_a = %3.6f;\n',sigma_a);
fprintf(fid, 'abar = %3.6f;\n', abar);
fprintf(fid, '\n');

fprintf(fid, 'initval;\n');
fprintf(fid, 'lambda = %3.12f;\n',log(lambdass));
fprintf(fid, 'c = %3.12f;\n',log(css));
fprintf(fid, 'w = %3.12f;\n',log(wss));
fprintf(fid, 'y = %3.12f;\n',log(yss));
fprintf(fid, 'h = %3.12f;\n',log(hss));
fprintf(fid, 'i = %3.12f;\n',log(iss));
fprintf(fid, 'k = %3.12f;\n',log(kss));
fprintf(fid, 'r = %3.12f;\n',log(rss));
fprintf(fid, 'a = %3.12f;\n',abar);
fprintf(fid, 'end;\n');
fprintf(fid, '\n');

fprintf(fid, 'model;\n');
fprintf(fid, '1/exp(c) = exp(lambda);\n');
fprintf(fid, 'gamma = (bbeta/gamma) * (exp(c)/exp(c(+1))) *((1-alpha) * exp(a(+1)) * exp(k)^(-alpha) * exp(h(+1))^alpha + (1-delta));\n');
fprintf(fid, 'psi/(1-exp(h)) = (1/exp(c)) * (alpha * exp(a) * exp(k)^(1-alpha) * exp(h)^(alpha-1));\n');
fprintf(fid, 'exp(a)*k(-1)^(1-alpha)*exp(h)^(alpha) + (1-delta) * exp(k(-1)) - exp(k) * gamma = exp(c);\n');
fprintf(fid, 'a = abar + rho_a * a(-1) + sigma_a * eps;\n');
fprintf(fid, 'exp(y) = exp(a)*(exp(k(-1)))^(1-alpha) * exp(h)^(alpha);\n');
fprintf(fid, 'exp(y) = exp(c) + exp(i);\n');
fprintf(fid, 'exp(w) = alpha * exp(a) * exp(k(-1))^(1-alpha) * exp(h)^(alpha-1);\n');
fprintf(fid, 'exp(r) = (1-alpha) * exp(a) * exp(k(-1))^(-alpha) * exp(h)^alpha;\n');
fprintf(fid, 'end;\n');

fprintf(fid, '\n');
fprintf(fid, 'shocks;\n');
fprintf(fid, 'var eps = 1; \n');
fprintf(fid, 'end;\n');
fprintf(fid, 'steady;\n');
fprintf(fid,'options_.TeX = 1;\n');
fprintf(fid, 'stoch_simul(order=1,graph,hp_filter=1600, irf = 60);\n');
fclose(fid);

%% ===========================================================
% Solve Model using Dynare
%=============================================================
 addpath c:\dynare\4.5.7\matlab    
 dynare SLMR.mod

%%=======================
%  Customized output
%%========================
 printu=1;  %this NEEDs to be updated, indexes have shifted

 if printu==1 
 sdout=(diag(oo_.var(1:9,1:9)))'.^.5;
 siy=sdout(3)/sdout(2); 
 scy=sdout(5)/sdout(2);
 shy=sdout(6)/sdout(2);
 swy=sdout(8)/sdout(2);
 
 disp( '         std/stdy') 
 disp([M_.endo_names(1,:) num2str(sdout(1),3) ])
 disp([ M_.endo_names(2,:) num2str(sdout(2),3)  ])
 disp([M_.endo_names(3,:) num2str(siy,3) ])
 disp([ M_.endo_names(5,:) num2str(scy,2)  ])
 disp([ M_.endo_names(6,:) num2str(shy,2)  ])
 disp([ M_.endo_names(8,:) num2str(swy,2)  ])

 si=sdout(3); 
 sc=sdout(5);
 sh=sdout(6);
 sw=sdout(8);
 
  
 disp( '         std') 
 disp([M_.endo_names(1,:) num2str(sdout(1),3) ])
 disp([ M_.endo_names(2,:) num2str(sdout(2),3)  ])
 disp([M_.endo_names(3,:) num2str(si,3) ])
 disp([ M_.endo_names(5,:) num2str(sc,2)  ])
 disp([ M_.endo_names(6,:) num2str(sh,2)  ])
 disp([ M_.endo_names(8,:) num2str(sw,2)  ])
 end
 
% 
 checku=0;   %set to for checks on steady state
 if checku==1


ssve=[1
yss     
iss     
css        
kss
hss       
rss
wss
lambdass];    

ssvem=exp(oo_.steady_state(1:9));

[ssve  ssvem]
 end;