clear;clc;
% VFI
% This routine outlines the solution procedure for a simple consumer problem using 
% some popular variants of basic value function iteration technique

% Set parameter values
r = 1.04;
sigma = 1/2;
ybar = 1;
beta = 0.95;

% Choose solution technique
% Alternatives are
% discrete(1)
%technique = discrete1

% Set grid parameters
nk =  50 + 1;
kbar = 5;
kmin = 0;
kmax = 6;
kstep = (kmax - kmin)/(nk - 1);

k = [kmin: kstep: kmax]';                       % Grid for current wealth
kp = k;                                         % Grid for next period's wealth

ny = 9;
[P, y] = quadNorm(ny, 1, 0.1, 0.95);             % Generate transition matrix and y grid
P = P';
y = exp(y);
% set parameters for size of chain
nsim = 10000;
income_state = zeros(1,nsim);
income_state(1)=5; % initial state to start chain

% preallocate memory to store y realizations chosen by Markov chain
y_markov = zeros(1, nsim);
y_markov(1) = y(income_state(1));

% Markov chains
for i=2:nsim
    % select row that specifies the distribution associated with current state 
    markov_dist = P(income_state(i - 1), :); 
    
    % calculate cumulative distribution of switching to any state
    cumulative_distribution = cumsum(markov_dist);
    
    % randomly generate probability threshold
    q = rand();
    
    % select new state
    income_state(i) = find(cumulative_distribution > q, 1);
    y_markov(i) = y(income_state(i));
end

% calculate summary statistics of our markov chain
mean(y_markov)
autocorr(y_markov);
sqrt(var(y_markov))
y_sigma = mean(y_markov + sqrt(var(y_markov)));

% Compute accurate solution to test quality of the procedure
ho = ybar/(1 - 1/r);
co = (r*k + ho)*(1 - (beta^sigma)*(r^(sigma - 1)));
Ve = ((co.^(1 - 1/sigma))/(1 - ((r*beta)^sigma)*beta) - 1/(1 - beta))/(1 - 1/sigma);

% Compute momentary utility function  
c = zeros(nk, nk, ny);
U = zeros(nk, nk, ny);

% this loop calculates all consumptions implied by k's in the k grid 
% if the implied consumption is negative, change it to 0.001
for i = 1:nk
    for j = 1:ny
        c(:,i, j) = max((r*k(i) + y(j))*ones(size(kp)) - kp, 0.001); 
    end
end  

U = (c.^(1 - 1/sigma) - 1 )/(1 - 1/sigma);

% Initial guess for the value function   
V0 = zeros(nk, ny);
%BV0 = zeros(nk, nk, ny);
err = 1;

t12 = clock;

% Start iterations to determine optimal value function
while err > 0.0001
   
  for i=1:ny
      V1(:,i) = max(U(:,:,i) + beta*V0(:,i)*ones(1,nk));    % for every income level y_i pick the max along the k grid
  end
 
  err = max(max(abs(V1 - V0)))
  
  V0 = V1;
  %pause
end


e2 = etime(clock, t12)

figure(2)
plot(k, V0(:,5), k, V0(:, 6))
legend('ybar', 'ybar + 1*sd')

