clear;clc;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for problem 2 of problem set 2. Discretizes AR(1) process and  %
% compute VFI solutions for optimal consumption and savings processes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameter values
r = 1.025725
%r = 1.01
sigma = 0.35;
gamma = 2;
ybar = 1;
beta = 0.97;
rho = 0.6;
var  = sigma^2 * (1 - rho^2)^(1/2);
b = -1; % borrowing constraint
var_process  = var/(1-rho^2);
mean_process = 0;
S = 2500; 
H = 500;
% Set grid parameters
nk =  500 + 1;
kbarstate = 0.5 * (nk - 1) + 1;
kmin = b;
kmax = 25;
kstep = (kmax - kmin)/(nk - 1);

k = [kmin: kstep: kmax]';                        % Grid for current wealth
kp = k;                                    % Grid for next period's wealth


rate_seq =[0, 3, r];
suppl_seq = [H*kmin-1, kmax*H+1];
alpha = 3;
tic;
toc_seq = 0;

    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Markov simulation of y process %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ngridpoints = 9;
ny = ngridpoints;

[P, logy] = quadNorm(ny, mean_process, 0.9, rho);            % Generate transition matrix and y grid
P = P';

y = exp(logy);

n = 25000;
% set parameters for size of chain
nsim = n;
income_state = zeros(1,nsim);
income_state(1) = 5; % initial state to start chain

% preallocate memory to store y realizations chosen by Markov chain
y_markov = zeros(1, nsim);
y_markov(1) = y(income_state(1));

% Markov chains
for i = 2:nsim
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

% plot markov chain simulation

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Value function iteration to compute consumption and savings decisions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

U = (c.^(1 - gamma) - 1 )/(1 - gamma);

% Initial guess for the value function
V0 = zeros(nk, ny);
kpopt_ind = zeros(nk, ny); % matrix to store the index of the optimal kp for each ny
err = 1;

t12 = clock;

% Start iterations to determine optimal value function
while err > 0.0001
    
    for i=1:ny
        % use transition matrix to generate distribution for each y
        dist = P(i, :);
        
        % Calculate the expected value of V0
        EV0 = dist(1) * V0(:, 1);                           % calculate the first element in the sum to get expected value of V0
        for n = 2:ny                                        % start at two since we took care of first element outside of loop
            EV0 = EV0 + dist(n) * V0(:, n);
        end
        
        [V1(:, i), kpopt_ind(:, i)] = max(U(:, :, i) + beta * EV0 * ones(1, nk));      % for every income level y_i pick the max along the k grid
    end
    
    err = max(max(abs(V1 - V0)));
    
    V0 = V1;
    %pause
end


e2 = etime(clock, t12)

    figure(1)
    plot(k, V0(:,5), k, V0(:, 9))
    legend('ybar', 'yhigh')
    xlabel('k');
    ylabel('Value function (R = 1.01)');
    saveas(gcf,sprintf('../Write_up/value_functions_gridpoints%d.png',nk));

% calculate policy function for kp
kp = k(kpopt_ind);

% calculate consumption process
c_policy = zeros(nk, ny);

% calculate consumption
for i = 1:ny
    c_policy(:, i) = y(i)*ones(nk, 1) + r * k - kp(:, i);
end

    figure(2)
    subplot(2, 1, 1);
    plot(k, kp(:, 5),k, kp(:, 9));
    title('kp policy function (R = 1.01)');
    xlabel('k') ;
    ylabel('kp');
    legend('kp(k, ybar)', 'kp(k, yhigh')

    subplot(2, 1, 2);
    plot(k, c_policy(:, 5), k, c_policy(:, 9));
    title('Consumption (R = 1.01)');
    xlabel('k') ;
    ylabel('c');
    legend('c(ybar)', 'c(yhigh)');
  saveas(gcf,sprintf('../Write_up/c_and_k_gridpoints%d.png',nk));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate  economy with given y process and policy functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nhouseholds = 1000;
n = 2500;

%pre-allocate memory for each household
demand_mat = zeros(nhouseholds, (n-500)); % subtract 500 to get far enough away from our initial state assumption

% simulate y process
nsim = n;
income_state = zeros(nhouseholds, nsim);
income_state(:, 1) = 3; % initial state to start chain

% preallocate memory to store y realizations chosen by Markov chain
y_markov = zeros(nhouseholds, nsim);
y_markov(:, 1) = y(income_state(:,1));
kp_sim = zeros(nhouseholds, nsim);
c_sim = zeros(nhouseholds, nsim);

% Markov chains
for i = 2:nsim
    % select row that specifies the distribution associated with current state
    markov_dist = P(income_state(:, i - 1), :); % select distribution for each household in each iteration (i.e. should be a mtrix where each row corresponds to the transition probabilities for that housheold
    
    % calculate cumulative distribution of switching to any state
    cumulative_distribution = cumsum(markov_dist, 2);
    
    % randomly generate probability threshold
    q = rand(1, nhouseholds);
    
    % select new state
    for j = 1:nhouseholds
        index_next_state = find(cumulative_distribution(j, :) > q(j), 1);
        income_state(j, i) = index_next_state;
    end
    
    y_markov(:, i) = y(income_state(:,i));
end

% use simulated y process to generate simulated consumption and savings decisions
k_state = zeros(nhouseholds, nsim);
k_state(:, 1) = 25; % pick an arbitrary initial k state

for i = 1:nsim
    ky_state = [k_state(:, i) income_state(:, i)];
    
    for j = 1:length(ky_state)
        kp_sim(j, i) = kp(ky_state(j, 1), ky_state(j, 2));   % calculate optimal kp for given y and k states
    end
    
    for j = 1:size(kp_sim, 1)
        k_state(j, (i+1)) = find(k == kp_sim(j, i));  %update state of k we are in
    end
    
end

% drop first 500 observations
y_markov = y_markov(:, 501:size(y_markov, 2));
kp_sim = kp_sim(:, 501:size(kp_sim, 2));
aggregate_demand = mean(kp_sim, 1);
sum(aggregate_demand)


figure(3)
subplot(1, 1, 1);
plot(aggregate_demand)
title('Aggregate Demand for Assets R = R*');
ylabel('AD assets') ;
xlabel('period of simulation');
saveas(gcf,'../Write_up/aggregate_demand_equilibrium.png');
