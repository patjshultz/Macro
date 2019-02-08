clear;clc;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for problem 2 of problem set 1. Discretizes AR(1) process and  %
% compute VFI solutions for optimal consumption and savings processes % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set parameter values
r = 1.04;
sigma = 1/2;
ybar = 1;
beta = 0.95;

% Set grid parameters
nk =  1000 + 1;
kbarstate = 0.5 * (nk - 1) + 1;
kbar = 5;
kmin = 4;
kmax = 6;
kstep = (kmax - kmin)/(nk - 1);

k = [kmin: kstep: kmax]';                       % Grid for current wealth
kp = k;                                         % Grid for next period's wealth

markov_summary_stats = zeros(3, 6);                    % set up matrix to store summary statistics
ar1_summary_stats = zeros(3, 6);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Markov simulation of y process %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
col = 1;

for ngridpoints = [5 9]
    ny = ngridpoints;
    
    [P, logy] = quadNorm(ny, 1, 0.1, 0.95);            % Generate transition matrix and y grid
    P = P';
    
    %y = exp(logy);
    y = logy
    
    for n = [1000 5000 100000]
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
        
        % calculate summary statistics of our markov chain and save in
        % summary stats matrix
        markov_summary_stats(1, col) = mean(y_markov);
        [ac] = autocorr(y_markov, 'NumLags', 1);
        markov_summary_stats(2, col) =  ac(2);
        markov_summary_stats(3, col) = sqrt(var((y_markov)));
        
        % Simulate AR(1) process
        pd = makedist('Normal');
        truncated_dist = truncate(pd,-4,4);
        epsilon = random(pd, [1, nsim]);
        logy_sim = zeros(1,nsim);
        logy_sim(1,1) = 1;
        
        for i = 2:nsim
            logy_sim(1,i) = 0.05 + 0.95*logy_sim(1,i-1) + 0.1*epsilon(1,i);
        end
        
        % plot AR(1) process
        subplot(2,1,1)
        plot(logy_sim)
        title('Sample Path')
        subplot(2,1,2)
        autocorr(logy_sim)
        
        ar1_summary_stats(1, col) = mean(logy_sim);
        [ac] = autocorr(logy_sim, 'NumLags', 1);
        ar1_summary_stats(2, col) =  ac(2);
        ar1_summary_stats(3, col) = sqrt(var((logy_sim)));
        col = col + 1;
    end
    
end


y = exp(logy);
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

U = (c.^(1 - 1/sigma) - 1 )/(1 - 1/sigma);

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
    
    err = max(max(abs(V1 - V0)))
    
    V0 = V1;
    %pause
end


e2 = etime(clock, t12)

figure(2)
plot(k, V0(:,5), k, V0(:, 6))
legend('ybar', 'ybar + 1*sd')
saveas(gcf, 'value_functions.png')

% calculate policy function for kp
kp = k(kpopt_ind);

% calculate consumption process
c_policy = zeros(nk, ny);

% calculate consumption
for i = 1:ny
    c_policy(:, i) = y(i)*ones(nk, 1) + r * k - kp(:, i)
end

subplot(2, 1, 1);
plot(k, kp(:, 5),k, kp(:, 6));
title('kp policy function');
xlabel('k') ;
ylabel('kp');
legend('kp(k, ybar)', 'kp(k, ybar + \sigma_y)')
mesh(kp)

subplot(2, 1, 2);
plot(k, c_policy(:, 5), k, c_policy(:, 6));
title('Consumption');
xlabel('k') ;
ylabel('c');
legend('c(ybar)', 'c(ybar + \sigma_y')
saveas(gcf, 'c_and_k.png')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate  economy with given y process and policy functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j = 1;
simulation_summary_stats = zeros(3, 2);

for n = [5500, 10500]
    % simulate y process
    nsim = n;
    income_state = zeros(1,nsim);
    income_state(1) = 5; % initial state to start chain
    
    % preallocate memory to store y realizations chosen by Markov chain
    y_markov = zeros(1, nsim);
    y_markov(1) = y(income_state(1));
    kp_sim = zeros(1, nsim);
    c_sim = zeros(1, nsim);
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
    
    % use simulated y process to generate simulated consumption and savings decisions
    k_state = 25; % pick an arbitrary initial k state
    for i = 1:nsim
        y_state = income_state(i);                           % update y state
        kp_sim(i) = kp(k_state, y_state);                    % calculate optimal kp for given y and k states
        c_sim(i) = y(y_state) + r * k(k_state) - kp_sim(i) ;     % calculate optimal c for given y and k states
        k_state = find(k == kp_sim(i));                      %update state of k we are in
    end
    
    % drop first 500 observations
    y_markov = y_markov(501:length(y_markov));
    kp_sim = kp_sim(501:length(kp_sim));
    c_sim = c_sim(501:length(c_sim));
    
    % calculate consumption growth 
    c_growth = log(c_sim(2:length(c_sim))) - log(c_sim(1:length(c_sim)-1));
    
    subplot(2, 2, 1);
    plot(y_markov);
    title('Simulated y process')
    
    subplot(2, 2, 2);
    plot(kp_sim);
    title('Simulated kp process');
    
    subplot(2,2,3);
    plot(c_sim);
    title('Simulated c process');
    
    subplot(2,2,4);
    plot(c_growth);
    title('Simulated c growth');
    saveas(gcf, 'simulation.png')

    simulation_summary_stats(1, j) = mean(c_growth);
    simulation_summary_stats(2, j) = std(c_growth);
    acf = autocorr(c_growth, 'NumLags', 1);
    simulation_summary_stats(3, j) = acf(2);
    
    j = j+ 1;
end
