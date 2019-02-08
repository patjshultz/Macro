clear;clc;
% VFI
% This routine outlines the solution procedure for a simple consumer problem using 
% some popular variants of basic value function iteration technique

% Set parameter values
r = 1.04;
sigma = 1/2;
ybar = 1;
beta = 0.95;

% Set grid parameters
nk =  50 + 1;
kbarstate = 0.5 * (nk - 1) + 1;
kbar = 5;
kmin = 4;
kmax = 6;
kstep = (kmax - kmin)/(nk - 1);

k = [kmin: kstep: kmax]';                       % Grid for current wealth
kp = k;                                         % Grid for next period's wealth

summary_stats = zeros(3, 6);                    % set up matrix to store summary statistics
col = 1;

for ngridpoints = [5 9]
    ny = ngridpoints;
    
    %[logy,P,d]=tauchen1(ny, 1, 0.95 , 0.1, 4); % choose 4 std dev for grid of 9 to have each point at one std dev

    [P, logy] = quadNorm(ny, 1, 0.1, 0.99);            % Generate transition matrix and y grid
    P = P';
    
    y = exp(logy);
    % y = logy to check if rho = 0.95
    
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

        % calculate summary statistics of our markov chain

        summary_stats(1, col) = mean(y_markov);                    % store mean
        [ac] = autocorr(y_markov, 'NumLags', 1);
        summary_stats(2, col) =  ac(2);      % store autocorr
        summary_stats(3, col) = sqrt(var((y_markov)));            % store var
        col = col + 1;
    end

end


% model = arima('Constant',0.05,'AR',{0.95},'Variance',sqrt(0.1));
% rng('default')
% ysim = simulate(log(model),1000);
% 
% figure
% plot(ysim)
% xlim([0,1000])
% title('Simulated AR(1) Process')

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
    % use transition matrix to generate distribution for each y
    dist = P(i, :);                                    
    
    % Calculate the expected value of V0
    EV0 = dist(1) * V0(:, 1);                           % calculate the first element in the sum to get expected value of V0
    for n = 2:ny                                        % start at two since we took care of first element outside of loop
        EV0 = EV0 + dist(n) * V0(:, n);
    end
    
    V1(:, i) = max(U(:, :, i) + beta * EV0 * ones(1, nk));      % for every income level y_i pick the max along the k grid
  end
 
  err = max(max(abs(V1 - V0)))
  
  V0 = V1;
  %pause
end


e2 = etime(clock, t12)

figure(2)
plot(k, V0(:,5), k, V0(:, 6))
legend('ybar', 'ybar + 1*sd')

