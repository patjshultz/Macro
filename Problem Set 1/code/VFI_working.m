%% Patrick's & Felix' Code
%Task (1a) 

ny = 9 % number of grid points
bar_lny = 1 % unconditional mean of lny
sigma_eps = 0.1 % unconditional std of lny
rho_lny = 0.95 % autocorrelation coefficient

[grid_lny,P,d]=tauchen1(ny,bar_lny,rho_lny,sigma_eps,4); % choose 4 std dev for grid of 9 to have each point at one std dev
    %[P, grid_lny] = Quadnorm(points, bar_lny, sigma_lny, rho_lny)
    %P = P' % to comply with lecture notation where rows sum up to 1

%Task (1b)
chain_length = 10000;
chain = zeros(1,chain_length);
chain(1)=grid_lny((ny+1)/2); % initial state to start chain

% preallocate memory to store y realizations chosen by Markov chain
y_markov = zeros(1, chain_length);
y_markov(1) = grid_lny(chain(1));

for i=2:chain_length
    % select column that specifies the distribution associated with current state 
    distribution = P(chain(i-1), :); 
    % calculate cumulative distribution of switching to any state
    cumulative_distribution = cumsum(distribution);
    % randomly generate threshold
    j = rand();
    chain(i) = find(cumulative_distribution>j,1);
    y_markov(i) = grid_lny(chain(i));
end

figure(1)
subplot(2,1,1)
plot(y_markov)
subplot(2,1,2)
autocorr(y_markov)
mean(y_markov)
sqrt(var(y_markov))

%% Joao's Code (adjusted)
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
nk = 6 + 1;
kbar = 5;
kmin = 4;
kmax = 6;
kstep = (kmax - kmin)/(nk - 1);

k = [kmin: kstep: kmax]';                       % Grid for current wealth
kp = k;                                         % Grid for next period's wealth

y = exp(grid_lny);
yp = y;

% Compute momentary utility function  
c = zeros(nk, nk, ny);
U = zeros(nk, nk, ny);

for j = 1:ny
for i = 1:nk
   c(:,i,j) = max( (r*k(i) + y(j))*ones(size(kp)) - kp, 0.001); %max to avoid 0 or neg c ("o/w agent dies")
end  
end

U = ( c.^(1 - 1/sigma) - 1 )/(1 - 1/sigma);

% Method 1: Discrete state space
% Initial guess for the value function

V0 = zeros(nk,ny);

err = 1;

% Start iterations to determine optimal value function
while err > 0.0001
    for i=1:ny
    V1(:,i) = max(U(:,:,i) + beta*V0(:,i)*ones(1,7));    % for every income level y_i pick the max along the k grid
    end
    
    temp = (V1 - V0).^2;
    err = sum(temp(:))/numel(V1);
    V0 = V1;
end


figure(2)
subplot(2,1,1)
plot(V1(:,0.5*(nk+1)),k)
subplot(2,1,2)
plot(V1(:,0.5*(nk+1)+1),k) %Tauchen set such that each pin of discretization
%corresponds to 1 std dev, therefore +1 row from the mean is =1 std dev