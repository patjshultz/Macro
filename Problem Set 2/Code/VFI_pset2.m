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

tic;
toc_seq = 0;
while abs(rate_seq(end)-rate_seq(end-1)) > 0.0005
    % equivalent result for the following while cond. (H*kstep := max dev. dep. on the economy):
    % (abs(suppl_seq(end)-suppl_seq(end-1)) > H*kstep) || (suppl_seq(end)*suppl_seq(end-1)>0)
    
    
    % Compute momentary utility function
    c = zeros(nk, nk, ny);
    U = zeros(nk, nk, ny);
    
    for j = 1:ny
        for i = 1:nk
            c(:,i,j) = max( (r*k(i) + y(j))*ones(size(kp)) - kp, 0.001);
        end
    end
    
    U = (c.^(1 - 1/sigma) - 1 )/(1 - 1/sigma);
    
    % Initial guess for the value function
    
    V0 = zeros(nk,ny);
    init = repmat(kp,1,size(y,2))+repmat(y,size(kp,1),1);
    V1 = (init.^(1 - 1/sigma) - 1 )/(1 - 1/sigma);
    
    kp_opt_ind = zeros(nk,ny);
    
    err = 1;
    
    % Start iterations to determine optimal value function
    while err > 0.0001
        for i=1:ny
            [V1(:,i), kp_opt_ind(:,i)] = max(U(:,:,i) + beta*V0(:,:)*P(i,:)'*ones(1,nkp));
            % for every income level y_i pick the max kp along the k grid
            % use transition matrix to find expected next periods value fct. which
            % is effectively a probability weighted average value function (expected value)
        end
        % for each level of y_i we add last round's optimal kp picks (which is
        % V0, so it's mutliplied by beta etc.) which is constant per level of k
        % (thus upscaled in the kp dimension
        temp = (V1 - V0).^2;
        err = sum(temp(:))/numel(V1);
        V0 = V1; %we now have a matrix of kp values on a k*y grid
    end
    
    nybar = 0.5*(ny+1);
    
    if r == 1.01
        figure(1)
        subplot(2,1,1)
        plot(k,V1(:,nybar))
        title(['\mu_y =' num2str(round(y(nybar),3))])
        xlabel('k')
        ylabel('V(\cdot)')
        subplot(2,1,2)
        plot(k,V1(:,end))
        title(['y = ' num2str(round(y(end),3))])
        xlabel('k')
        ylabel('V(\cdot)')
        %saveas(figure(1),[pwd '/02_Graphs/II_1a_TwoLevels.png'])
    end
    
    
    %Task (2c)
    %% Extract c* and kp*
    
    kp_opt = kp*ones(1,ny);
    kp_opt = kp_opt(kp_opt_ind);
    
    c_opt = transpose(max( (ones(ny,1)*r*k' + y'*ones(1,nk)) - kp_opt', 0.001));
    
    if r == 1.01
        figure(2)
        subplot(2,1,1)
        plot(k,c_opt(:,nybar), k,c_opt(:,end))
        title('Optimal Consumption')
        legend(strcat('\mu_y = ', num2str(round(y(nybar),1))), ...
            strcat('y = ', num2str(round(y(end),1))),'Location','SouthEast')
        xlabel('k')
        ylabel('c')
        subplot(2,1,2)
        plot(k,kp_opt(:,nybar), k,kp_opt(:,end))
        title('Optimal Investment Policy')
        legend(strcat('\mu_y = ', num2str(round(y(nybar),1))), ...
            strcat('y = ', num2str(round(y(end),1))),'Location','SouthEast')
        xlabel('k')
        ylabel('kp')
        
        saveas(figure(2),[pwd '/02_Graphs/II_a_kpAndcPerK.png'])
    end
    
    %Task (2b)
    %% Simulate the economy
    chain_length = S;
    chain = zeros(H,chain_length);
    chain(:,1) = ones(H,1)*(ny+1)/2; % initial state number to start chain
    %chain(1:3,1) = [1 4 9]
    
    % preallocate memory to store y realizations chosen by Markov chain
    y_markov = zeros(H, chain_length);
    grid_lny_H = repmat(grid_lny,H,1);
    y_markov(:,1) = grid_lny_H(1,chain(:,1));
    
    for i=2:chain_length
        distribution = P(chain(:,i-1), :);
        cumulative_distribution = cumsum(distribution,2);
        j = rand(H,1);
        chain(:,i) = max(1,L-sum(cumulative_distribution>j,2));
        % L the number of states in y, sum() gives row-wise no of cols in cdf
        % matrix P where j is exceeded, so L-sum() is col-number in P, max is
        % to handle cases where even first row exceeds j - avoiding col=0
        y_markov(:,i) = grid_lny_H(1,chain(:,i));
    end
    y_chain = exp(y_markov);
    k_chain = zeros(H,S+1); % requires one column more as policy concerns kp
    k_chain(:,1) = repmat(k(0.5*(nk+1)),H,1);
    c_chain = zeros(H,S);
    
    %k_chain(1:3,1) = k(1:3);
    
    for j = 1:S
        [lia, loc] = ismember(k_chain(:,j),k);
        [lib, locy] = ismember(y_chain(:,j),y);
        k_chain(:,j+1) = diag(kp_opt(loc,locy));
        
        c_chain(:,j) = (max( (r*k_chain(:,j) + y_chain(:,j) ) - k_chain(:,j+1), 0.001));
    end
    % cut out first 20% of obs for stable economy (0.2*S=500, from 501
    % onwards and accounting for kp being shifted by +1 makes +2)
    
    k_stable = k_chain(:,round(0.20*S,0)+2:end);
    
    supply = sum(k_stable,1);
    suppl_seq = [suppl_seq mean(supply)];
    
    if suppl_seq(alpha)*suppl_seq(alpha-1)<0 %eq. in the middle of n and n-1
        r = 0.5*(rate_seq(alpha)+rate_seq(alpha-1)); %avg
    elseif suppl_seq(alpha)*suppl_seq(alpha-2)<0 %eq. in the middle of n and n-2
        r = 0.5*(rate_seq(alpha)+rate_seq(alpha-2)); %avg
    elseif suppl_seq(alpha)<0 %avg with last opposing sign rate
        r = 0.5*(rate_seq(alpha)+rate_seq(find(suppl_seq>0,1,'last')));
    else r = 0.5*(rate_seq(alpha)+rate_seq(find(suppl_seq<0,1,'last')));
    end
    
    rate_seq = [rate_seq r];
    alpha = alpha +1
    toc1 = toc
    toc_seq = [toc_seq toc1];
end
eq_rate = r;

