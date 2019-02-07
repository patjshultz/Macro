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

% Compute accurate solution to test quality of the procedure
ho = ybar/(1 - 1/r);
co = (r*k + ho)*(1 - (beta^sigma)*(r^(sigma - 1)));
Ve = ((co.^(1 - 1/sigma))/(1 - ((r*beta)^sigma)*beta) - 1/(1 - beta))/(1 - 1/sigma);

% Compute momentary utility function  
c = zeros(nk, nk);
U = zeros(nk, nk);

for i = 1:nk
   c(:,i) = max( (r*k(i) + ybar)*ones(size(kp)) - kp, 0.001);
end  

U = ( c.^(1 - 1/sigma) - 1 )/(1 - 1/sigma);

% Method 1: Discrete state space
% Case 1: Linear function of wealth
% Initial guess for the value function 
V0 = zeros(nk,1);

err = 1;

t11 = clock;

% Start iterations to determine optimal value function
while err > 0.0001
	V1 = max(U + beta*V0*ones(1, nk));    
   
   err = sqrt((V1 - V0')*(V1' - V0));
   
   V0 = V1';
   
end

figure(1)
e1 = etime(clock, t11)
subplot(2,2,1)
plot(k, V0, k, Ve)

% Case 2: Same functional form as the objective function
% Initial guess for the value function   
V0 = (k.^(1 - 1/sigma))/(1 - beta);

err = 1;

t12 = clock;

% Start iterations to determine optimal value function
while err > 0.0001
   U + beta*V0*ones(1, nk);

   V1 = max(U + beta*V0*ones(1, nk));   
   
   err = sqrt((V1 - V0')*(V1' - V0));
   
   V0 = V1';
   %pause
end

e2 = etime(clock, t12)
subplot(2,2,2)
plot(k, V0, k, Ve)

% Method 2: Polynomial Approximation (Second Order)
% Case 1: Ordinary polynomials, linear function of wealth
b0 = [0 1 0 0]';
K = [ones(size(kp)) kp kp.^2 kp.^3];
V0 = K*b0;
X = inv(K'*K)*K';

err = 1;

t21 = clock;

% Start iterations to determine optimal value function
while err > 0.0001
   
	V1 = max(U + beta*V0*ones(1, nk));    
   
   b1 = X*V1';
   err = sqrt((b1 - b0)'*(b1 - b0));
   
   V0 = K*b1;
   b0 = b1;
end

e3 = etime(clock, t21)

subplot(2,2,3)

plot(k, V0, k, Ve)

% Case 2: Chebyshev polynomials, linear function of wealth
b0 = [0 1 0 0]';
T0 = ones(size(kp));
T1 = kp;
T2 = 2*kp.*T1 - T0;						% General recursion form for high order terms
T3 = 2*kp.*T2 - T1;						% General recursion form for high order terms
KT = [T0 T1 T2 T3];
V0 = KT*b0;
X = inv(KT'*KT)*KT';

err = 1;

t22 = clock;

% Start iterations to determine optimal value function
while err > 0.0001
   
  V1 = max(U + beta*V0*ones(1, nk));    
   
   b1 = X*V1';
   err = sqrt((b1 - b0)'*(b1 - b0));
   
   V0 = KT*b1;
   b0 = b1;
end

e4 = etime(clock, t22)

subplot(2,2,4)

plot(k, V0, k, Ve)

%i = find(U + beta*V0*ones(1, nk) == max(U + beta*V0*ones(1, nk)))
% compute optimal decision rules
policy = (U + beta*V0*ones(1, nk) == (max(U + beta*V0*ones(1, nk)))'*ones(1,nk));
Kprime = policy*k;

figure(2)
plot(k, Kprime)