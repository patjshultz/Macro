function F = modelsteady_GHH(delta, alpha, bbeta, hss, gamma, nu)

% Write down equations to compute steady state values for the deterministic model
% Make sure to solve for 9 variables (rss iss yss css wss abar psi
% lambdass)

rss = 1/bbeta - (1 - delta);
kss = (rss / ((1 - alpha) * hss^alpha))^(-1/alpha);
iss = (gamma - 1 + delta) * kss;
yss = kss^(1-alpha) * hss^(alpha);
css = yss - iss;
wss = alpha * (kss^(1-alpha)) * hss ^(alpha-1);
abar = 0;
psi =  alpha * (yss)/(hss^(1+nu));
lambdass = 1/(css - psi * (hss^(1.5)/1.5));


% Define output vector
%     1   2   3   4   5   6   7     8     9        
F = [rss iss yss css kss wss abar psi  lambdass]';
