function F = modelsteady_tbd(delta, alpha, bbeta, hss, gamma, gy)

% Write down equations to compute steady state values for the deterministic model
% Make sure to solve for 9 variables (rss iss yss css wss abar psi
% lambdass)

rss = 1/bbeta - (1 - delta);
kss = (rss / ((1 - alpha) * hss^alpha))^(-1/alpha);
iss = (gamma - 1 + delta) * kss;
yss = kss^(1-alpha) * hss^(alpha);
gss = gy * yss;
css = yss - iss - gss;
wss = alpha * (kss^(1-alpha)) * hss ^(alpha-1);
abar = 0;
psi =  alpha * (yss/css)*((1-hss)/hss);
lambdass = 1/css;


% Define output vector
%     1   2   3   4   5   6   7     8     9      10  
F = [rss iss yss css kss wss abar psi  lambdass gss]';
