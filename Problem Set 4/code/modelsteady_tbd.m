function F = modelsteady(delta, gamma, alpha, bbeta, kss, hss, ies)

% Write down equations to compute steady state values for the deterministic model
% Make sure to solve for 9 variables (rss iss yss css wss abar theta eta
% lambdass)

rss = 1/bbeta - (1 - delta);
%missing remaining equations

% Define output vector
%     1   2   3   4   5    6    7    8      9      
F = [rss iss yss css wss abar theta eta lambdass]';
