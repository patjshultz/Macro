%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: [Z,Zprob] = mstauchen(N,mu,rho,sigma,Pi,m)
%
% This function produces a discrete markov chain whose sample paths
% approximate those of the markov switching process:
%           z(t+1) = mu_s + rho_s*z(t) + sigma_s*eps(t+1)
% where eps are normally distributed with variance one. This procedure is an
% implementation of George Tauchen's algorithm described in Ec. Letters 20 
% (1986) 177-181. This code is partly based on the function tauchen.m 
% written by Martin Flodén. Modifications to the case of a markov switching
% model were added. 
%
% INPUTS:
%   N:      scalar, number of nodes for Z
%   mu:     vector, intercepts of process
%   rho:    vector, persistence parameters
%   sigma:  vector, std. devs. of epsilons
%   Pi:     transition matrix for the regimes
%   m:      max +- std. devs.
%
% OUTPUTS:
%   Z:      (N*r)*1 vector, nodes for Z. First N elements refer to
%            regime 1, second N elements to regime 2, etc.
%   Zprob:  (N*r)*(N*r) matrix, transition probabilities 
%
%
% Author:   Jan Hannes Lang
% Date:     15.6.2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Z,Zprob] = mstauchen(N,mu,rho,sigma,Pi,m)

% Number of regimes
r     = size(Pi,1);

% Total number of grid point
NT    = N*r;  

% Vector of nodes in each regimes
Z     = zeros(N,1);

% Shell for joint transition probabilities of regimes and process realizations
Zprob = zeros(NT,NT);

% Maximum grid point
Z(N)  = max(mu./(1-rho) + m*sigma./sqrt(1 - rho.^2));

% Minimum grid point
Z(1)  = min(mu./(1-rho) - m*sigma./sqrt(1 - rho.^2));

% Equidistant space between each grid point
zstep = (Z(N) - Z(1)) / (N - 1);

% Relabel the intercept term
a     = mu;

% Create the remaining grid points
for i=2:(N-1)
    Z(i) = Z(1) + zstep * (i - 1);
end 

% Replicate the grid vector r times
Z     = repmat(Z,r,1);

% Loop to fill the cells of the transition matrix
for i = 1:r                 % Loops over the current regime
    for j = 1:r             % Loops over the future regime
        for k = 1:N         % Loops over the current grid point
            for l = 1:N     % Loops over the future grid point
                % Fill first grid point cells
                if l == 1
                    Zprob(k+N*(i-1),l+N*(j-1)) = Pi(i,j)*cdf_normal((Z(l+N*(j-1)) - a(j) - rho(j)*Z(k+N*(i-1)) + zstep/2) / sigma(j));
                % Fill last grid point cells
                elseif l == N
                    Zprob(k+N*(i-1),l+N*(j-1)) = Pi(i,j)*(1 - cdf_normal((Z(l+N*(j-1)) - a(j) - rho(j)*Z(k+N*(i-1)) - zstep/2) / sigma(j)));
                % Fill all other grid point cells
                else
                    Zprob(k+N*(i-1),l+N*(j-1)) = Pi(i,j)*(cdf_normal((Z(l+N*(j-1)) - a(j) - rho(j)*Z(k+N*(i-1)) + zstep/2) / sigma(j)) - ...
                                 cdf_normal((Z(l+N*(j-1)) - a(j) - rho(j)*Z(k+N*(i-1)) - zstep/2) / sigma(j)));
                end
            end
        end
    end
end

% Sub function to create value of a normal distribution
function c = cdf_normal(x)
    c = 0.5 * erfc(-x/sqrt(2));