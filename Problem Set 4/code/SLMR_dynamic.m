function [residual, g1, g2, g3] = SLMR_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(9, 1);
T10 = 1/exp(y(7));
T20 = params(1)/params(5)*exp(y(7))/exp(y(13));
T30 = (1-params(3))*exp(y(12))*exp(y(6))^(-params(3));
T33 = exp(y(8))^params(3);
T34 = T30*T33;
T47 = params(3)*exp(y(3))*exp(y(6))^(1-params(3));
T49 = exp(y(8))^(params(3)-1);
T50 = T47*T49;
T51 = T10*T50;
T55 = exp(y(3))*y(2)^(1-params(3));
T56 = T33*T55;
T74 = exp(y(3))*exp(y(2))^params(3);
T75 = exp(y(8))^(1-params(3));
T76 = T74*T75;
T85 = params(3)*exp(y(3))*exp(y(2))^(1-params(3));
T86 = T49*T85;
T92 = (1-params(3))*exp(y(3))*exp(y(2))^(-params(3));
T93 = T33*T92;
lhs =T10;
rhs =exp(y(11));
residual(1)= lhs-rhs;
lhs =params(5);
rhs =T20*(T34+1-params(4));
residual(2)= lhs-rhs;
lhs =params(2)/(1-exp(y(8)));
rhs =T51;
residual(3)= lhs-rhs;
lhs =T56+(1-params(4))*exp(y(2))-params(5)*exp(y(6));
rhs =exp(y(7));
residual(4)= lhs-rhs;
lhs =y(3);
rhs =params(7)*y(1)+params(8)*x(it_, 1);
residual(5)= lhs-rhs;
lhs =exp(y(4));
rhs =T76;
residual(6)= lhs-rhs;
lhs =exp(y(4));
rhs =exp(y(7))+exp(y(5));
residual(7)= lhs-rhs;
lhs =exp(y(10));
rhs =T86;
residual(8)= lhs-rhs;
lhs =exp(y(9));
rhs =T93;
residual(9)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(9, 14);

  %
  % Jacobian matrix
  %

T149 = exp(y(8))*getPowerDeriv(exp(y(8)),params(3),1);
T159 = exp(y(8))*getPowerDeriv(exp(y(8)),params(3)-1,1);
  g1(1,7)=(-exp(y(7)))/(exp(y(7))*exp(y(7)));
  g1(1,11)=(-exp(y(11)));
  g1(2,12)=(-(T20*T34));
  g1(2,6)=(-(T20*T33*(1-params(3))*exp(y(12))*exp(y(6))*getPowerDeriv(exp(y(6)),(-params(3)),1)));
  g1(2,7)=(-(T20*(T34+1-params(4))));
  g1(2,13)=(-((T34+1-params(4))*params(1)/params(5)*(-(exp(y(7))*exp(y(13))))/(exp(y(13))*exp(y(13)))));
  g1(2,8)=(-(T20*T30*T149));
  g1(3,3)=(-T51);
  g1(3,6)=(-(T10*T49*params(3)*exp(y(3))*exp(y(6))*getPowerDeriv(exp(y(6)),1-params(3),1)));
  g1(3,7)=(-(T50*(-exp(y(7)))/(exp(y(7))*exp(y(7)))));
  g1(3,8)=(-(params(2)*(-exp(y(8)))))/((1-exp(y(8)))*(1-exp(y(8))))-T10*T47*T159;
  g1(4,3)=T56;
  g1(4,2)=(1-params(4))*exp(y(2))+T33*exp(y(3))*getPowerDeriv(y(2),1-params(3),1);
  g1(4,6)=(-(params(5)*exp(y(6))));
  g1(4,7)=(-exp(y(7)));
  g1(4,8)=T55*T149;
  g1(5,1)=(-params(7));
  g1(5,3)=1;
  g1(5,14)=(-params(8));
  g1(6,3)=(-T76);
  g1(6,4)=exp(y(4));
  g1(6,2)=(-(T75*exp(y(3))*exp(y(2))*getPowerDeriv(exp(y(2)),params(3),1)));
  g1(6,8)=(-(T74*exp(y(8))*getPowerDeriv(exp(y(8)),1-params(3),1)));
  g1(7,4)=exp(y(4));
  g1(7,5)=(-exp(y(5)));
  g1(7,7)=(-exp(y(7)));
  g1(8,3)=(-T86);
  g1(8,2)=(-(T49*params(3)*exp(y(3))*exp(y(2))*getPowerDeriv(exp(y(2)),1-params(3),1)));
  g1(8,8)=(-(T85*T159));
  g1(8,10)=exp(y(10));
  g1(9,3)=(-T93);
  g1(9,2)=(-(T33*(1-params(3))*exp(y(3))*exp(y(2))*getPowerDeriv(exp(y(2)),(-params(3)),1)));
  g1(9,8)=(-(T92*T149));
  g1(9,9)=exp(y(9));

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],9,196);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],9,2744);
end
end
end
end
