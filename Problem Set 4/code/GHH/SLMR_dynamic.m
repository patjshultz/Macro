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
T18 = exp(y(7))-params(2)*exp(y(8))^(1+params(9))/(1+params(9));
T29 = params(1)/params(5)*exp(y(13))/exp(y(11));
T39 = (1-params(3))*exp(y(12))*exp(y(6))^(-params(3));
T40 = exp(y(8))^params(3);
T41 = T39*T40;
T55 = params(3)*exp(y(3))*exp(y(2))^(1-params(3));
T57 = exp(y(8))^(params(3)-1);
T58 = T55*T57;
T61 = exp(y(3))*y(2)^(1-params(3));
T62 = T40*T61;
T79 = exp(y(3))*exp(y(2))^params(3);
T80 = exp(y(8))^(1-params(3));
T81 = T79*T80;
T94 = (1-params(3))*exp(y(3))*exp(y(2))^(-params(3));
T95 = T40*T94;
lhs =1/T18;
rhs =exp(y(11));
residual(1)= lhs-rhs;
lhs =params(5);
rhs =T29*(T41+1-params(4));
residual(2)= lhs-rhs;
lhs =params(2)*exp(y(8))^params(9);
rhs =T58;
residual(3)= lhs-rhs;
lhs =T62+(1-params(4))*exp(y(2))-params(5)*exp(y(6));
rhs =exp(y(7));
residual(4)= lhs-rhs;
lhs =y(3);
rhs =params(7)*y(1)+params(8)*x(it_, 1);
residual(5)= lhs-rhs;
lhs =exp(y(4));
rhs =T81;
residual(6)= lhs-rhs;
lhs =exp(y(4));
rhs =exp(y(7))+exp(y(5));
residual(7)= lhs-rhs;
lhs =exp(y(10));
rhs =T58;
residual(8)= lhs-rhs;
lhs =exp(y(9));
rhs =T95;
residual(9)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(9, 14);

  %
  % Jacobian matrix
  %

T108 = (-(T57*params(3)*exp(y(3))*exp(y(2))*getPowerDeriv(exp(y(2)),1-params(3),1)));
T140 = exp(y(8))*getPowerDeriv(exp(y(8)),params(3),1);
T149 = T55*exp(y(8))*getPowerDeriv(exp(y(8)),params(3)-1,1);
  g1(1,7)=(-exp(y(7)))/(T18*T18);
  g1(1,8)=params(2)*exp(y(8))*getPowerDeriv(exp(y(8)),1+params(9),1)/(1+params(9))/(T18*T18);
  g1(1,11)=(-exp(y(11)));
  g1(2,12)=(-(T29*T41));
  g1(2,6)=(-(T29*T40*(1-params(3))*exp(y(12))*exp(y(6))*getPowerDeriv(exp(y(6)),(-params(3)),1)));
  g1(2,8)=(-(T29*T39*T140));
  g1(2,11)=(-((T41+1-params(4))*params(1)/params(5)*(-(exp(y(11))*exp(y(13))))/(exp(y(11))*exp(y(11)))));
  g1(2,13)=(-(T29*(T41+1-params(4))));
  g1(3,3)=(-T58);
  g1(3,2)=T108;
  g1(3,8)=params(2)*exp(y(8))*getPowerDeriv(exp(y(8)),params(9),1)-T149;
  g1(4,3)=T62;
  g1(4,2)=(1-params(4))*exp(y(2))+T40*exp(y(3))*getPowerDeriv(y(2),1-params(3),1);
  g1(4,6)=(-(params(5)*exp(y(6))));
  g1(4,7)=(-exp(y(7)));
  g1(4,8)=T61*T140;
  g1(5,1)=(-params(7));
  g1(5,3)=1;
  g1(5,14)=(-params(8));
  g1(6,3)=(-T81);
  g1(6,4)=exp(y(4));
  g1(6,2)=(-(T80*exp(y(3))*exp(y(2))*getPowerDeriv(exp(y(2)),params(3),1)));
  g1(6,8)=(-(T79*exp(y(8))*getPowerDeriv(exp(y(8)),1-params(3),1)));
  g1(7,4)=exp(y(4));
  g1(7,5)=(-exp(y(5)));
  g1(7,7)=(-exp(y(7)));
  g1(8,3)=(-T58);
  g1(8,2)=T108;
  g1(8,8)=(-T149);
  g1(8,10)=exp(y(10));
  g1(9,3)=(-T95);
  g1(9,2)=(-(T40*(1-params(3))*exp(y(3))*exp(y(2))*getPowerDeriv(exp(y(2)),(-params(3)),1)));
  g1(9,8)=(-(T94*T140));
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
