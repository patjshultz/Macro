function [residual, g1, g2, g3] = SLMR_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 9, 1);

%
% Model equations
%

T18 = exp(y(5))-params(2)*exp(y(6))^(1+params(9))/(1+params(9));
T25 = params(1)/params(5);
T35 = (1-params(3))*exp(y(1))*exp(y(4))^(-params(3));
T36 = exp(y(6))^params(3);
T37 = T35*T36;
T47 = params(3)*exp(y(1))*exp(y(4))^(1-params(3));
T49 = exp(y(6))^(params(3)-1);
T50 = T47*T49;
T53 = exp(y(1))*y(4)^(1-params(3));
T54 = T36*T53;
T70 = exp(y(1))*exp(y(4))^params(3);
T71 = exp(y(6))^(1-params(3));
T72 = T70*T71;
lhs =1/T18;
rhs =exp(y(9));
residual(1)= lhs-rhs;
lhs =params(5);
rhs =T25*(T37+1-params(4));
residual(2)= lhs-rhs;
lhs =params(2)*exp(y(6))^params(9);
rhs =T50;
residual(3)= lhs-rhs;
lhs =T54+exp(y(4))*(1-params(4))-params(5)*exp(y(4));
rhs =exp(y(5));
residual(4)= lhs-rhs;
lhs =y(1);
rhs =y(1)*params(7)+params(8)*x(1);
residual(5)= lhs-rhs;
lhs =exp(y(2));
rhs =T72;
residual(6)= lhs-rhs;
lhs =exp(y(2));
rhs =exp(y(5))+exp(y(3));
residual(7)= lhs-rhs;
lhs =exp(y(8));
rhs =T50;
residual(8)= lhs-rhs;
lhs =exp(y(7));
rhs =T37;
residual(9)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(9, 9);

  %
  % Jacobian matrix
  %

T94 = T36*(1-params(3))*exp(y(1))*exp(y(4))*getPowerDeriv(exp(y(4)),(-params(3)),1);
T101 = (-(T49*params(3)*exp(y(1))*exp(y(4))*getPowerDeriv(exp(y(4)),1-params(3),1)));
T123 = exp(y(6))*getPowerDeriv(exp(y(6)),params(3),1);
T132 = T47*exp(y(6))*getPowerDeriv(exp(y(6)),params(3)-1,1);
  g1(1,5)=(-exp(y(5)))/(T18*T18);
  g1(1,6)=params(2)*exp(y(6))*getPowerDeriv(exp(y(6)),1+params(9),1)/(1+params(9))/(T18*T18);
  g1(1,9)=(-exp(y(9)));
  g1(2,1)=(-(T25*T37));
  g1(2,4)=(-(T25*T94));
  g1(2,6)=(-(T25*T35*T123));
  g1(3,1)=(-T50);
  g1(3,4)=T101;
  g1(3,6)=params(2)*exp(y(6))*getPowerDeriv(exp(y(6)),params(9),1)-T132;
  g1(4,1)=T54;
  g1(4,4)=exp(y(4))*(1-params(4))+T36*exp(y(1))*getPowerDeriv(y(4),1-params(3),1)-params(5)*exp(y(4));
  g1(4,5)=(-exp(y(5)));
  g1(4,6)=T53*T123;
  g1(5,1)=1-params(7);
  g1(6,1)=(-T72);
  g1(6,2)=exp(y(2));
  g1(6,4)=(-(T71*exp(y(1))*exp(y(4))*getPowerDeriv(exp(y(4)),params(3),1)));
  g1(6,6)=(-(T70*exp(y(6))*getPowerDeriv(exp(y(6)),1-params(3),1)));
  g1(7,2)=exp(y(2));
  g1(7,3)=(-exp(y(3)));
  g1(7,5)=(-exp(y(5)));
  g1(8,1)=(-T50);
  g1(8,4)=T101;
  g1(8,6)=(-T132);
  g1(8,8)=exp(y(8));
  g1(9,1)=(-T37);
  g1(9,4)=(-T94);
  g1(9,6)=(-(T35*T123));
  g1(9,7)=exp(y(7));
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],9,81);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],9,729);
end
end
end
end
