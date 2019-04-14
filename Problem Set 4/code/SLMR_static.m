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

T10 = 1/exp(y(5));
T16 = params(1)/params(5);
T26 = (1-params(3))*exp(y(1))*exp(y(4))^(-params(3));
T29 = exp(y(6))^params(3);
T30 = T26*T29;
T41 = params(3)*exp(y(1))*exp(y(4))^(1-params(3));
T43 = exp(y(6))^(params(3)-1);
T44 = T41*T43;
T48 = exp(y(1))*y(4)^(1-params(3));
T49 = T29*T48;
T65 = exp(y(1))*exp(y(4))^params(3);
T66 = exp(y(6))^(1-params(3));
T67 = T65*T66;
lhs =T10;
rhs =exp(y(9));
residual(1)= lhs-rhs;
lhs =params(5);
rhs =T16*(T30+1-params(4));
residual(2)= lhs-rhs;
lhs =params(2)/(1-exp(y(6)));
rhs =T10*T44;
residual(3)= lhs-rhs;
lhs =T49+exp(y(4))*(1-params(4))-params(5)*exp(y(4));
rhs =exp(y(5));
residual(4)= lhs-rhs;
lhs =y(1);
rhs =y(1)*params(7)+params(8)*x(1);
residual(5)= lhs-rhs;
lhs =exp(y(2));
rhs =T67;
residual(6)= lhs-rhs;
lhs =exp(y(2));
rhs =exp(y(5))+exp(y(3));
residual(7)= lhs-rhs;
lhs =exp(y(8));
rhs =T44;
residual(8)= lhs-rhs;
lhs =exp(y(7));
rhs =T30;
residual(9)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(9, 9);

  %
  % Jacobian matrix
  %

T90 = T29*(1-params(3))*exp(y(1))*exp(y(4))*getPowerDeriv(exp(y(4)),(-params(3)),1);
T96 = T43*params(3)*exp(y(1))*exp(y(4))*getPowerDeriv(exp(y(4)),1-params(3),1);
T117 = exp(y(6))*getPowerDeriv(exp(y(6)),params(3),1);
T128 = T41*exp(y(6))*getPowerDeriv(exp(y(6)),params(3)-1,1);
  g1(1,5)=(-exp(y(5)))/(exp(y(5))*exp(y(5)));
  g1(1,9)=(-exp(y(9)));
  g1(2,1)=(-(T16*T30));
  g1(2,4)=(-(T16*T90));
  g1(2,6)=(-(T16*T26*T117));
  g1(3,1)=(-(T10*T44));
  g1(3,4)=(-(T10*T96));
  g1(3,5)=(-(T44*(-exp(y(5)))/(exp(y(5))*exp(y(5)))));
  g1(3,6)=(-(params(2)*(-exp(y(6)))))/((1-exp(y(6)))*(1-exp(y(6))))-T10*T128;
  g1(4,1)=T49;
  g1(4,4)=exp(y(4))*(1-params(4))+T29*exp(y(1))*getPowerDeriv(y(4),1-params(3),1)-params(5)*exp(y(4));
  g1(4,5)=(-exp(y(5)));
  g1(4,6)=T48*T117;
  g1(5,1)=1-params(7);
  g1(6,1)=(-T67);
  g1(6,2)=exp(y(2));
  g1(6,4)=(-(T66*exp(y(1))*exp(y(4))*getPowerDeriv(exp(y(4)),params(3),1)));
  g1(6,6)=(-(T65*exp(y(6))*getPowerDeriv(exp(y(6)),1-params(3),1)));
  g1(7,2)=exp(y(2));
  g1(7,3)=(-exp(y(3)));
  g1(7,5)=(-exp(y(5)));
  g1(8,1)=(-T44);
  g1(8,4)=(-T96);
  g1(8,6)=(-T128);
  g1(8,8)=exp(y(8));
  g1(9,1)=(-T30);
  g1(9,4)=(-T90);
  g1(9,6)=(-(T26*T117));
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
