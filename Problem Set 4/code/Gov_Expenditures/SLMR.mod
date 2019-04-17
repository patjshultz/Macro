var a y i k c h r w lambda; 
varexo eps ;

parameters bbeta psi alpha delta gamma abar rho_a sigma_a; 

bbeta = 0.990000;
psi = 3.074648;
alpha = 0.666667;
delta = 0.025000;
gamma = 1.008000;
rho_a= 0.600000;
sigma_a = 0.006600;
abar = 0.000000;

initval;
lambda = 0.677635969756;
c = -0.677635969756;
w = 0.719991433117;
y = -0.301659814415;
h = -1.427116355640;
i = -1.461994449481;
k = 1.949253268034;
r = -3.349525371118;
a = 0.000000000000;
end;

model;
1/exp(c) = exp(lambda);
gamma = (bbeta/gamma) * (exp(c)/exp(c(+1))) *((1-alpha) * exp(a(+1)) * exp(k)^(-alpha) * exp(h(+1))^alpha + (1-delta));
psi/(1-exp(h)) = (1/exp(c)) * (alpha * exp(a) * exp(k)^(1-alpha) * exp(h)^(alpha-1));
exp(a)*k(-1)^(1-alpha)*exp(h)^(alpha) + (1-delta) * exp(k(-1)) - exp(k) * gamma = exp(c);
a = abar + rho_a * a(-1) + sigma_a * eps;
exp(y) = exp(a)*(exp(k(-1)))^(1-alpha) * exp(h)^(alpha);
exp(y) = exp(c) + exp(i);
exp(w) = alpha * exp(a) * exp(k(-1))^(1-alpha) * exp(h)^(alpha-1);
exp(r) = (1-alpha) * exp(a) * exp(k(-1))^(-alpha) * exp(h)^alpha;
end;

shocks;
var eps = 1; 
end;
steady;
options_.TeX = 1;
stoch_simul(order=1,graph,hp_filter=1600, irf = 60);
