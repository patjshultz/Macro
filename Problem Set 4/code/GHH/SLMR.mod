var a y i k c h r w lambda; 
varexo eps ;

parameters bbeta psi alpha delta gamma abar rho_a sigma_a nu; 

bbeta = 0.990000;
psi = 4.193558;
alpha = 0.666667;
delta = 0.025000;
gamma = 1.008000;
rho_a= 0.600000;
sigma_a = 0.006600;
abar = 0.000000;
nu = 0.500000;

initval;
lambda = 1.719757541289;
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
1/(exp(c) - psi *(exp(h)^(1+nu)/(1+nu))) = exp(lambda);
gamma = (bbeta/gamma) * (exp(lambda(+1))/exp(lambda)) *((1-alpha) * exp(a(+1)) * exp(k)^(-alpha) * exp(h)^alpha + (1-delta));
psi * exp(h)^nu = alpha * exp(a) * exp(k(-1))^(1-alpha)*exp(h)^(alpha-1);
exp(a)*k(-1)^(1-alpha)*exp(h)^(alpha) + (1-delta) *exp(k(-1)) - exp(k) *gamma = exp(c);
a = rho_a * a(-1) + sigma_a * eps;
exp(y) = exp(a)*(exp(k(-1)))^(alpha) * exp(h)^(1-alpha);
exp(y) = exp(c)+exp(i);
exp(w) = alpha * exp(a) * exp(k(-1))^(1-alpha) * exp(h)^(alpha-1);
exp(r) = (1-alpha) * exp(a) * exp(k(-1))^(-alpha)*exp(h)^alpha;
end;

shocks;
var eps = 1; 
end;
steady;
stoch_simul(order=1,graph,hp_filter=1600);
