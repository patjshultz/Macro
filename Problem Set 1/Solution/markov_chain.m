clear variables
clc

p=params_func();


N_simul=100000;
chain_pos=zeros(N_simul,1);

chain_pos(1)=ceil(p.Ny/2);

for i=2:N_simul
    chain_pos(i)=sum(rand >= cumsum([0, p.Py(chain_pos(i-1),:)]));
end

chain=log(p.y_grid(chain_pos));

disp(['MEAN    ','Theory=', num2str(p.mu_rho), ' Numerical=',num2str(mean(chain))])
disp(['STD     ','Theory=', num2str(sqrt(p.sigma_y^2/(1-p.rho_y^2))), ' Numerical=',num2str(std(chain))])
disp(['AC(1)   ','Theory=', num2str(p.rho_y), ' Numerical=',num2str(corr(chain(1:end-1),chain(2:end)))])
