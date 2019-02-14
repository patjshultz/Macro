function p=params_func()
    p=struct;
    
    %risk aversion
    p.gamma=2; 
    
    %time discounting
    p.discounting=0.95;
    
    %interest rate
    p.R=1.04;
    
    %income process
    p.rho_y=0.95;
    p.sigma_y=0.1;
    p.mu_rho=0.05/(1-p.rho_y);
    p.Ny=9;
    baseSigma=(0.5 + p.rho_y/4)*p.sigma_y +(0.5 - p.rho_y/4)*p.sigma_y/sqrt(1-p.rho_y^2);
    [p.y_grid,p.Py] = tauchenhussey(p.Ny,p.mu_rho,p.rho_y,p.sigma_y,baseSigma);
%     disp(p.y_grid)
    p.y_grid=exp(p.y_grid);
    
    %grid for capital
    p.K_min=0;
    p.K_max=100;
    p.Nk=1001;
    p.K_grid=linspace(p.K_min,p.K_max,p.Nk);
    
    %precision and maxeval
    p.precision=10^-7;
    p.maxeval=1500;
    

end

    
