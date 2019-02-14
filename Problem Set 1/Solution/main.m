clear variables
clc

p=params_func();

%% VFI
policy.K_prime=zeros(p.Nk,p.Ny);
Pos=repmat(reshape(1:p.Nk,p.Nk,1),1,p.Ny);


Y_full=repmat(reshape(p.y_grid,1,p.Ny),p.Nk,1,p.Nk);
K_prime_full=repmat(reshape(p.K_grid,1,1,p.Nk),p.Nk,p.Ny,1);
K_full=repmat(reshape(p.K_grid,p.Nk,1),1,p.Ny,p.Nk);
C_full=Y_full+p.R*K_full-K_prime_full;
Util=C_full.^(1-p.gamma)/(1-p.gamma);
Util(C_full<0)=-10^10;

V=zeros(p.Nk,p.Ny);

tic
for l=1:p.maxeval
    V_prev=V;
    EV=V*p.Py';
    V_new=Util+p.discounting*repmat(reshape(EV',1,p.Ny,p.Nk),p.Nk,1,1);
% %     V_new=Util+p.discounting*reshape(EV',1,p.Ny,p.Nk);
    [V,Pos]=max(V_new,[],3);

   
    diff=max(max(abs(V-V_prev)));
%     disp([l,diff/p.precision])
    if diff<p.precision
        disp('Convergence achieved')
        break
    end
    if l==p.maxeval
        disp('MaxEval achieved, no convergence')
    end
end
toc
policy.V=V;
policy.K_prime=p.K_grid(Pos);
mat=reshape(1:p.Nk*p.Ny,p.Nk,p.Ny);
policy.C=C_full((Pos-1)*p.Nk*p.Ny+mat);


save policy policy


% %%
%graphs
figure(1)
subplot(2,2,1)
plot(p.k_grid,V(:,1),p.k_grid,V(:,3),p.k_grid,V(:,5),'Linewidth',2)
legend({'Low Z','Medium Z', 'High Z'},'Fontsize',13,'Location','northwest')
title('V(k)','Fontsize',13)
subplot(2,2,2)
plot(k_grid,k_grid,'black')
hold on
plot(k_grid,Kap(:,1),k_grid,Kap(:,3),k_grid,Kap(:,5),'Linewidth',2)
hold off
title('k^\prime(k)','Fontsize',13)
subplot(2,2,3)
plot(k_grid,Inv(:,1),k_grid,Inv(:,3),k_grid,Inv(:,5),'Linewidth',2)
title('i(k)','Fontsize',13)
subplot(2,2,4)
plot(k_grid,Div(:,1),k_grid,Div(:,3),k_grid,Div(:,5),'Linewidth',2)
title('d(k)','Fontsize',13)
