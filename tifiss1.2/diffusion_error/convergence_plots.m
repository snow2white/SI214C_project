%CONVERGENCE_PLOTS Plotting convergence of energy error estimators
%    TIFISS scriptfile: DJS; 5 March 2017.
% Copyright (c) 2017 David Silvester, Alex Bespalov, Leonardo Rocchi

%% test problem 3: square domain, regular solution
%  coeff: a(x,y) = 1 + C*fx*fy
%    rhs: f(x,y) = (2-x^2-y^2)/8 
%
dof = [9, 49, 225, 961, 3969, 16129, 65025];
%err_hat = [1.0150e-02, 5.9369e-03, 3.1746e-03, 1.6395e-03, 8.3307e-04, 4.1993e-04, 2.1082e-04];
err_hat = [9.7550e-03, 5.8061e-03, 3.1187e-03, 1.6128e-03, 8.1986e-04, 4.1333e-04, 2.0752e-04];
err_bub = [1.2442e-02, 7.7281e-03, 4.3326e-03, 2.3015e-03, 1.1876e-03, 6.0354e-04, 3.0428e-04];

figure
loglog(dof,err_hat,'bo-') 
hold on
loglog(dof,err_bub,'r-sq')
hold on
loglog(dof,0.1*dof.^(-1/2),'k')  
xlabel('degree of freedom','Fontsize',17)
ylabel('estimated energy error','Fontsize',16)
legend('p/w linear bubbles','quadratic bubbles','N^{-1/2}')
title('test problem 3','Fontsize',18)
grid on
hold off



%% test problem 5: L-shaped domain, singular solution
%  singular solution: u(r,theta) = r^(2/3)sin((pi + 2theta)/2) with
%  coeff: a(x,y) = 1 - (x^2 + y^2)/4
%    rhs: f(x,y) = u(r,theta)/3
%
dof = [5, 33, 161, 705, 2945, 12033, 48641];
%err_hat = [3.2083e-01, 2.1347e-01, 1.3377e-01, 8.2106e-02, 5.0324e-02, 3.1016e-02, 1.9237e-02];
 err_hat = [3.3078e-01, 2.1279e-01, 1.3229e-01, 8.1236e-02, 4.9914e-02, 3.0839e-02, 1.9163e-02];
err_bub = [3.4129e-01, 2.4299e-01, 1.6149e-01, 1.0422e-01, 6.6436e-02, 4.2118e-02, 2.6630e-02];

figure
loglog(dof,err_hat,'bo-')  
hold on
loglog(dof,err_bub,'r-sq')
hold on
loglog(dof,0.85*dof.^(-1/3),'k')
xlabel('degree of freedom','Fontsize',17)
ylabel('estimated energy error','Fontsize',16)
legend('p/w linear bubbles','quadratic bubbles','N^{-1/3}')
title('test problem 5','Fontsize',18)
grid on
hold off


