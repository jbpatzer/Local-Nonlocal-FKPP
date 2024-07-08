function [zout, Uout] = solvespike(C,w)
%% Find the solution of (3.4) that satisfies (3.6)
% C and w the parameters in (3.4) and (3.6)

a = sqrt(pi)/w; % initial guess of a

la = fzero(@(la) rhs(la,C,w),log(a));
a = exp(la); 

%Plotting
z1 = -1; while U(z1,a,C)>1e-6,z1 = z1-1; end
z2 = 1; while U(z2,a,C)>1e-6,z2 = z2+1; end
zout = linspace(z1,z2,5000);
Uout = U(zout,a,C);
pp = spline(zout,Uout);
zout = zout-fminbnd(@(z) -ppval(pp,z),z1,z2);
plot(zout,Uout,'linewidth',1.4)
grid on
xlabel('$\hat{z}$','interpreter','latex')
ylabel('$\bar{U}$','interpreter','latex')
title(sprintf('C = %2.2e, w = %2.2e, a = %2.2e',C,w,a))
drawnow
%%
function f = rhs(la,C,w)
% Equation (3.6)
f = w-integral(@(z) U(z,exp(la),C),-Inf,Inf);

%%
function f = U(z,a,C)
% The solution (3.4)
f = exp(-z.^2)./(a+C*erfc(z));