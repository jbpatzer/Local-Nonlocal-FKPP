function [zout, Uout] = solvespike(C,w,clearflag)
a = sqrt(pi)/w;
figure(1)
if clearflag, clf, end
hold on
la = fzero(@(la) rhs(la,C,w),log(a));
a = exp(la); 
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
f = w-integral(@(z) U(z,exp(la),C),-Inf,Inf);

%%
function f = U(z,a,C)
f = exp(-z.^2)./(a+C*erfc(z));