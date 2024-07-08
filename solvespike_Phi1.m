function [sol,w,umax] = solvespike_Phi1(c0,clearflag)
z1 = -10; z2 = -z1;
sol = ode45(@(z,y) rhs(z,y,c0),[z1 z2],[z1 1 0 0],odeset('AbsTol',1e-14,'RelTol',5e-14));
z = linspace(z1,z2,1000);
figure(1)
if clearflag,clf, end
L = deval(sol,z,1);
w = deval(sol,z2,4);
pp = spline(z,L);
[~,Lmax] = fminbnd(@(z) -ppval(pp,z),z1,z2);
umax = exp(-Lmax);
% subplot(1,2,1)
% plot(z,L,'linewidth',1.4)
% grid on
% xlabel('$\bar{z}$','interpreter','latex')
% ylabel('$\bar{L}$','interpreter','latex')
% title(sprintf('c_0-1 = %1.2e, w = %1.2e, u_{max} = %2.2f',c0-1,w,umax))
% hold on
% subplot(1,2,2)
plot(z,exp(L),'linewidth',1.4)
grid on
xlabel('$\bar{z}$','interpreter','latex')
ylabel('$U$','interpreter','latex')
hold on
drawnow

%%
function f = rhs(~,y,c0)
f = [ y(2); y(3); exp(y(1))*(-1+(y(3)+y(2)^2)/c0); exp(y(1))];