function [sol,c0,umax] = solvespike_Phi1(c0hat)
z1 = -10; z2 = -z1;
sol = ode45(@(z,y) rhs(z,y,c0hat),[z1 z2],[z1 1 0 0],odeset('AbsTol',1e-14,'RelTol',5e-14));

% Plotting
z = linspace(z1,z2,1000);
figure(1)
L = deval(sol,z,1);
w = deval(sol,z2,4);
c0 = c0hat/w;
pp = spline(z,L);
[~,Lmax] = fminbnd(@(z) -ppval(pp,z),z1,z2);
umax = exp(-Lmax);
plot(z,exp(L),'linewidth',1.4)
grid on
xlabel('$\bar{z}$','interpreter','latex')
ylabel('$U$','interpreter','latex')
hold on
drawnow

%%
function f = rhs(~,y,c0hat)
f = [ y(2); y(3); exp(y(1))*(-1+(y(3)+y(2)^2)/c0hat); exp(y(1))];