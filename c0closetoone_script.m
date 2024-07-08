z1 = -100; z2 = 10;
sol = ode45(@(z,y) [y(2); y(3); exp(z)*(-1+2*y(2)+y(3))],[z1 z2],[exp(z1) exp(z1) exp(z1)],...
    odeset('AbsTol',1e-14,'RelTol',1e-14));
z = linspace(z1,z2,2000);
ep = 10.^(-1:-1:-100);
M = length(ep);
umax = zeros(1,M); w = zeros(1,M);
for k = 1:M
    figure(1),clf
    L = z+ep(k)*deval(sol,z,1);
    L(isnan(L)) = -Inf;
    u = exp(L);
    plot(z,u,log(log(1/ep(k)))*[1 1],[0 max(u)],'k--','linewidth',1.4)
    grid on
    xlabel('z'), ylabel('U')
    drawnow
    pp = spline(z,u);
    w(k) = integral(@(z) ppval(pp,z),z1,z2);
    [~,mu] = fminbnd(@(z) -ppval(pp,z),z1,z2);
    umax(k) = -mu;
end
figure(2), clf
plot(z,-deval(sol,z,1)./exp(exp(z)),'linewidth',1.4)
grid on
ylim([0 1])
