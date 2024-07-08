function [minU,minz,z,C,zplot,U] = Phi2asymp(Nz,alpha,z0)
%% Solve for the leading order asymptotic solution with c<<1 and alpha = O(1) for phi = Phi_1
if isempty(z0)
    disp('Need a guess')
    minU = []; minz = []; z = []; C = [];
    return
end

options = optimoptions('fsolve','display','off','MaxFunctionEvaluations',1e6,'MaxIterations',1e6,'OptimalityTolerance',1e-10,'FunctionTolerance',1e-12);

if Nz==1
    C = solvelinear(Nz,alpha,0);
    z = 0;
else
    dz0 = -diff(z0);
    logdz = fsolve(@(logdz) rhs(logdz,Nz,alpha),log(dz0),options);
    dz = exp(logdz);
    z = [0 -cumsum(dz)];
    z = real(z);
    C = solvelinear(Nz,alpha,z);
end
dz = 0.01;
zplot = linspace(5, 0,ceil(5/dz));
U = zeros(size(zplot));
W = WA(zplot,C(1:8));
for k = 1:Nz-1
    znew = linspace(zplot(end), z(2*k), ceil(zplot(end)-z(2*k))/dz);
    W0 = WB(znew,C(8*k-7:8*k),alpha);
    W = [W W0];
    U = [U (1-(1-alpha)*W0)/alpha];
    zplot = [zplot znew]; %#ok<*AGROW>
    znew = linspace(zplot(end), z(2*k+1), ceil(zplot(end)-z(2*k+1))/dz);
    W0 = WA(znew,C(8*k+1:8*k+8));
    W = [W W0];
    U = [U zeros(size(znew))];
    zplot = [zplot znew];
end
znew = zplot(end):-dz:-5;
W0 = WB(znew,C(end-7:end),alpha);
W = [W W0];
U = [U (1-(1-alpha)*W0)/alpha];
zplot = [zplot znew];
figure(1),clf
plot(zplot,U,'k',zplot,1+100*(W-1),'b','linewidth',1.4)
grid on
xlabel('z'), ylabel('U, 1+100(W-1)')
title(sprintf('\\alpha = %1.3e',alpha))
axis([zplot(end) 0 0 5])
drawnow
U0 = fliplr((1-(1-alpha)*W0)/alpha);
z0 = fliplr(znew);
[~,k] = min(U0); pp0 = spline(z0,U0);
[minz,minU] = fminbnd(@(z) ppval(pp0,z),z0(max(1,k-1)),z0(min(k+1,length(z0))));

%%
function f = rhs(logdz,Nz,alpha)
dz = exp(logdz);
z = [0 -cumsum(dz)];
C = solvelinear(Nz,alpha,z);
f = zeros(2*Nz-2,1);
for k = 1:Nz-1
    f(2*k-1) = WB(z(2*k),C(8*k-7:8*k),alpha)-1/(1-alpha);
    f(2*k) = z(2*k+1)-z(2*k)+sqrt(1-alpha)*integral(@(z) sqrt(WA(z,C(8*k+1:8*k+8))),z(2*k+1),z(2*k),'AbsTol',1e-14);
end