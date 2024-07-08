function [minU,minz,z,zg,U,zW,ppW] = asymp_inteqn(Nz,alpha,z0,ker,zmin)
if isempty(z0)
    disp('Need a guess')
    minU = []; minz = []; z = []; zg = []; zW=[]; ppW=[];
    return
end

options = optimoptions('fsolve','display','off','MaxFunctionEvaluations',1e6,'MaxIterations',1e6,'OptimalityTolerance',1e-10,'FunctionTolerance',1e-12);

if Nz==1
    [U,zg,zW,ppW] = solvelinear_inteqn(alpha,0,ker,zmin);
    % C = solvelinear(Nz,alpha,0);
    z = 0;
else
    dz0 = -diff(z0);
    logdz = fsolve(@(logdz) rhs(logdz,alpha,ker,zmin),log(dz0),options);
    dz = exp(logdz);
    z = [0 -cumsum(dz)];
    z = real(z);
    [U,zg,zW,ppW] = solvelinear_inteqn(alpha,z,ker,zmin);
end
W = (1-alpha*U{1})/(1-alpha); zW1 = zg{1};
Uplot = [U{1} 0]; zU = [zg{1} zg{1}(end)];
for k = 1:Nz-1
    W = [W ppval(ppW{k},zW{k}) (1-alpha*U{k+1})/(1-alpha)];
    zW1 = [zW1 zW{k} zg{k+1}];
    Uplot = [Uplot U{k+1} 0]; %#ok<*AGROW>
    zU = [zU zg{k+1} zg{k+1}(end)];
end
figure(3),clf
% plot(zU,Uplot,'k',zW1,1+100*(W-1),'b','linewidth',1.4)
plot(zU,Uplot,'k','linewidth',1.4)
grid on
xlabel('z'), ylabel('U, 1+100(W-1)')
title(sprintf('\\alpha = %1.3e',alpha))
ylim([0 5])
drawnow
[~,k] = min(U{1}); pp0 = spline(zg{1},U{1});
[minz,minU] = fminbnd(@(z) ppval(pp0,z),zU(max(1,k-1)),zU(min(k+1,length(zU))));

%%
function f = rhs(logdz,alpha,ker,zmin)
dz = exp(logdz);
z = [0 -cumsum(dz)];
% C = solvelinear(Nz,alpha,z);
[U,~,~,ppW] = solvelinear_inteqn(alpha,z,ker,zmin);
Nz = length(U);
f = zeros(2*Nz-2,1);
for k = 1:Nz-1
    % f(2*k-1) = WB(z(2*k),C(8*k-7:8*kR),alpha)-1/(1-alpha);
    f(2*k-1) = U{k+1}(1);%ppval(ppW{k},z(2*k))-1/(1-alpha);
    % f(2*k) = z(2*k+1)-z(2*k)+sqrt(1-alpha)*integral(@(z) sqrt(WA(z,C(8*k+1:8*k+8))),z(2*k+1),z(2*k),'AbsTol',1e-14);
    f(2*k) = z(2*k+1)-z(2*k)+sqrt(1-alpha)*integral(@(z) sqrt(ppval(ppW{Nz-k},z)),z(2*k+1),z(2*k),'AbsTol',1e-14);
end