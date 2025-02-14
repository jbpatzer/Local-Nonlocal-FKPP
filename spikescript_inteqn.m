%% Find the leading order asymptotic solution of the TW problem with c << 1 and alpha = O(1) with phi = Phi_infty
% Performs adaptive continuation in alpha
alpha = 0.025; fac = 0.99; fac1 = 0.95; alphamin = 1e-10;
y = []; z = 0; Nz = 1; alphaspike = []; count = 1; ker = Inf; zmin = 20;
zout = cell(1); Umax = cell(1); dzout = cell(1);
alpha1 = []; 
outcount = 1;
format compact
format short g
warning('off','MATLAB:nearlySingularMatrix')
while alpha>alphamin
    save dum Nz alpha z ker zmin
    [minU,~,z,zg,U,zW,ppW] = asymp_inteqn(Nz,alpha,z,ker,zmin);
    if minU<0
        a0 = fsolve(@(alpha) asymp_inteqn(Nz,alpha,z,ker,zmin),alpha,optimset('display','off')); 
        alphaspike(count) = a0; 
        count = count+1;
        disp(alphaspike)
        alpha = a0;
        [~,minz,z,~,U] = asymp_inteqn(Nz,alpha,z,ker,zmin);
        alpha1 = [alpha1 alpha];
        zout{outcount} = [z minz];
        Umax{outcount} = [];
        for k = 1:Nz
            Umax{outcount} = [Umax{outcount} max(U{k})];
        end
        outcount = outcount+1;
        alpha = fac1*alpha;
        [~,~,~,zg,U,zW,ppW] = asymp_inteqn(Nz,alpha,z,ker,zmin);
        z3 = minz;
        ppU = spline(zg{1},U{1});
        while ppval(ppU,z3)<0
            z3 = 0.999*z3;
        end
        z1 = fzero(@(z1) ppval(ppU,z1),[z3/0.999 z3]);
    try
        z3 = minz;
        while (1-alpha*ppval(ppU,z3)>0)&&(z3>minz-0.5)
            z3 = 1.01*z3;
        end
        z3 = z3/1.01;
        z2 = fzero(@(z2) z2-z1+integral(@(z) sqrt((1-alpha*ppval(ppU,z))),z2,z1,'AbsTol',1e-14),[z3 z1-1e-8],optimset('display','off'));
    catch
        z3 = minz;
        while (1-alpha*ppval(ppU,z3)<0)&&(z3>minz-0.5)
            z3 = 1.01*z3;
        end
        z3 = z3/1.01;
        z2 = fzero(@(z2) z2-z1+integral(@(z) sqrt((1-alpha*ppval(ppU,z))),z2,z1,'AbsTol',1e-14)/sqrt(1-alpha),[z3 z1-1e-8],optimset('display','off'));
    end
        Nz = Nz+1;
        [~,~,z] = asymp_inteqn(Nz,alpha,[z z1 z2],ker,zmin);
    else
        alpha1 = [alpha1 alpha]; %#ok<*AGROW>
        Umax{outcount} = [];
        for k = 1:Nz
            Umax{outcount} = [Umax{outcount} max(U{k})];
        end
        zout{outcount} = z; %#ok<*SAGROW>
        outcount = outcount+1;
        alpha = alpha*fac;
    end
end

return
%%
k = 1;
while all(isreal(zout{k}))&&(k<length(zout))
    k = k+1;
end
zout = zout(1:k);

figure(2), clf
N0 = length(zout{end});
for N = N0:-1:2
    k = length(zout);
    zplot = zout{k}(N);
    while (k>1)&&(length(zout{k-1})>=N)
        k = k-1;
        zplot = [zplot zout{k}(N)];
    end
    if mod(N,2)==1
        zplot = [zplot zout{k-1}(N-1)];
    else
        k = k+1;
    end
    semilogx(alpha1(k-1:end),fliplr(zplot),'k','linewidth',1.4)
    hold on
end
semilogx(alpha1([1 end]),[0 0],'k','linewidth',1.4)
xlabel('\alpha'), ylabel('z')
grid on
axis tight
hold off

%%
figure(3), clf
N0 = length(Umax{end});
for N = N0:-1:1
    k = length(Umax);
    Uplot = Umax{k}(N);
    while (k>1)&&(length(Umax{k-1})>=N)
        k = k-1;
        Uplot = [Uplot Umax{k}(N)];
    end
    if N>1
        Uplot = [Uplot 0];
    else
        k = k+1;
    end
    loglog(alpha1(k-1:end),fliplr(Uplot),'k','linewidth',1.4)
    hold on
end
loglog(alpha1, 1./sqrt(alpha1),'b--','linewidth',1.4)
xlabel('\alpha'), ylabel('U_{max}')
grid on
axis tight
hold off

%%
figure(4), clf
for k = 1:length(zout)
    dzout{k}  = -diff(zout{k});
end
N0 = length(dzout{end});
for N = N0:-1:1
    k = length(dzout);
    dzplot = dzout{k}(N);
    while (k>1)&&(length(dzout{k-1})>=N)
        k = k-1;
        dzplot = [dzplot dzout{k}(N)];
    end
    k = k+1;
    loglog(alpha1(k-1:end),fliplr(dzplot),'k','linewidth',1.4)
    hold on
end
loglog(alpha1, sqrt(alpha1),'b--','linewidth',1.4)
xlabel('\alpha'), ylabel('\Deltaz')
grid on
axis tight
hold off

save spikedata alpha1 zout Umax alphaspike dzout