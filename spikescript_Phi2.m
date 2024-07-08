alpha = 0.0065; fac = 0.99; fac1 = 0.95; alphamin = 5e-8;
y = []; z = 0; Nz = 1; alphaspike = []; count = 1;
zout = cell(1); Umax = cell(1); dzout = cell(1);
alpha1 = []; 
outcount = 1;
format compact
format short g
warning('off','MATLAB:nearlySingularMatrix')
while alpha>alphamin
    [minU,~,z,C] = Phi2asymp(Nz,alpha,z);
    if minU<0
        a0 = fsolve(@(alpha) Phi2asymp(Nz,alpha,z),alpha,optimset('display','off')); 
        alphaspike(count) = a0; 
        count = count+1;
        disp(alphaspike)
        alpha = a0;
        [~,minz,z,C] = Phi2asymp(Nz,alpha,z);
        alpha1 = [alpha1 alpha];
        zout{outcount} = [z minz];
        Umax{outcount} = (1-(1-alpha)*WB(0,C(1:8),alpha))/alpha;
        for k = 2:Nz
            Umax{outcount} = [Umax{outcount} (1-(1-alpha)*WB(z(2*k-1),C(8*k-7:8*k),alpha))/alpha];
        end
        Umax{outcount} = [Umax{outcount} 0];
        outcount = outcount+1;
        alpha = fac1*alpha;
        [~,~,~,C] = Phi2asymp(Nz,alpha,z);
        z3 = minz;
        while WB(z3,C(end-7:end),alpha)>1/(1-alpha)
            z3 = 0.999*z3;
        end
        z1 = fzero(@(z1) WB(z1,C(end-7:end),alpha)-1/(1-alpha),[z3/0.999 z3]);
        try
            z3 = minz;
            while (WB(z3,C(end-7:end),alpha)>0)&&(z3>minz-0.5)
                z3 = 1.01*z3;
            end
            z3 = z3/1.01;
            z2 = fzero(@(z2) z2-z1+sqrt(1-alpha)*integral(@(z) sqrt(WB(z,C(end-7:end),alpha)),z2,z1,'AbsTol',1e-14),[z3 z1-1e-8],optimset('display','off'));
        catch
            z3 = minz;
            while (WB(z3,C(end-7:end),alpha)<0)&&(z3>minz-0.5)
                z3 = 1.01*z3;
            end
            z3 = z3/1.01;
            z2 = fzero(@(z2) z2-z1+sqrt(1-alpha)*integral(@(z) sqrt(WB(z,C(end-7:end),alpha)),z2,z1,'AbsTol',1e-14),[z3 z1-1e-8],optimset('display','off'));
        end
        Nz = Nz+1;
        [~,~,z] = Phi2asymp(Nz,alpha,[z z1 z2]);
    else
        alpha1 = [alpha1 alpha]; %#ok<*AGROW>
        Umax{outcount} = (1-(1-alpha)*WB(0,C(1:8),alpha))/alpha;
        for k = 2:Nz
            Umax{outcount} = [Umax{outcount} (1-(1-alpha)*WB(z(2*k-1),C(8*k-7:8*k),alpha))/alpha];
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