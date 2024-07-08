function [z,L,f,fval,w] = TW(z0,L0,zinit,ker,c,D0,alpha,minflag)
% z0, L0 = initial guess of solution. Take z0 = L0 = [] if not available.
% zinit = initial, coarse grid
% ker = kernel Phi_ker
% c = wavespeed
% diffusion parameter is D = D0*c^2
% alpha = nonlinearity parameter
% calculate local minimum behind trough minflag - used to calculate figs 10
% and 12

persistent uint

format compact
g = 0.5*(1-1/sqrt(3)); %Gaussian quadrature points
D = D0*c^2;

dzinit = diff(zinit);
dz0pp = spline(zinit,[dzinit(1) 0.5*(dzinit(1:end-1)+dzinit(2:end)) dzinit(end)]);
% spline of initial grid spacing

L1 = zinit(1); L2 = zinit(end); %end points of grid

if D==0
    kplus = -1/c;
else
    kplus = (-c+sqrt(c^2-4*D))/2/D;
end
% L ~ kplus*z as z -> infinity

z = z0; L = L0;
if isempty(z)
    z = zinit;
end
dz = diff(z);
zmid = 0.5*(z(1:end-1)+z(2:end));
N = length(z); 

if isempty(L)
    %Heaviside initial guess
    u = ones(1,N); u(z>0) = 1e-12;
    L = log(u);
    uint = sum(0.5*dz.*(u(1:N-1)+u(2:N)));
end

if isempty(uint)
    u = exp(L); dz = diff(z);
    uint = sum(0.5*dz.*(u(1:end-1)+u(2:end)));
end

% Fix wavefront by specifying area under u(z).

tol = 1e-3; tol0 = 0.1;
if isempty(z0)
    Nmax = 2*length(zinit);
else
    Nmax = 2*length(z0);
end

Nnew = 0;
options = optimoptions('fsolve','SpecifyObjectiveGradient',true,'display','off',...
    'MaxIter',100,'FunctionTolerance',tol^2,'StepTolerance',1e-16);

[F,M1,M2,dmat2] = setup(z,zmid,dz,N,alpha,g,ker); % Set up constant matrices
[L,fval] = fsolve(@(L) rhs(L,F,M1,M2,c,D,alpha,dz,N,uint,dmat2,kplus),L,options); %Solve discretised equations
nfv = norm(fval);
disp(nfv)

regridflag = true; regridcount = 0; regridcountmax = 8; %counters for regridding
maxposold = []; frontposold = []; %arrays for finding local maxima and fronts

%Regridding loop
while (nfv<1e12)&&(regridcount<regridcountmax)&&(regridflag||(nfv>tol0))&&(Nnew<=Nmax)
    u = exp(L);
    %Find spikes
    [~,maxpos] = findpeaks(u,z,'MinPeakProminence',0.5);

    %Find fronts
    dudz = diff(u)./diff(z); zmid = 0.5*(z(1:end-1)+z(2:end));
    [~,frontpos] = findpeaks(-dudz,zmid,'MinPeakHeight',10);

    dudzpp = spline(zmid,dudz);

    %Regrid
    [znew,Lnew,Nnew,dznew] = setup_grid(maxpos,frontpos,dudzpp,L1,L2,z,L,dz0pp);

    if (length(maxpos)~=length(maxposold))||(length(frontpos)~=length(frontposold))
        shiftflag = true;
    else
        dpos = [abs(maxpos-maxposold) abs(frontpos-frontposold)];
        shiftflag = any(dpos>10*dznew);
    end
    frontposold = frontpos; maxposold = maxpos;

    if ((shiftflag)||(((abs(Nnew-N)>5)&&(nfv>tol^2))||((regridcount==0)&&(N~=Nnew))||(nfv>tol)))&&(Nnew<=Nmax)
        %Solve on new grid
        disp(['Regridding, old N = ' num2str(N), ', new N = ' num2str(Nnew)])

        L = Lnew; z = znew; N = Nnew;
        dz = diff(z);
        zmid = 0.5*(z(1:end-1)+z(2:end));

        [F,M1,M2,dmat2] = setup(z,zmid,dz,Nnew,alpha,g,ker);
        [L, fval] = fsolve(@(L) rhs(L,F,M1,M2,c,D,alpha,dz,N,uint,dmat2,kplus),L,options);
        nfv = norm(fval);
        disp(nfv)
        regridcount = regridcount+1;
    else
        regridflag = false;
    end
end

u = exp(L);
w = -F/(1-alpha) + u(1:N-1)*M1'+u(2:N)*M2';

%Plotting
if (nfv<tol0)&&(Nnew<=Nmax)
    figure(1), clf
    grid on
    xlabel('z'), ylabel('W')
    subplot(5,1,3)
    plot(z,u,'k.-')
    axis([L1 L2 0 1.1*max(u)])
    grid on
    xlabel('z'), ylabel('U')
    subplot(5,1,4)
    semilogy(z(1:end-1),diff(z),'k','linewidth',1.4)
    xlim([L1 L2])
    grid on
    xlabel('z'), ylabel('dz')
    subplot(5,1,2)
    plot(z,L,'k','linewidth',1.4)
    p = findpeaks(-L(1:end-10));
    axis([L1 L2 -max(p) 2])
    ylabel('L')
    grid on
    subplot(5,1,1)
    plot(z,u,'k','linewidth',1.4)
    axis([L1 L2 0 2])
    title(['D_0 = ' num2str(D0) ', c = ' num2str(c) ', \alpha = ' num2str(alpha) ', N = ' num2str(N)])
    grid on
    ylabel('U')
    subplot(5,1,5)
    if isempty(frontpos)
        Nbehind = length(zmid);
    else
        Nbehind = length(zmid(zmid<frontpos(end)));
    end
    axis([L1 L2 1+1.1*max(abs(w(1:Nbehind)-1))*[-1 1]])
    plot(zmid,w-1,'b','linewidth',1.4)
    axis([L1 L2 1.1*min(w-1) 1.1*max(w-1)])
    grid on
    ylabel('w-1')
    drawnow
end

if ~isempty(minflag)
    f0 = 1e-2; %f = 0 when local minimum is f0. Used to calculate figs 10 and 12
    p = fliplr(-findpeaks(-L(z<0),'MinPeakHeight',0));
    if length(p)>=minflag
        f = exp(p(minflag))-f0;
    else
        f = 1-f0;
    end
else
    f = [];
end

if N>=Nmax, fval = 1e12; end

%%
function [f,J] = rhs(L,F,M1,M2,c,D,alpha,dz,N,uint,dmat2,kplus)
f = zeros(1,N); u = exp(L);
Lz = (L(2:N)-L(1:N-1))./dz;
Lzz = L*dmat2';
f(1:N-1) = c*Lz+D*(Lzz(1:N-1)+Lz.^2)+F+...
    1-(1-alpha)*u(1:N-1)*M1'-(1-alpha)*u(2:N)*M2'-0.5*alpha*(u(1:N-1)+u(2:N));
f(N-1) = (Lz(N-1)-kplus); %Far field condition
f(N) = 0.5*sum((u(1:N-1)+u(2:N)).*dz)-uint; %Area condition

% Construct Jacobian
J = zeros(N);
J(1:N-2,1:N) = D*dmat2(1:N-2,1:N);
dum = (c+2*D*Lz)./dz;
J(1:N-1,1:N-1) = J(1:N-1,1:N-1)-(1-alpha)*repmat(u(1:N-1),N-1,1).*M1-spdiags(dum',0,N-1,N-1);
J(1:N-1,2:N) = J(1:N-1,2:N)-(1-alpha)*repmat(u(2:N),N-1,1).*M2+spdiags(dum',0,N-1,N-1);
for k = 1:N-2
    J(k,k) = J(k,k)-0.5*alpha*u(k);
    J(k,k+1) = J(k,k+1)-0.5*alpha*u(k+1);
end
J(N-1,:) = [zeros(1,N-2) -1 1]/dz(N-1);
J(N,:) = 0.5*[dz(1)*u(1) (dz(2:N-1)+dz(1:N-2)).*u(2:N-1) dz(N-1)*u(N)];

%%
function [F,M1,M2,dmat2] = setup(z,zmid,dz,N,alpha,g,ker)
F = -(1-alpha)*(0.5-kernel(zmid-z(1),ker,true)); %Convolution phi*(1-u) for z>L2
z1 = z(1:N-1)+g*dz; z2 = z(1:N-1)+(1-g)*dz;

%Ml and Mr are the matrices that give the convolution
% at the quadrature points.
Mleft = zeros(N-1); Mright = zeros(N-1);
for i = 1:N-1
    Mleft(i,:) = 0.5*dz.*kernel(zmid(i)-z1,ker,false);
    Mright(i,:) = 0.5*dz.*kernel(zmid(i)-z2,ker,false);
end
[~,dmat2] = setupdmat(dz);
M1 = (1-g)*Mleft+g*Mright;
M2 = g*Mleft+(1-g)*Mright;

%%
function [znew,Lnew,Nnew,dzpos] = setup_grid(zmax,zdumax,dudzpp,L1,L2,z,L,dz0pp)
%Regrids based on value of u at the spikes and their width, and slope of fronts
N = length(L);
if isempty(zmax)&&isempty(zdumax)
    znew = z; Lnew = L; Nnew = N;
else
    u = exp(L); ppu = spline(z,u);
    Nz1 = length(zmax); Nz2 = length(zdumax); Nz = Nz1+Nz2;
    dznew = zeros(1,Nz);
    dzpp = cell(1,Nz);
    fac = 1.04; Nfac = 8;
    active = false(1,Nz);
    for i = 1:Nz1
        % For each spike make a spline for a new grid if spike is narrow enough
        dz0 = ppval(dz0pp,zmax(i));
        umax = ppval(ppu,zmax(i));
        zp0 = zmax(i); while ppval(ppu,zp0)>0.5*umax, zp0 = zp0+dz0; end
        zp = fzero(@(z) ppval(ppu,z)-0.5*umax,[zp0-dz0 zp0]);
        dznew(i) = (zp-zmax(i))/(Nfac+1);
        if dznew(i)<ppval(dz0pp,zmax(i))
            active(i) = true;
            znew = zmax(i)+linspace(-Nfac*dznew(i),Nfac*dznew(i),2*Nfac+1);
            if length(znew)>30000
                Nnew = 1e12; znew = []; Lnew = []; dzpos = zeros(1,Nz);
                return
            end
            dz = dznew(i);
            while znew(end)<L2
                znew = [znew znew(end)+dz]; dz = min(ppval(dz0pp,znew(end)),fac*dz);
            end
            dz = dznew(i);
            while znew(1)>L1
                znew = [znew(1)-dz znew]; dz = min(ppval(dz0pp,znew(1)),fac*dz); %#ok<*AGROW>
            end
            dz1 = diff(znew);
            dz1 = [dz1(1) 0.5*(dz1(1:end-1)+dz1(2:end)) dz1(end)];
            dzpp{i} = spline(znew,smooth(znew,dz1,11));
        else
            dzpp{i} = dz0pp;
        end
    end

    Nfac2 = 8;
    for i = 1:Nz2
        % For each front make a spline for a new grid if front is narrow enough
        dumax = ppval(dudzpp,zdumax(i)); umax = ppval(ppu,zdumax(i));
        dznew(Nz1+i) = abs(umax/dumax)/(Nfac2+1);
        if dznew(Nz1+i)<ppval(dz0pp,zdumax(i))
            active(Nz1+i) = true;
            znew = zdumax(i)+linspace(-2*Nfac2*dznew(Nz1+i),2*Nfac2*dznew(Nz1+i),4*Nfac2+1);
            if length(znew)>30000
                Nnew = 1e12; znew = []; Lnew = []; dzpos = zeros(1,Nz);
                return
            end
            dz = dznew(Nz1+i);
            while znew(end)<L2
                znew = [znew znew(end)+dz]; dz = min(ppval(dz0pp,znew(end)),fac*dz);
            end
            dz = dznew(Nz1+i);
            while znew(1)>L1
                znew = [znew(1)-dz znew]; dz = min(ppval(dz0pp,znew(1)),fac*dz); %#ok<*AGROW>
            end
            dz1 = diff(znew);
            dz1 = [dz1(1) 0.5*(dz1(1:end-1)+dz1(2:end)) dz1(end)];
            dzpp{Nz1+i} = spline(znew,smooth(znew,dz1,11));
        else
            dzpp{Nz1+i} = dz0pp;
        end
    end
    dzpos = dznew;

    % Final new grid
    znew = L1;
    while znew(end)<L2
        dz = ppval(dzpp{1},znew(end));
        for i = 2:Nz
            if active(i)
                dz = min(dz,ppval(dzpp{i},znew(end)));
            end
        end
        znew = [znew znew(end)+dz];
    end
    dznew = diff(znew); ppdznew = spline(znew,smooth(znew,[dznew(1) dznew],25));
    znew = L1;
    while znew(end)<L2
        znew = [znew znew(end)+ppval(ppdznew,znew(end))];
    end
    Lnew = spline(z,L,znew);
    Nnew = length(Lnew);
end