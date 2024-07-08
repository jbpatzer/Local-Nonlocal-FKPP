function [U,zg,zW,ppW] = solvelinear_inteqn(alpha,z,ker,zmin)
g = 0.5*(1-1/sqrt(3)); %Gaussian quadrature points
dz0 = 0.1; Nmin = 10; NW = 100;
z = [z z(end)-zmin]; %#ok<*NASGU>
Nz = length(z)/2;
zg = cell(1,Nz); dzg = cell(1,Nz);
Mg = cell(1,Nz); F = cell(1,Nz);
N = zeros(1,Nz);
for k = 1:Nz
    if (k>=Nz-1)
        zg{Nz-k+1} = z(2*k-1); dz = 0.01*dz0;
        while zg{Nz-k+1}(1)>z(2*k)
            zg{Nz-k+1} = [zg{Nz-k+1}(1)-dz zg{Nz-k+1}];
            dz = min(dz0,1.05*dz);
        end
        zg{Nz-k+1} = zg{Nz-k+1}(end)+(zg{Nz-k+1}(end)-z(2*k))*(zg{Nz-k+1}-zg{Nz-k+1}(end))/(zg{Nz-k+1}(end)-zg{Nz-k+1}(1));
    else
        zg{Nz-k+1} = linspace(z(2*k),z(2*k-1),max(Nmin,ceil((z(2*k-1)-z(2*k))/dz0)));
    end
    N(Nz-k+1) = length(zg{Nz-k+1});
    dzg{Nz-k+1} = diff(zg{Nz-k+1});
    F{Nz-k+1} = 0.5-kernel(zg{Nz-k+1}-z(end),ker,true);
end
Ntot = sum(N); zgtot = cell2mat(zg);
for k = 1:Nz
    z1 = zg{Nz-k+1}(1:N(Nz-k+1)-1)+g*dzg{Nz-k+1}; z2 = zg{Nz-k+1}(1:N(Nz-k+1)-1)+(1-g)*dzg{Nz-k+1};
    Mleft = zeros(Ntot,N(Nz-k+1)-1); Mright = zeros(Ntot,N(Nz-k+1)-1);
    for i = 1:Ntot
        Mleft(i,:) = dzg{Nz-k+1}.*kernel(zgtot(i)-z1,ker,false);
        Mright(i,:) = dzg{Nz-k+1}.*kernel(zgtot(i)-z2,ker,false);
    end
    Mg{Nz-k+1} = zeros(Ntot,N(Nz-k+1));
    Mg{Nz-k+1}(:,1:N(Nz-k+1)-1) = 0.5*((1-g)*Mleft+g*Mright);
    Mg{Nz-k+1}(:,2:N(Nz-k+1)) = Mg{Nz-k+1}(:,2:N(Nz-k+1))+0.5*(g*Mleft+(1-g)*Mright);
end
U0 = (((1-alpha)*[Mg{:}]+alpha*speye(Ntot))\(ones(Ntot,1)-(1-alpha)*cell2mat(F)'))';

U = cell(1,Nz);
N0 = 1;
for k = 1:Nz
    U{k} = U0(N0:N0+length(zg{k})-1);
    N0 = N0 + length(zg{k});
end

if Nz>1
    Nzw = Nz-1;
    zW = cell(1,Nzw); ppW = cell(1,Nzw);
    for i = 1:Nzw
        zW{Nzw-i+1} = linspace(z(2*i+1),z(2*i),NW); 
        W = 0.5-kernel(zW{Nzw-i+1}-z(end),ker,true);
        for k = 1:Nz
            for j = 1:NW
                z1 = zg{k}(1:N(k)-1)+g*dzg{k}; z2 = zg{k}(1:N(k)-1)+(1-g)*dzg{k};
                U1 = (1-g)*U{k}(1:N(k)-1)+g*U{k}(2:N(k));
                U2 = g*U{k}(1:N(k)-1)+(1-g)*U{k}(2:N(k));
                W(j) = W(j) + 0.5*sum(dzg{k}.*(kernel(zW{Nzw-i+1}(j)-z1,ker,false).*U1+kernel(zW{Nzw-i+1}(j)-z2,ker,false).*U2));
            end
        end
        ppW{Nzw-i+1} = spline(zW{Nzw-i+1},W);
    end
else
    ppW = []; zW = [];
end