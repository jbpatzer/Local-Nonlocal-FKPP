function [z0,L0] = followpath(z1,L1,c1,c2,a1,a2,dsmax,zmin,zmax,dz0,ker,savename,count)
%% Follow a path in alpha c parameter space finding TW solutions of the mixed local-nonlocal Fisher-KPP equation
% z1, L1 initial solution
% continuation from (a1,c1) to (a2,c2)
% dsmax = maximum step length in log a, log c parameter space
% initial grid is zinit = zmin:dz0:zmax
% ker specifies kernel
% save to files savename_count
% continue from point 'count' on the curve (count = [] to begin)

format compact
format short g

tol = 1e-2;
dsmin = 1e-6;
ds = dsmax;
D0 = 0.25;
zinit = zmin:dz0:zmax;
umax = cell(1); weight = cell(1); maxpos = cell(1); aout = []; cout = [];

if ~isempty(count)
    load([savename '_' num2str(count)],'z1','L1','s','ds')
else
    s = 1;
    count = 1;
    if (isempty(L1))&&(D0>0)
        [z1,L1] = TW([],[],zinit,ker,c1,0,a1,[]);
    end
    c = c1; alpha = a1;
    [z1,L1] = TW(z1,L1,zinit,ker,c,D0,alpha,[]);
    if ~isempty(savename)
        [umax{count},weight{count},maxpos{count}] = analyse_spikes(z1,L1);
        aout(count) = alpha; cout(count) = c;
        save([savename '_' num2str(count)],'z1','L1','c','alpha','s','ds','count')
        save(savename,'weight','umax','maxpos','aout','cout')
    end
end

f = 0;
while ((s>1e-14)||(norm(f)>tol))&&(ds>dsmin)
    if norm(f)<tol
        z0 = z1; L0 = L1;
        s = max(0,s-ds);
        ds = min(dsmax,1.05*ds);
    else
        s = s+ds; ds = 0.75*ds;
        s = s-ds;
    end
    disp([s ds])
    c = 10^(s*log10(c1)+(1-s)*log10(c2));
    alpha= 10^(s*log10(a1)+(1-s)*log10(a2));
    try
        [z1,L1,~,f] = TW(z0,L0,zinit,ker,c,D0,alpha,[]);
    catch
        f = 1e12;
    end
    drawnow
    if (~isempty(savename))&&(norm(f)<tol)
        count = count+1;
        [umax{count},weight{count},maxpos{count}] = analyse_spikes(z1,L1);
        aout(count) = alpha; cout(count) = c; %#ok<AGROW>
        save([savename '_' num2str(count)],'z1','L1','c','alpha','s','ds','count')
        save(savename,'weight','umax','maxpos','aout','cout')
    end
end
z0 = z1; L0 = L1;