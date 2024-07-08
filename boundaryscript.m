function boundaryscript(ker,cmin,alphamin,z0,z1,dz0,Nm,Nm0,savefile)
%% Find boundaries between regions with different numbers of spikes
%% Arc-length continuation in (log c, log alpha)-space.
% ker = kernel
% cmin = smallest value of c
% alphamin = smallest value of alpha
% initial grid z0:dz0:z1
% Nm = number of regions
% Nm0>1 for continuation
% save to savefile
D0 = 0.25;
zinit = z0:dz0:z1;
ds = 0.025;
if Nm0>1
    load([savefile '_data'],'cbound','alphabound','cstart','zstart','Lstart')
else
    cbound = cell(1,Nm); alphabound = cell(1,Nm);
    cstart = 1e-2;
    if D0>0
        [zstart,Lstart] = TW([],[],zinit,ker,cstart,0,alphamin,[]);
    else
        zstart = []; Lstart = [];
    end
end
for k = Nm0:Nm
    cbound{k} = []; alphabound{k} = alphamin; %#ok<*SAGROW>
end
figure(2), clf
options = optimset('display','iter');
for minflag = Nm0:Nm
    if alphabound{minflag}(end)==alphamin
        z0 = zstart; L0 = Lstart;
        c0 = cstart;
    end
    dsv = [0 ds];
    fac = 0.9;

    %First pass, increment c
    while (c0>=cmin)&&(abs(dsv(1))<abs(dsv(2)))
        [z1,L1,f] = TW(z0,L0,zinit,ker,c0,D0,alphabound{minflag}(end),minflag);
        if alphabound{minflag}(end)==alphamin
            z0 = z1; L0 = L1;
        end
        if f>0
            while (f>0)&&(c0>=cmin)
                c0 = c0*fac;
                [z1,L1,f] = TW(z0,L0,zinit,ker,c0,D0,alphabound{minflag}(end),minflag);
                if alphabound{minflag}(end)==alphamin
                    z0 = z1; L0 = L1;
                end
                disp([c0 f])
            end
            c1 = [c0 c0/fac];
        else
            while f<=0
                c0 = c0/fac;
                [z1,L1,f] = TW(z0,L0,zinit,ker,c0,D0,alphabound{minflag}(end),minflag);
                if alphabound{minflag}(end)==alphamin
                    z0 = z1; L0 = L1;
                end
                disp([c0 f])
            end
            c1 = [c0*fac c0];
        end
        if alphabound{minflag}(end)==alphamin
            zstart = z0; Lstart = L0;
            cstart = c0;
            if c0<cmin
                break
            end
        end
        cbound{minflag} = [cbound{minflag} fzero(@(c) fmin(z0,L0,zinit,ker,c,D0,alphabound{minflag}(end),minflag),c1,options)];
        [z0,L0] = TW(z0,L0,zinit,ker,cbound{minflag}(end),D0,alphabound{minflag}(end),minflag);
        figure(2), clf
        for j = 1:minflag-1
            loglog(alphabound{j},cbound{j},'k-','linewidth',1.4)
            hold on
        end
        loglog(alphabound{minflag},cbound{minflag},'k-','linewidth',1.4)
        xlabel('\alpha'),ylabel('c')
        grid on
        hold off
        drawnow

        if length(alphabound{minflag})>1
            dsv = [log10(cbound{minflag}(end)/cbound{minflag}(end-1)) log10(alphabound{minflag}(end)/alphabound{minflag}(end-1))];
            ds0 = norm(dsv);
            c0 = 10^(log10(cbound{minflag}(end))+dsv(1)*ds/ds0);
            alphabound{minflag} = [alphabound{minflag} 10^(log10(alphabound{minflag}(end))+dsv(2)*ds/ds0)];
        else
            alphabound{minflag}(2) = alphabound{minflag}(1)*10^ds;
            c0 = cbound{minflag}(1);
        end
        fac = 0.995;
        options = optimset('display','iter','TolX',1e-4*cbound{minflag}(end));
    end

        % Second pass increment alpha
    cbound{minflag} = [cbound{minflag} c0]; alpha0 = alphabound{minflag}(end);
    alphabound{minflag} = alphabound{minflag}(1:end-1);
    if cbound{minflag}(end)>=cmin
        while (cbound{minflag}(end)>cmin)&&...
                (abs(alphabound{minflag}(end)-alphabound{minflag}(end-1))>0.01*ds*alphabound{minflag}(end))
            % disp(abs(alphabound{minflag}(end)-alphabound{minflag}(end-1))/0.1/ds/alphabound{minflag}(end))
            options = optimset('display','iter','TolX',1e-4*alphabound{minflag}(end));
            [~,~,f] = TW(z0,L0,zinit,ker,cbound{minflag}(end),D0,alpha0,minflag);
            if f>0
                while f>0
                    alpha0 = alpha0*fac;
                    [~,~,f] = TW(z0,L0,zinit,ker,cbound{minflag}(end),D0,alpha0,minflag);
                    disp([alpha0 f])
                end
                alpha1 = [alpha0 alpha0/fac];
            else
                while f<=0
                    alpha0 = alpha0/fac;
                    [~,~,f] = TW(z0,L0,zinit,ker,cbound{minflag}(end),D0,alpha0,minflag);
                    disp([alpha0 f])
                end
                alpha1 = [alpha0*fac alpha0];
            end
            alphabound{minflag} = [alphabound{minflag} fzero(@(alpha) fmin(z0,L0,zinit,ker,cbound{minflag}(end),D0,alpha,minflag),alpha1,options)];
            [z0,L0] = TW(z0,L0,zinit,ker,cbound{minflag}(end),D0,alphabound{minflag}(end),minflag);
            figure(2), clf
            for j = 1:minflag-1
                loglog(alphabound{j},cbound{j},'k-','linewidth',1.4)
                hold on
            end
            loglog(alphabound{minflag},cbound{minflag},'k-','linewidth',1.4)
            xlabel('\alpha'),ylabel('c')
            grid on
            hold off
            drawnow

            dsv = [log10(cbound{minflag}(end)/cbound{minflag}(end-1)) log10(alphabound{minflag}(end)/alphabound{minflag}(end-1))];
            ds0 = norm(dsv);
            alpha0 = 10^(log10(alphabound{minflag}(end))+dsv(2)*ds/ds0);
            cbound{minflag} = [cbound{minflag} 10^(log10(cbound{minflag}(end))+dsv(1)*ds/ds0)];
        end
    end
    cbound{minflag} = cbound{minflag}(1:end-1);
    save([savefile '_data'],'cbound','alphabound','cstart','zstart','Lstart')

    figure(2)
    savefig([savefile '_boundaries'])
end


%%
function f = fmin(z,L,zinit,ker,c,D0,alpha,minflag)
[~,~,f] = TW(z,L,zinit,ker,c,D0,alpha,minflag);