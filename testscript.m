N = 250;
alpha0 = logspace(log10(0.05),-10,N);
L1 = -logspace(log10(15),log10(75));
umax = cell(1,N); weight = cell(1,N); maxpos = cell(1,N);
c = 2e-8;
savename = 'Infdata/fullBVP';
for k = 17:N
    alpha = alpha0(k); 
    [z,L] = followpath([],[],1e-2,c,alpha,alpha,0.2,L1(k),1,0.05,Inf,[],[]);
    % [z,L] = followpath(z,L,5e-8,c,alpha,alpha,0.2,L1(k),1,0.1,Inf,[],[]);
    [umax{k},weight{k},maxpos{k}] = analyse_spikes(z,L);
    save([savename '_' num2str(k)],'z','L','c','alpha')
    save(savename,'weight','umax','maxpos','alpha')
end