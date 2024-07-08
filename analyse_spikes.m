function [umax,weight,maxpos] = analyse_spikes(z,L)
u = exp(L);
[Lmax,maxpos] = findpeaks(L,z,'MinPeakProminence',0.5);
[~,minpos] = findpeaks(-L,z,'MinPeakProminence',0.5);
umax = exp(Lmax);
minpos = [minpos z(end)];
if length(minpos)<length(maxpos)+1
    maxpos = maxpos(end-length(minpos)+2:end);
    umax = umax(end-length(minpos)+2:end);
end
ppu = spline(z,u); ppL = spline(z,L);
Npeaks = length(umax);
umax = zeros(1,Npeaks); weight = zeros(1,Npeaks);
for k = 1:Npeaks
    try
        z1 = fzero(@(z) ppval(ppL,z)-log(1e-10),[minpos(k) maxpos(k)]);
    catch
        z1 = minpos(k);
    end
    try
        z2 = fzero(@(z) ppval(ppL,z)-log(1e-10),[maxpos(k) minpos(k+1)]);
    catch
        z2 = minpos(k+1);
    end
    [maxpos(k),Lmax] = fminbnd(@(z) -ppval(ppL,z),z1,z2);
    umax(k) = exp(-Lmax);
    weight(k) = integral(@(z) ppval(ppu,z),z1,z2);
end