%%
function [A,umax,w] = analyse_spikes2(x,u)
[~,ipeaks] = findpeaks(u,'MinPeakHeight',1);
xmax = zeros(1,2); A = zeros(1,2); umax = zeros(1,2);
if length(ipeaks)>3
    for k = 1:2
        i0 = ipeaks(end-k-1);
        iplus = i0;
        while (u(iplus)>u(iplus+1))&&(iplus<length(x)-1), iplus = iplus+1; end
        iminus = i0;
        while (u(iminus)>u(iminus-1))&&(iminus>2), iminus = iminus-1; end
        if (iminus>1)&&(iplus<length(x))
            pp = spline(x(iminus:iplus),u(iminus:iplus));
            [xmax(k),umax(k)] = fminbnd(@(x) -ppval(pp,x),x(iminus), x(iplus)); umax(k) = -umax(k);
            A(k) = integral(@(x) ppval(pp,x),x(iminus),x(iplus));
        end
    end
end
umax = umax(2); A = A(2); w = xmax(1)-xmax(2);


%%
% function w = findwavelength(x,u)
% w = []; umin = 2;
% if max(u)>umin
%     N1 = floor(0.125*length(x));
%     N2  = ceil(0.875*length(x));
%     [~,locs] = findpeaks(u(N1:N2),x(N1:N2),'MinPeakHeight',umin);
%     if length(locs)>5
%         w = diff(locs);
%         w = mean(setdiff(w,min(w)));
%     end
% end

