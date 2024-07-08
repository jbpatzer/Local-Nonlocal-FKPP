function f = kernel(y,N,Iflag)
persistent k0 r1 r2
if isempty(k0), k0 = 1/2/sqrt(pi); end
if isempty(r2), r2 = sqrt(2); end
if isempty(r1), r1 = 0.25*r2; end
if Iflag
    switch N
        case -1
            k = 0.97;
            f = 0.5*(1-(1+0.5*k*y).*exp(-y));
        case 1
            f = 0.5*(1-exp(-abs(y)));
        case 2
            f = 0.5*(1-(1+y/r2).*exp(-r2*y));
        case 3
            f = 0.5*(1-(1+y/2).*exp(-y));
        case Inf
            f = 0.5*erf(0.5*y);
    end
else
    switch N
        case -1
            k = 0.97;
            f = 0.5*(1+0.5*k*(abs(y)-1)).*exp(-abs(y));
        case 1
            f = 0.5*exp(-abs(y));
        case 2
            f = (r1+0.5*abs(y)).*exp(-r2*abs(y));
        case 3 % Psi_2
            f = 0.25*(1+abs(y)).*exp(-abs(y));
        case Inf
            f = k0*exp(-0.25*y.^2);
    end
end