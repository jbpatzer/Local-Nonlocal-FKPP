function wB = WB(z,C,alpha)
kp = sqrt(sqrt(1/alpha)+1); km = sqrt(sqrt(1/alpha)-1);
ex = exp(kp*z); cz = cos(km*z); sz = sin(km*z);
wB = 1 + ex.*(C(5)*cz + C(6)*sz) + (C(7)*cz+ C(8)*sz)./ex;