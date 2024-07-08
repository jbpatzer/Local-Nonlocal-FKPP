function wA= WA(z,C)
ex = exp(sqrt(2)*z);
wA = (C(1) + C(2)*z).*ex+ (C(3)+C(4)*z)./ex;