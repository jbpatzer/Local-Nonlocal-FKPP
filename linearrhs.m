function f = linearrhs(y,Nz,alpha,z)
[A,B,C,D,E,F,G,H] = extract(y,Nz);
f = zeros(1,8*Nz-4);
[WA0, WAz, WAzz, WAzzz] = WA(0,A(1),B(1),0,0);
[WB0, WBz, WBzz, WBzzz] = WB(0,E(1),F(1),G(1),H(1),alpha);
f(1) = WA0-WB0; f(2) = WAz-WBz;
f(3) = WAzz-WBzz; f(4) = WAzzz-WBzzz;
for k = 2:Nz
    [WA0, WAz, WAzz, WAzzz] = WA(z(2*k-2),A(k),B(k),C(k),D(k));
    [WB0, WBz, WBzz, WBzzz] = WB(z(2*k-2),E(k-1),F(k-1),G(k-1),H(k-1),alpha);
    f(8*k-11) = WA0-WB0; f(8*k-10) = WAz-WBz;
    f(8*k-9) = WAzz-WBzz; f(8*k-8) = WAzzz-WBzzz;
    [WA0, WAz, WAzz, WAzzz] = WA(z(2*k-1),A(k),B(k),C(k),D(k));
    [WB0, WBz, WBzz, WBzzz] = WB(z(2*k-1),E(k),F(k),G(k),H(k),alpha);
    f(8*k-7) = WA0-WB0; f(8*k-6) = WAz-WBz;
    f(8*k-5) = WAzz-WBzz; f(8*k-4) = WAzzz-WBzzz;
end