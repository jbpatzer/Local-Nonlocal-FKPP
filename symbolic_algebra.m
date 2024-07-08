C = sym('c',[1 8]);
F = sym('f', [1 4]);
M = sym('m',[4 8]);
syms a z km kp
wA = (C(1) + C(2)*z)*exp(sqrt(2)*z) + (C(3)+C(4)*z)*exp(-sqrt(2)*z);
wB = 1 + exp(kp*z)*(C(5)*cos(km*z) + C(6)*sin(km*z)) + exp(-kp*z)*(C(7)*cos(km*z)+ C(8)*sin(km*z));
F(1) = wA-wB+1;
F(2) = diff(wA)-diff(wB);
F(3) = diff(wA,2)-diff(wB,2);
F(4) = diff(wA,3)-diff(wB,3);

Cnew = sym(zeros(1,8));
for j = 1:8
    Cnew(j) = 1;
    for k = 1:4
        M(k,j) = subs(F(k),C,Cnew);
    end
    Cnew(j) = 0;
end
M = simplify(M,10);

r2 = sqrt(2); exA = exp(r2*z); exmA = 1/exA; 
exB = exp(kp*z); exmB = 1/exB; cz = cos(km*z); sz = sin(km*z);
M0 = [[exA, z*exA,exmA,z*exmA, -exB*cz,-exB*sz,-exmB*cz,-exmB*sz]
[r2*exA,exA*(r2*z + 1), -r2*exmA,-exmA*(r2*z - 1),-exB*(kp*cz - km*sz),-exB*(km*cz + kp*sz),exmB*(kp*cz + km*sz),-exmB*(km*cz - kp*sz)]
[2*exA,2*exA*(z + r2),2*exmA,2*exmA*(z - r2),exB*(km^2*cz - kp^2*cz + 2*km*kp*sz), -exB*(kp^2*sz - km^2*sz + 2*km*kp*cz), -exmB*(kp^2*cz - km^2*cz + 2*km*kp*sz),exmB*(km^2*sz - kp^2*sz + 2*km*kp*cz)]
[2*r2*exA, 2*exA*(r2*z + 3), -2*r2*exmA, -2*exmA*(r2*z - 3), -exB*(kp^3*cz + km^3*sz - 3*km^2*kp*cz - 3*km*kp^2*sz), exB*(km^3*cz - kp^3*sz - 3*km*kp^2*cz + 3*km^2*kp*sz), exmB*(kp^3*cz - km^3*sz - 3*km^2*kp*cz + 3*km*kp^2*sz), exmB*(km^3*cz + kp^3*sz - 3*km*kp^2*cz - 3*km^2*kp*sz)]]; 

simplify(M-M0)

clear all

%%
syms a z km kp A B C D
wA = (A+B*z)*exp(-sqrt(2)*z);
wB = 1 + exp(kp*z)*(C*cos(km*z) + D*sin(km*z));
wAz = diff(wA); wAzz = diff(wAz); wAzzz = diff(wAzz); 
wBz = diff(wB); wBzz = diff(wBz); wBzzz = diff(wBzz); 
z = 0;
wA = simplify(subs(wA,z));
wB = simplify(subs(wB,z));
wAz = simplify(subs(wAz,z));
wBz = simplify(subs(wBz,z));
wAzz = simplify(subs(wAzz,z));
wBzz = simplify(subs(wBzz,z));
wAzzz = simplify(subs(wAzzz,z));
wBzzz = simplify(subs(wBzzz,z));
S = solve(wA==wB, wAz==wBz, wAzz==wBzz,wAzzz==wBzzz,A,B,C,D);
wB0 = simplify(subs(wB,{C,D},{S.C,S.D}),10);
wB0 = simplify(subs(wB0,{km,kp},{sqrt(sqrt(1/a)-1),sqrt(sqrt(1/a)+1)}),10);
disp(wB0)

alpha = logspace(-10,0,1000);
w0 = zeros(size(alpha));
for k = 1:1000
    w0(k) = double(subs(wB0,a,alpha(k)));
end
loglog(alpha,(1-(1-alpha).*w0)./alpha,alpha,1./sqrt(alpha),'--')
    