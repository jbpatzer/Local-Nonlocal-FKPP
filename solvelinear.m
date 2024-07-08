function C = solvelinear(Nz,alpha,z)
N = 8*Nz;
M = zeros(N); b = zeros(N,1);
M(1,1) = 1; M(2,2) = 1;
M(3:6,1:8) = bcs(0,alpha);
b(3) = 1;
for k = 2:Nz
    M(8*k-9:8*k-6,[8*k-7:8*k-4 8*k-11:8*k-8]) = bcs(z(2*k-2),alpha);
    M(8*k-5:8*k-2,8*k-7:8*k) = bcs(z(2*k-1),alpha);

    b(8*k-9) = 1; b(8*k-5) = 1;
end
M(N-1,N-1) = 1; M(N,N) = 1;
C = (M\b)';

%%
function M = bcs(z,a)
kp = sqrt(sqrt(1/a)+1); km = sqrt(sqrt(1/a)-1);
r2 = sqrt(2); exA = exp(r2*z); exmA = 1/exA; 
exB = exp(kp*z); exmB = 1/exB; cz = cos(km*z); sz = sin(km*z);
M = [[exA, z*exA,exmA,z*exmA, -exB*cz,-exB*sz,-exmB*cz,-exmB*sz]
[r2*exA,exA*(r2*z + 1), -r2*exmA,-exmA*(r2*z - 1),-exB*(kp*cz - km*sz),-exB*(km*cz + kp*sz),exmB*(kp*cz + km*sz),-exmB*(km*cz - kp*sz)]
[2*exA,2*exA*(z + r2),2*exmA,2*exmA*(z - r2),exB*(km^2*cz - kp^2*cz + 2*km*kp*sz), -exB*(kp^2*sz - km^2*sz + 2*km*kp*cz), -exmB*(kp^2*cz - km^2*cz + 2*km*kp*sz),exmB*(km^2*sz - kp^2*sz + 2*km*kp*cz)]
[2*r2*exA, 2*exA*(r2*z + 3), -2*r2*exmA, -2*exmA*(r2*z - 3), -exB*(kp^3*cz + km^3*sz - 3*km^2*kp*cz - 3*km*kp^2*sz), exB*(km^3*cz - kp^3*sz - 3*km*kp^2*cz + 3*km^2*kp*sz), exmB*(kp^3*cz - km^3*sz - 3*km^2*kp*cz + 3*km*kp^2*sz), exmB*(km^3*cz + kp^3*sz - 3*km*kp^2*cz - 3*km^2*kp*sz)]]; 
