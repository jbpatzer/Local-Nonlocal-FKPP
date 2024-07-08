function [dmat, dmat2] = setupdmat(ds)
%Build three point first and second derivative matrices

n = 1+length(ds);
dmat2 = spalloc(n,n,3*n);

for j = 2:n-1
    ds1 = ds(j-1); ds2 = ds(j);
    dmat2(j,j-1) = 2/(ds1*(ds1 + ds2));  %#ok<*SPRIX>
    dmat2(j,j) = -2/(ds1*ds2);
    dmat2(j,j+1) = 2/(ds2*(ds1 + ds2));
end
dmat2(n,n-2:n) = dmat2(n-1,n-2:n);
dmat2(1,1:3) = dmat2(2,1:3);

dmat = spalloc(n,n,3*n);
dmat(1,1) = -(2*ds(1)+ds(2))/(ds(1)+ds(2))/ds(1);
dmat(1,2) = 1/ds(1)+1/ds(2);
dmat(1,3) = -ds(1)/ds(2)/(ds(1)+ds(2));
for j = 2:n-1
    ds1 = ds(j-1); ds2 = ds(j);
    dmat(j,j-1) = -ds2/ds1/(ds1+ds2);
    dmat(j,j) = 1/ds1-1/ds2;
    dmat(j,j+1) = ds1/(ds1+ds2)/ds2;
end
ds1 = ds(n-2); ds2 = ds(n-1);
dmat(n,n-2) = ds2/(ds1+ds2)/ds1;
dmat(n,n-1) = -(ds1+ds2)/ds1/ds2;
dmat(n,n) = (ds1+2*ds2)/(ds1+ds2)/ds2;