function [A,B,C,D,E,F,G,H] = extract(y,Nz)
% A = y(1:Nz); B = y(Nz+1:2*Nz);
% E = y(2*Nz+1:3*Nz); F = y(3*Nz+1:4*Nz);
% C = [0 y(4*Nz+1:5*Nz-1)]; D = [0 y(5*Nz:6*Nz-2)];
% G = [y(6*Nz-1:7*Nz-3) 0]; H = [y(7*Nz-2:8*Nz-4) 0];
A = y(1:8:8*Nz-7); B = y(2:8:8*Nz-6); C = y(3:8:8*Nz-5); D = y(4:8:8*Nz-4);
E = y(5:8:8*Nz-3); F = y(6:8:8*Nz-2); G = y(7:8:8*Nz-1); H = y(8:8:8*Nz);
% if length(y)>8*Nz
%     z = [0 -exp([y(8*Nz:10*Nz-3)])];
% else
%     z = [];
% end