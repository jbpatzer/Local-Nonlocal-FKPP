# Local-Nonlocal FKPP
 
MATLAB code for the paper *Slow travelling wave solutions of the mixed local-nonlocal Fisher-KPP equation can have multiple sharp wavefronts* by Jia Yin Lee and John Billingham, submitted to *Proceedings of the Royal Society A*. Descriptions of the main subroutines follow below.

## TW.m
Solves for travelling wave solutions of the mixed local-nonlocal Fisher-KPP equation, which satisfy (1.3) subject to (1.4)

## followpath.m
Perform continuation in (alpha,c) space, solving using TW.m

To reproduce the three solutions in Figure 1
followpath([],[],1e-3,1e-10,0,0,0.1,-2,0.1,0.01,1,[],[]);
followpath([],[],1e-3,1e-10,0,0,0.1,-5,0.1,0.01,2,[],[]);
followpath([],[],1e-2,1e-10,0,0,0.1,-100,2,0.15,Inf,[],[]);

Similarly for Figures 2 and 3