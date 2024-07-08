# Local-Nonlocal FKPP
 
MATLAB code for the paper *Slow travelling wave solutions of the mixed local-nonlocal Fisher-KPP equation can have multiple sharp wavefronts* by Jia Yin Lee and John Billingham, submitted to *Proceedings of the Royal Society A*. Descriptions of the main subroutines follow below.

## TW.m
Solves for travelling wave solutions of the mixed local-nonlocal Fisher-KPP equation, which satisfy (1.3) subject to (1.4)

## followpath.m
Perform continuation in $(\alpha,c)$-space, solving using TW.m


To reproduce the three solutions in Figure 1

followpath([],[],1e-3,1e-10,0,0,0.1,-2,0.1,0.01,1,[],[]);

followpath([],[],1e-3,1e-10,0,0,0.1,-5,0.1,0.01,2,[],[]);

followpath([],[],1e-2,1e-10,0,0,0.1,-100,2,0.15,Inf,[],[]);


Continuation in alpha from the solutions calculated above gives the solutions shown in Figure 2.

For example 

[z,L] = followpath([],[],1e-3,1e-10,0,0,0.1,-2,0.1,0.01,1,[],[]);

then the solution with $\alpha = $ alpha can be found using

followpath(z,L,1e-10,1e-10,1e-16,alpha,1,-2,0.1,0.01,1,[],[]);


For Figure 3, the solution with $\alpha = $ alpha can be found using

followpath([],[],1e-2,1e-9,alpha,alpha,0.025,-40,2,0.15,Inf,[],[]);


## solvespike.m
Find the asymptotic solution (3.4) that satisfies (3.6) (Figure 4).

## solvespike_Phi1.m
Solve (3.15) subject to (3.16), and find c0 (Figures 5 to 7).

## spikescript_Phi2.m
Running this uses the formulation described in Section 4(c) to find the leading order asymptotic solution for $\phi = \Phi_2$ (Figures 8, 10 and dashed lines in Figure 11)

## spikescript_inteqn.m
Running this uses the formulation described in Section 4(d) to find the leading order asymptotic solution for $\phi = \Phi_{\infty}$ (dashed lines in Figures 12)

## boundaryscript.m
Find the boundaries between regions with different numbers of spikes. 

Figure 10 

boundaryscript(2,1e-10,1e-10,-5,0.1,0.01,5,1,'Phi1')

Figure 12 

boundaryscript(Inf,1e-9,1e-9,-40,2,0.1,15,1,'PhiInf')