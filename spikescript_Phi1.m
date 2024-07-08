c0hat = 1+logspace(-14,4,1800);
w = zeros(size(c0hat));
umax0 = zeros(size(c0hat));
for k = 1:length(c0hat)
    [~,w(k),umax0(k)] = solvespike_Phi1(c0hat(k),true);
end
K = 1./sqrt(w);
figure(2),clf
semilogx(c0hat-1,K,'k','linewidth',1.4)
grid on
xlabel('$\hat{c}_0-1$','interpreter','latex')
ylabel('$K$','interpreter','latex')
figure(3),clf
loglog(c0hat-1,c0hat.*K.^2,'k','linewidth',1.4)
grid on
xlabel('$\hat{c}_0-1$','interpreter','latex')
ylabel('$c_0 = \hat{c}_0K^2$','interpreter','latex')

ppumax0 = spline(c0hat,umax0); 
ppK = spline(c0hat,K);

c0 = logspace(-1.5,2,300);
k = zeros(size(c0));
umax = zeros(size(c0));
for i = 1:length(c0)
    c1 = fzero(@(c0hat) c0(i)-ppval(ppK,c0hat)^2*c0hat,[c0hat(1) c0hat(end)]);
    umax(i) = ppval(ppK,c1)^3*ppval(ppumax0,c1);
end

figure(4),clf
loglog(c0,umax,'k',[0.02 c0],[0.02 c0].^(0.5),'--k',[c0(1) c0(end)],0.283*[1 1],'k-.','linewidth',1.4)
grid on
xlabel('$c_0$','interpreter','latex')
ylabel('$\bar{U}_{max}$','interpreter','latex'), ylim([0.15 0.3])
legend('numerical','$c_0^{1/2}$','0.283','interpreter','latex')