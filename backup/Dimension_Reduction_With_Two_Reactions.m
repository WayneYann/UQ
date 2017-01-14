%% 
% Assume two reactions are independent
clear;clc;close all;

mu = log([1E9, 5E8]);
sigma = log([2, 1.5]);
tau_0 = 1./(exp(mu(1)) + exp(mu(2)));
k1 = 10.^(8:0.001:10);
k2 = k1;
tau = logspace(log10(1.E-5/(max(k1)+ max(k2))), log10(1.E5/(min(k1)+min(k2))),...
    1000);
tau_cdf = tau;

y1 = lognpdf(k1,mu(1),sigma(1));
y2 = lognpdf(k1,mu(2),sigma(2));

figure
h1 = semilogx (k1, y1);
hold all;
h2 = semilogx (k2, y2);
hold all;

set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
% set(gca, 'ZScale', 'log')
grid on;
xlabel('k'); ylabel('pdf');
legend('k1','k2','Location','best');
legend('boxoff');
% ylim([10^-50,1])
% zlabel('p')

%% pdf of k2 + k2
fun_cdf=@(x,tau) (1-logncdf(1./(tau)-x,mu(1),sigma(1))).*lognpdf(x,mu(2),sigma(2));
fun_pdf=@(x,tau) lognpdf(1./tau-x,mu(1),sigma(1)).*lognpdf(x,mu(2),sigma(2));
cdf_tau = @(tau) integral(@(x)fun_cdf(x,tau),0,Inf,'AbsTol',1e-40,...
    'RelTol',1e-30);
pdf_tau = @(tau) integral(@(x)fun_pdf(x,tau),0,Inf,'AbsTol',1e-40,...
    'RelTol',1e-30);
cdf_tau(tau_0)

parfor i=1:length(tau)
    tau_cdf(i) = cdf_tau(tau(i));
end

tau_pdf = diff(tau_cdf)./diff(tau);
tau_pdf = [tau_pdf(1), tau_pdf];

figure
semilogx(tau, tau_cdf);
grid on;
xlabel('tau'); ylabel('cdf');
tau_mean = sum(diff(tau_cdf).*log(tau(1:length(tau)-1)))
tau_std = sum(diff(tau_cdf).*...
    (log(tau(1:length(tau)-1))-tau_mean).^2).^0.5
tau_mid = log(tau_0)

tau_lognpdf = lognpdf(tau,tau_mean,tau_std); 

figure
plot(log10(tau), tau_pdf, 'o');
hold all;
plot(log10(tau), tau_lognpdf, '--');
grid on;
xlabel('tau'); ylabel('pdf');
xlim( [log10(1./(max(k1)+ max(k2))), log10(1./(min(k1)+min(k2)))] );
legend('boxoff');





