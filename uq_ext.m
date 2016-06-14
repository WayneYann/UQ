%%uq.m
% post_ignout;
close all;
clear;clc;
debug_flag = 0;
m = 2; 
N_train = 100;
N_para = 1E4;
index = [13,  29];
sigma = [1.5, 2];
COMP = 'H2:2,O2:1,N2:3.76';

disp('generate samples');
X0 = generate_sample(m, N_train, sigma, debug_flag);
dlmwrite('samples.txt', X0(:,1:m),'delimiter','\t','precision','%.6f');
dlmwrite('samples_index.txt', index(1:m),'delimiter','\t');

disp('evaluate samples');
EXT = zeros(N_train, 2);
% system('run_samples.py');
P_in = 1.0;
T_in = 700;
Omega_init = 1./(2E-4);

for i =1:N_train
    EXT(i,:) = psr_extinction( P_in, T_in, COMP, Omega_init, index, X0(i,:) );
end

disp('train samples');
% the extinction residence time is also denoted as IDT
IDT0 = EXT(:,1);
% regenrate incase start from here
X0 = generate_sample(m, N_train, sigma, debug_flag);

figure(1);
subplot(2,1,1);
histfit(log(IDT0));
[mu_IDT, sigma_IDT] = normfit(log(IDT0));
xlabel('log (IDT)');
ylabel('pdf');
xlim([mu_IDT - 3*sigma_IDT, mu_IDT + 3*sigma_IDT]);
hold all;

%Train ANN
hiddenLayerSize = 20;
net = fitnet(hiddenLayerSize);
net = train(net, X0', IDT0');

disp('estimate new samples');
X1 = generate_sample(m, 2000, sigma, debug_flag);
IDT1 = abs( net(X1')' );
[mu_IDT, sigma_IDT] = normfit(log(IDT1));

subplot(2,1,2);
h = histfit(log(IDT1));
xlabel('log (IDT)');
ylabel('pdf');
xlim([mu_IDT - 3*sigma_IDT, mu_IDT + 3*sigma_IDT]);
hold all;

disp('evaluate para space');
para = rand(N_para, m);
para = (para-0.5)./0.5;
dist= zeros(N_para,2); %mu and sigma
X_TEMP = X1;

for i =1:N_para
    for j =2:m
         X_TEMP(:,j) = exp( log(X1(:,1))/log(sigma(1))*log(sigma(j))*para(i,j) );
    end
    X_TEMP(:,1) = exp( log(X1(:,1))/log(sigma(1))*log(sigma(1))*para(i,1) );
    IDT_TEMP = abs( net(X_TEMP')' );
    [mu_IDT_TEMP, sigma_IDT_TEMP] = normfit(log(IDT_TEMP));
    dist(i,:) = [mu_IDT_TEMP-mu_IDT,  sigma_IDT_TEMP-sigma_IDT];
end

%Train ANN
hiddenLayerSize = 20;
netp = fitnet(hiddenLayerSize);
netp = train(netp, para', dist');

figure;
semilogy(X0(:,1), IDT0,'o');
xlabel('k1','FontSize',14); 
ylabel('IDT','FontSize',14);

figure;
x1 = 0:0.05:1;
y1 = -1:0.05:1;
[x2, y2]= meshgrid(x1,y1);
z2 = netp([x2(:),y2(:)]')';
p1 = surf(x1,y1,reshape(z2(:,1), length(y1),length(x1)));
hold all;
p2 = surf(x1,y1,reshape(z2(:,2), length(y1),length(x1)));
hold all;
% scatter3(para(:,1),para(:,2),dist(:,2), 12, dist(:,2),'filled');
xlabel('a1'); ylabel('a2'); zlabel('dist sigma');
xlim([0,1]);
colormap(jet);
colorbar;
hold all;
p = surf(x1,y1,x2*0);
set(p,'facealpha',0.3)
set(p1,'facealpha',0.8)
set(p2,'facealpha',0.8)

figure;
scatter(para(:,1), dist(:,1),'o');
hold all;
scatter(para(:,1), dist(:,2),'+');
h_leg = legend('mu','sigma');
set(h_leg,'Location', 'Best','FontSize',14);
legend 'boxoff';
xlabel('a1','FontSize',14); 
ylabel('dist','FontSize',14);
xlim([0,1]);