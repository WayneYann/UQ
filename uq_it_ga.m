%%uq.m
% post_ignout;
close all;
clear;clc;
% global case_origin sigma

%% ----------------------------------------------------------------
bool_new_sample = 1;
N_train = 2000;
sup = 1; % sup bound for the optimization

p = 1;
T_list = 1045;
% T_list = 1000;
COMP = 'H2:2,O2:1,N2:3.76';
%% ----------------------------------------------------------------

nrxn = 33;
% Number for traing
N_estimate = N_train*10;
% Number for evaluating samples
index = [1:nrxn];
m = length(index); 
sigma = readrateuq(index);

%% disp('Generate samples ...');
if bool_new_sample
    disp('Generate samples ... ');
    X0 = generate_sample(m, N_train, sup*sigma);
    dlmwrite('data/samples.txt', X0(:,1:m),'delimiter','\t','precision','%.10f');
    dlmwrite('data/samples_index.txt', index(1:m),'delimiter','\t');
end
X1 = generate_sample(m, N_estimate, sigma);
%% 
disp('Evaluate ignition samples ... ');
source = 'data/samples_out_ign.txt';
source_sample = 'data/samples.txt';
case_origin.T_list = T_list;
case_origin.idt_net_list = cell(1, length(T_list));
case_origin.idt =  zeros(1,length(T_list));
case_origin.mean = zeros(1,length(T_list));
case_origin.std = zeros(1,length(T_list));

case_surrogate.mean = zeros(1,length(T_list));
case_surrogate.std = zeros(1,length(T_list));

for i=1:length(T_list)
    T = T_list(i);
    destination = [ 'data/samples_out_ign_', num2str(p),'_' , num2str(T),'_K.txt' ];
    destination_sample = [ 'data/samples_', num2str(p),'_', num2str(T),'_K.txt' ];
    if bool_new_sample
    %The py script will run idt simulation
        system( ['python run_samples.py '  num2str(p), ' ', num2str(T), ' ', COMP]);  
        copyfile( source, destination);
        copyfile( source_sample, destination_sample);
    else
        X0 = dlmread(destination_sample);
    end
    %  The ignition delay is always handled in log scale
    
    IDT0 = log10( dlmread( destination ) );
    net = train_ann_idt(X0, IDT0);
    IDT1 = net(X1')';
    
    case_origin.idt_net_list{i} = net;
    case_origin.mean(i) = mean(IDT1);
    case_origin.std(i) = std(IDT1);
    
    idt = idt_hp( p, T, COMP );
    case_origin.idt(i) = log10(idt);
        
    figure(1)
    [f, xi] = ksdensity( IDT0 );
    disp( 'Origin ANN datasets' );
    disp( ['T= ', num2str(T)  ,' Mean=',num2str(case_origin.mean(i)), ' std=',num2str(case_origin.std(i))  ] );
    disp( ['Origin IDT, log10 idt = ', num2str( log10(idt))] );
    h = plot(xi, f, '-' ,'DisplayName', ['Origin samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold on;
    
    [f, xi] = ksdensity( IDT1 );
    h = plot(xi, f, '--' ,'DisplayName', ['ANN samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
end
xlabel(' log10 IDT ');
ylabel(' pdf ');
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', COMP, ' phi=1']);
% pause;


%% Optimization
disp( 'Optimizing the surrogate subspace' );
ConstraintFunction = [];
nvars = nrxn;
LB = -sup*ones(1,nvars); 
UB =  -LB;
options = gaoptimset('Display','iter');
options.FitnessLimit = 1.e-6;
options.Generations = 200;
r = generate_sample(1, 1000, 2);
% myPool = parpool();
if 0 %length( T_list ) > 1
    ObjectiveFunction = @multiobjective;
    options = gaoptimset('Display','iter', 'PlotFcn',@gaplotpareto);
    options.FitnessLimit = 1.e-6;
    [x,fval] = gamultiobj( ObjectiveFunction, nvars, [], [], [], [], LB, UB, ConstraintFunction, options );
    disp( ['T= ', num2str(T)  ,' Mean=',num2str(case_origin.mean(i)), ' std=',num2str(case_origin.std(i))  ] );
    disp('Surrogate space:')
    disp('Coeficients :');
    disp(x);
    disp('Target function :');
    disp(fval);
    [L1, L2] = size(x);
    vector = ones(length(T_list),1)*1/sqrt(2);
    dis_fval_vector = sum(fval.*fval, 2) - (fval*vector).^2;
    [dis, index] = min(dis_fval_vector);
    x_opt = x(index, :);
    [case_surrogate.mean, case_surrogate.std] = surrogate_evaluate(x_opt);
else
    for i=1:length(T_list)
        T = T_list(i);
        %% Use ga 
        ObjectiveFunction = @(x)surrogate(x, i, case_origin, sigma, r);
    %     options.PopulationSize = 1000;
        options.CrossoverFcn = {@crossoverintermediate};
        [x,fval] = ga( ObjectiveFunction, nvars, [], [], [], [], LB, UB, ConstraintFunction, options );
        if fval > 1.e10
            %% use Global Search
            gs = GlobalSearch;
            ms = MultiStart;
            % ms.UseParallel = 'always';
            x0 = 0.9*ones(1,nrxn);
            options = optimset('Algorithm','sqp', 'Display', 'iter');
            options.MaxFunEvals = 5000;
            options.FitnessLimit = 1.e-6;

            problem = createOptimProblem('fmincon','objective', ...
                 ObjectiveFunction,'x0',x0, 'lb',LB,'ub',UB, 'options',options);
            [x1, fval1] = run( gs, problem );

            if fval > fval1
                x = x1;
            end
        end

        % delete(myPool);
        disp( ['T= ', num2str(T)  ,' Mean=',num2str(case_origin.mean(i)), ' std=',num2str(case_origin.std(i))  ] );
        disp('Surrogate space:')
        % disp('Coeficients :');
        % disp(x);
        disp('Target function :');
        disp(fval);
        x_opt = x;
        [case_surrogate.mean(i), case_surrogate.std(i)] = surrogate_evaluate_single(x_opt, i, case_origin, sigma,r);
    end
end

figure();
%%
h1 = plot( 1000./case_origin.T_list, case_origin.idt, 's', 'DisplayName', 'Origin IDT', 'LineWidth', 2 );
hold all;
e1=errorbar(1000./case_origin.T_list, case_origin.mean, case_origin.std, 'o', ...
    'DisplayName', 'Detailed Mech' , 'LineWidth', 2 );
hold all;
e2=errorbar(1000./case_origin.T_list, case_surrogate.mean, case_surrogate.std, '-^', ...
    'DisplayName', 'Surrogate', 'LineWidth', 1);
hold all;

xlabel( '1000/T (1/K)' );
ylabel( 'log10 IDT (s)' );
legend( 'Origin IDT', 'Mean & Std', 'Surrogate' );
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', COMP, ' phi=1']);


% 
% figure();
% [f, xi] = ksdensity( IDT0 );
% mu_IDT = mean(IDT0);
% sigma_IDT = std(IDT0);
% sigma_ori = sigma_IDT;
% disp( 'Origin training datasets' );
% disp( ['Mean=',num2str(mu_IDT), ' std=',num2str(sigma_IDT)  ] );
% h = plot(xi, f, 'o' ,'DisplayName', ['train samples'], 'LineWidth',2 );
% hold all;
% 
% %Train ANN
% hiddenLayerSize = 10;
% net = feedforwardnet(hiddenLayerSize,'trainlm');
% net.trainParam.epochs = 300;
% net.trainParam.goal = 1e-8;
% 
% % list_ns = [500:100:900, 1000:1000:N_train];
% list_ns = N_train;
% list_stat = zeros(length(list_ns),2);
% for i =1:length(list_ns)
%     ns = list_ns(i);
%     disp( ['ANN: ', num2str(ns), ' samples'] );
%     
%     net = train(net, X0(1:ns, :)', IDT0(1:ns)');
%     IDT1 = net(X1')';
%     [f, xi] = ksdensity( IDT1 );
%     h = plot(xi, f, 'DisplayName', ['ANN ', num2str(ns), ' samples'], 'LineWidth',2 );
%     mu_IDT = mean(IDT1);
%     sigma_IDT = std(IDT1);
%     list_stat(i, :) = [mu_IDT, sigma_IDT];
%     disp( ['Mean=',num2str(mu_IDT), ' std=',num2str(sigma_IDT)  ] );
%     xlabel('log (IDT)');
%     ylabel('pdf');
%     xlim([mu_IDT - 3*sigma_IDT, mu_IDT + 3*sigma_IDT]);
%     hold all;
% end
% legend show;
% 
% % Optimization
% disp( 'Optimizing the surrogate subspace' );
% ObjectiveFunction = @surrogate;
% ConstraintFunction = [];
% nvars = nrxn;
% LB = -1.0*ones(1,nvars); 
% UB =  1.0*ones(1,nvars);
% % options = gaoptimset('PlotFcns',{@gaplotbestf,@gaplotmaxconstr},'Display','iter');
% % options.TolFun = 1.e-5;
% % options = gaoptimset(options,'InitialPopulation',ones(1,nrxn));
% [x,fval] = ga( ObjectiveFunction, nvars, [], [], [], [], LB, UB, ConstraintFunction );
% 
% disp('Surrogate space:')
% disp('Coeficients :');
% disp(x);
% disp('Target function :');
% disp(fval);


%% PLOT the convergence rate
% figure();
% h = plot( list_ns, list_stat(:,2), 'LineWidth',2 );
% xlabel('No. Samples');
% ylabel('Std');



% 
% disp('Evaluate para space / Optimization process');
% para = rand(N_para, m);
% para = (para-0.5)./0.5;
% dist= zeros(N_para,2); %mu and sigma
% X_TEMP = X1;
% 
% for i =1:N_para
%     for j =2:m
%          X_TEMP(:,j) = exp( log(X1(:,1))/log(sigma(1))*log(sigma(j))*para(i,j) );
%     end
%     X_TEMP(:,1) = exp( log(X1(:,1))/log(sigma(1))*log(sigma(1))*para(i,1) );
%     IDT_TEMP = abs( net(X_TEMP')' );
%     [mu_IDT_TEMP, sigma_IDT_TEMP] = normfit(log(IDT_TEMP));
%     dist(i,:) = [mu_IDT_TEMP-mu_IDT,  sigma_IDT_TEMP-sigma_IDT];
% end
% 
% %Train ANN
% hiddenLayerSize = 20;
% netp = fitnet(hiddenLayerSize);
% netp = train(netp, para', dist');
% 
% figure;
% semilogy(X0(:,1), IDT0,'o');
% xlabel('k1','FontSize',14); 
% ylabel('IDT','FontSize',14);
% 
% figure;
% x1 = 0:0.05:1;
% y1 = -1:0.05:1;
% [x2, y2]= meshgrid(x1,y1);
% z2 = netp([x2(:),y2(:)]')';
% p1 = surf(x1,y1,reshape(z2(:,1), length(y1),length(x1)));
% hold all;
% p2 = surf(x1,y1,reshape(z2(:,2), length(y1),length(x1)));
% hold all;
% % scatter3(para(:,1),para(:,2),dist(:,2), 12, dist(:,2),'filled');
% xlabel('a1'); ylabel('a2'); zlabel('dist sigma');
% xlim([0,1]);
% colormap(jet);
% colorbar;
% hold all;
% p = surf(x1,y1,x2*0);
% set(p,'facealpha',0.3)
% set(p1,'facealpha',0.8)
% set(p2,'facealpha',0.8)
% 
% figure;
% scatter(para(:,1), dist(:,1),'o');
% hold all;
% scatter(para(:,1), dist(:,2),'+');
% h_leg = legend('mu','sigma');
% set(h_leg,'Location', 'Best','FontSize',14);
% legend 'boxoff';
% xlabel('a1','FontSize',14); 
% ylabel('dist','FontSize',14);
% xlim([0,1]);
