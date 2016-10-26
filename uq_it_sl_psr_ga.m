%% uq.m
close all;
clear;clc;

%% ----------------------------------------------------------------
bool_new_sample_idt = 0;
bool_new_sample_sl = 0;
bool_new_sample_ext = 0;
bool_new_sample = bool_new_sample_idt+ bool_new_sample_ext+ bool_new_sample_sl;
N_train = 2000;
sup = 1; % sup bound for the optimization
bool_multiobject = 1;

p = 1;
T_list_idt = 1000:100:1800;
%  T_list_idt = 1045;
T_list_SL = 300;
T_list_ext = 1045;
COMP = 'H2:2,O2:1,N2:3.76';

nrxn = 33;
index = [1:nrxn];
m = length(index); 

N_estimate = N_train*5; % Number for traing
sigma = readrateuq(index);
case_origin.T_list_idt = T_list_idt;
case_origin.T_list_SL = T_list_SL;
case_origin.T_list_ext = T_list_ext;
case_origin.idt_net_list = cell(1, length(T_list_idt));
case_origin.sl_net_list = cell(1, length(T_list_SL));
case_origin.ext_net_list = cell(1, length(T_list_ext));
case_origin.idt =  zeros(1,length(T_list_idt));
case_origin.sl =  zeros(1,length(T_list_SL));
case_origin.ext =  zeros(1,length(T_list_ext));
case_origin.idt_mean = zeros(1,length(T_list_idt));
case_origin.idt_std = zeros(1,length(T_list_idt));
case_origin.sl_mean = zeros(1,length(T_list_SL));
case_origin.sl_std = zeros(1,length(T_list_SL));
case_origin.ext_mean = zeros(1,length(T_list_ext));
case_origin.ext_std = zeros(1,length(T_list_ext));

% case_surrogate.x = zeros(length(T_list), m);
% case_surrogate.fval = zeros(1,length(T_list));
case_surrogate.idt_mean = zeros(1,length(T_list_idt));
case_surrogate.idt_std = zeros(1,length(T_list_idt));
case_surrogate.sl_mean = zeros(1,length(T_list_SL));
case_surrogate.sl_std = zeros(1,length(T_list_SL));
case_surrogate.ext_mean = zeros(1,length(T_list_ext));
case_surrogate.ext_std = zeros(1,length(T_list_ext));

case_surrogate.idt_mean = zeros(1,length(T_list_idt));
case_surrogate.idt_std = zeros(1,length(T_list_idt));
case_surrogate.sl_mean = zeros(1,length(T_list_SL));
case_surrogate.sl_std = zeros(1,length(T_list_SL));
case_surrogate.ext_mean = zeros(1,length(T_list_ext));
case_surrogate.ext_std = zeros(1,length(T_list_ext));

%% disp('Generate samples ...');
if bool_new_sample
    disp('Generate samples ... ');
    X0 = generate_sample(m, N_train, sup*sigma);
    dlmwrite('data/samples.txt', X0(:,1:m),'delimiter','\t','precision','%.10f');
    dlmwrite('data/samples_index.txt', index(1:m),'delimiter','\t');
end
X1 = generate_sample(m, N_estimate, sigma);
%% 
disp('Evaluate ignition & flame speed samples ... ');
source_sample = 'data/samples.txt';
source_idt = 'data/samples_out_idt.txt';   %ignition delay time
source_sl = 'data/samples_out_sl.txt';      %laminar flame speed
source_ext = 'data/samples_out_ext.txt';      %psr turning point
%% Ignition delay time
figure()
for i=1:length(T_list_idt)
    T = T_list_idt(i);
    destination_sample = [ 'data/samples_idt_', num2str(p),'_', num2str(T),'_K.txt' ];
    destination_idt = [ 'data/samples_out_idt_', num2str(p),'_' , num2str(T),'_K.txt' ];
    if bool_new_sample_idt
        % The py script will run idt  simulation
        system( ['python run_samples_idt_h2.py '  num2str(p), ' ', num2str(T), ' ', COMP]);  
        copyfile( source_sample, destination_sample);
        copyfile( source_idt, destination_idt);
    else
        X0 = dlmread(destination_sample);
    end
    
    % The ignition delay is always handled in log scale
    IDT0 = log10( dlmread( destination_idt ) );
    net_idt = train_ann_idt(log(X0), IDT0);
    IDT1 = net_idt(log(X1'))';
    case_origin.idt_net_list{i} = net_idt;
    
    case_origin.idt_mean(i) = mean(IDT1);
    case_origin.idt_std(i) = std(IDT1);
    idt = idt_hp( p, T, COMP );
    case_origin.idt(i) = log10(idt);
        
    [f, xi] = ksdensity( IDT0 );
    disp( 'Origin ANN datasets' );
    disp( ['T= ', num2str(T)  ,' Mean=',num2str(case_origin.idt_mean(i)), ' std=',num2str(case_origin.idt_std(i))  ] );
    disp( ['Origin IDT, log10 idt = ', num2str( log10(idt))] );
    plot(xi, f, '-' ,'DisplayName', ['Origin samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold on;
    [f, xi] = ksdensity( IDT1 );
    plot(xi, f, '--' ,'DisplayName', ['ANN samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
end
xlabel(' log10 IDT (s) ');
ylabel(' pdf ');
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', COMP]);

%% flame speed
figure()
for i=1:length(T_list_SL)
    T = T_list_SL(i);
    destination_sample = [ 'data/samples_sl_', num2str(p),'_', num2str(T),'_K.txt' ];
    destination_sl = [ 'data/samples_out_sl_', num2str(p),'_' , num2str(T),'_K.txt' ];
    if bool_new_sample_sl
        % The py script will run idt & flame speed simulation
        system( ['python run_samples_SL_h2.py '  num2str(p), ' ', num2str(T), ' ', COMP]);  
        copyfile( source_sample, destination_sample);
        copyfile( source_sl, destination_sl);
    else
        X0 = dlmread(destination_sample);
    end
    
    SL0 = dlmread( destination_sl );
    net_sl = train_ann_idt(log(X0), SL0);
    case_origin.sl_net_list{i} = net_sl;
    
    SL1 = net_sl(log(X1'))';
    case_origin.sl_mean(i) = mean(SL1);
    case_origin.sl_std(i) = std(SL1);
    
%     sl = sl( p, T, COMP );
%     case_origin.sl(i) = sl;
    
    [f, xi] = ksdensity( SL0 );
    disp( 'Origin ANN datasets SL' );
    disp( ['T= ', num2str(T)  ,' Mean=',num2str(case_origin.sl_mean(i)), ' std=',num2str(case_origin.sl_std(i))  ] );
%     disp( ['Origin IDT, log10 idt = ', num2str( sl)] );
    plot(xi, f, '-' ,'DisplayName', ['Origin samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold on;
    [f, xi] = ksdensity( SL1 );
    plot(xi, f, '--' ,'DisplayName', ['ANN samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
end
xlabel(' SL (m/s) ');
ylabel(' pdf ');
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', COMP]);

%% psr turning point
figure()
for i=1:length(T_list_ext)
    T = T_list_ext(i);
    destination_sample = [ 'data/samples_ext_', num2str(p),'_', num2str(T),'_K.txt' ];
    destination_ext = [ 'data/samples_out_ext_', num2str(p),'_' , num2str(T),'_K.txt' ];
    if bool_new_sample_ext
        [~, ~, ext0] = uq_ext(p, T, COMP, index, X0, destination_ext);
        copyfile( source_sample, destination_sample);
%         copyfile( source_ext, destination_ext);
    else
        X0 = dlmread(destination_sample);
    end
    
    EXT0 = log10(dlmread( destination_ext )); %extinction time at turning point
    net_ext = train_ann_idt(log(X0), EXT0); % use the same ANN with IDT
    case_origin.ext_net_list{i} = net_ext;
    
    EXT1 = net_ext(log(X1'))';
    case_origin.ext_mean(i) = mean(EXT1);
    case_origin.ext_std(i) = std(EXT1);
    
%     case_origin.ext(i) = log10(ext0);
    
    [f, xi] = ksdensity( EXT0 );
    disp( 'Origin ANN datasets EXT' );
    disp( ['T= ', num2str(T)  ,' Mean=',num2str(case_origin.ext_mean(i)), ' std=',num2str(case_origin.ext_std(i))  ] );
%     disp( ['Origin ext, log10 ext = ', num2str( case_origin.ext(i) )] );
    plot(xi, f, '-' ,'DisplayName', ['Origin samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold on;
    [f, xi] = ksdensity( EXT1 );
    plot(xi, f, '--' ,'DisplayName', ['ANN samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
end
xlabel(' extintion time (s) ');
ylabel(' pdf ');
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', COMP]);

%% Optimization
disp( 'Optimizing the surrogate subspace' );

%% Multi-objects
if bool_multiobject
    ConstraintFunction = [];
    nvars = nrxn;
    LB = -sup*ones(1,nvars); 
    UB =  -LB;
    options = gaoptimset('Display','iter');
%     options = gaoptimset('Display','iter', 'PlotFcn',@gaplotpareto);
    options.FitnessLimit = 1.e-6;
    options.Generations = 200;
    r = generate_sample(1, 1000, 2);
    
    ObjectiveFunction = @(x)multiobjective(x, case_origin, sigma, r);
    options.FitnessLimit = 1.e-6;
    [x,fval] = gamultiobj( ObjectiveFunction, nvars, [], [], [], [], LB, UB, ConstraintFunction, options );
    [dis, index] = min(sum(fval,2));
    disp('Surrogate space:')
    disp('Coeficients :');
    disp(x(index, :));
    disp('Target function :');
    disp(fval(index,:));
    x_opt = x(index, :);
    
    [case_surrogate.idt_mean, case_surrogate.idt_std] = surrogate_evaluate(x_opt, case_origin.idt_net_list, sigma, r);
    [case_surrogate.sl_mean, case_surrogate.sl_std] = surrogate_evaluate(x_opt, case_origin.sl_net_list, sigma, r);
    [case_surrogate.ext_mean, case_surrogate.ext_std] = surrogate_evaluate(x_opt, case_origin.ext_net_list, sigma, r);
end

%% plot ignition delay time
figure();
h1 = plot( 1000./case_origin.T_list_idt, case_origin.idt, 's', 'DisplayName', 'Origin IDT', 'LineWidth', 2 );
hold all;
e1=errorbar(1000./case_origin.T_list_idt, case_origin.idt_mean, case_origin.idt_std, '-o', ...
    'DisplayName', 'Detailed Mech' , 'LineWidth', 2 );
hold all;
e2=errorbar(1000./case_origin.T_list_idt, case_surrogate.idt_mean, case_surrogate.idt_std, '--^', ...
    'DisplayName', 'Surrogate', 'LineWidth', 2);
hold all;
xlabel( '1000/T (1/K)' );
ylabel( 'log10 IDT (s)' );
legend( 'Origin IDT', 'Origin Mean & Std', 'Surrogate Mean & Std' );
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', COMP]);

%% plot laminar flame speed
figure();
e1=errorbar(1000./case_origin.T_list_SL, case_origin.sl_mean, case_origin.sl_std, '-o', ...
    'DisplayName', 'Detailed Mech' , 'LineWidth', 2 );
hold all;
e2=errorbar(1000./case_origin.T_list_SL, case_surrogate.sl_mean, case_surrogate.sl_std, '--^', ...
    'DisplayName', 'Surrogate', 'LineWidth', 2);
hold all;
xlabel( '1000/T (1/K)' );
ylabel( 'SL (m/s)' );
legend( 'Origin Mean & Std', 'Surrogate Mean & Std' );
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', COMP]);

%% psr extinction time at turning point
figure();
e1=errorbar(1000./case_origin.T_list_ext, case_origin.ext_mean, case_origin.ext_std, '-o', ...
    'DisplayName', 'Detailed Mech' , 'LineWidth', 2 );
hold all;
e2=errorbar(1000./case_origin.T_list_ext, case_surrogate.ext_mean, case_surrogate.ext_std, '--^', ...
    'DisplayName', 'Surrogate', 'LineWidth', 2);
hold all;
xlabel( '1000/T (1/K)' );
ylabel( 'extintion time (s)' );
legend( 'Origin Mean & Std', 'Surrogate Mean & Std' );
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', COMP]);
