% unsetenv DISPLAY
% matlab >&! matlab.out << EOF
%% uq.m
clear;clc;

if ismac
    % Code to run on Mac plaform
elseif isunix
    % Code to run on Linux plaform
elseif ispc
    % Code to run on Windows platform
    close all;
end

%% ----------------------------------------------------------------
bool_new_sample_idt = 0;
bool_new_sample_sl = 0;
bool_new_sample_ext = 0;
bool_new_sample = bool_new_sample_idt+ bool_new_sample_ext+ bool_new_sample_sl;
N_train = 2000;
sup = 1; % sup bound for the optimization
bool_multiobject = 1;
bool_psr = 1;

p = 1;
% T_list_idt = 1050:100:1850;
T_list_idt = 1000:100:1800;
T_list_SL = 300;
T_list_ext = 300;
COMP = 'H2:2,O2:1,N2:3.76';

nrxn = 33;
index = [1:nrxn];
m = length(index); 

N_estimate = N_train*5; % Number for traing
sigma = readrateuq(index); % UF = 3 sigma

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
case_origin.idt_hbw = zeros(3,length(T_list_idt));
case_origin.idt_pdf = cell(1,length(T_list_idt));

case_origin.sl_mean = zeros(1,length(T_list_SL));
case_origin.sl_std = zeros(1,length(T_list_SL));
case_origin.sl_hbw = zeros(3,length(T_list_SL));
case_origin.sl_pdf = cell(1,length(T_list_SL));

case_origin.ext_mean = zeros(1,length(T_list_ext));
case_origin.ext_std = zeros(1,length(T_list_ext));
case_origin.ext_hbw = zeros(3,length(T_list_ext));
case_origin.ext_pdf = cell(1,length(T_list_ext));

case_surrogate = case_origin;

%% disp('Generate samples ...');
source_sample = fullfile('data', 'samples.txt');
source_index = fullfile('data', 'samples_index.txt');
if bool_new_sample
    disp(['Generate samples', char(10)]);
    X0 = generate_sample(m, N_train, sup*sigma);
    dlmwrite(source_sample, X0(:,1:m),'delimiter','\t','precision','%.10f');
    dlmwrite(source_index, index(1:m),'delimiter','\t');
end
X1 = generate_sample(m, N_estimate, sigma);
%% 
disp(['Evaluate ignition & flame speed samples']);
source_idt = fullfile('data', 'samples_out_idt.txt');   %ignition delay time
source_sl = fullfile('data','samples_out_sl.txt');      %laminar flame speed
source_ext = fullfile('data', 'samples_out_ext.txt');      %psr turning point

%% Ignition delay time
figure;
for i=1:length(T_list_idt)
    T = T_list_idt(i);
    destination_sample = fullfile( 'data', ['samples_idt_', num2str(p),'_', num2str(T),'_K.txt' ]);
    destination_idt = fullfile( 'data', ['samples_out_idt_', num2str(p),'_', num2str(T),'_K.txt' ]);
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
%     idt = idt_hp( p, T, COMP );
%     case_origin.idt(i) = log10(idt);
     
    disp( ['p= ', num2str(p), ' T= ', num2str(T)  ,' Mean=',num2str(case_origin.idt_mean(i)), ' std=',num2str(case_origin.idt_std(i))  ] );
%     disp( ['Origin IDT, log10 idt = ', num2str( log10(idt))] );
    
    [f, xi] = ksdensity( IDT0 );
    case_origin.idt_pdf{i} = [f; xi];
    plot(xi, f, '-' ,'DisplayName', ['Origin samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
    
    [f, xi] = ksdensity( IDT1 );
    plot(xi, f, '--' ,'DisplayName', ['ANN samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
    
    fig_idt = get(gcf, 'Number');
    % Get full width of half height
    loc = find(f > max(f)./2);
    case_origin.idt_hbw(:,i) = [ xi(loc(1)), xi(loc(end)), xi(f==max(f)) ];
    
end
xlabel(' log10 IDT (s) ');
ylabel(' pdf ');
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', COMP]);

%% flame speed
figure;
for i=1:length(T_list_SL)
    T = T_list_SL(i);
    destination_sample = fullfile( 'data', ['samples_sl_', num2str(p),'_', num2str(T),'_K.txt' ]);
    destination_sl = fullfile( 'data', ['samples_out_sl_', num2str(p),'_', num2str(T),'_K.txt' ]);
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
    
    
    disp( 'Origin ANN datasets SL' );
    disp( ['T= ', num2str(T)  ,' Mean=',num2str(case_origin.sl_mean(i)), ' std=',num2str(case_origin.sl_std(i))  ] );
%     disp( ['Origin IDT, log10 idt = ', num2str( sl)] );
    
    [f, xi] = ksdensity( SL0 );
    case_origin.sl_pdf{i} = [f; xi];
    plot(xi, f, '-' ,'DisplayName', ['Origin samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
    
    [f, xi] = ksdensity( SL1 );
    plot(xi, f, '--' ,'DisplayName', ['ANN samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
    
    fig_sl = get(gcf, 'Number');
    
    % Get full width of half height
    loc = find(f > max(f)./2);
    case_origin.sl_hbw(:,i) = [ xi(loc(1)), xi(loc(end)), xi(f==max(f)) ];
end
xlabel(' SL (m/s) ');
ylabel(' pdf ');
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', COMP]);

if bool_psr
%% psr turning point
figure;
for i=1:length(T_list_ext)
    T = T_list_ext(i);
    destination_sample = fullfile( 'data', ['samples_ext_', num2str(p),'_', num2str(T),'_K.txt' ]);
    destination_ext = fullfile( 'data', ['samples_out_ext_', num2str(p),'_', num2str(T),'_K.txt' ]);
    if bool_new_sample_ext
        [~, ~, ext0] = uq_ext(p, T, COMP, index, X0, destination_ext);
%         system( ['python run_samples_ext_h2.py ']);  
        copyfile( source_sample, destination_sample);
%         copyfile( source_ext, destination_ext);
    else
        X0 = dlmread(destination_sample);
    end
    
    EXT = log10(dlmread( destination_ext )); %extinction time at turning point
    EXT0 = EXT(:,1); %extinction time at turning point
    net_ext = train_ann_idt(log(X0), EXT0); % use the same ANN with IDT
    case_origin.ext_net_list{i} = net_ext;
    
    EXT1 = net_ext(log(X1'))';
    case_origin.ext_mean(i) = mean(EXT1);
    case_origin.ext_std(i) = std(EXT1);
    
%     case_origin.ext(i) = log10(ext0);
    
    disp( 'Origin ANN datasets EXT' );
    disp( ['T= ', num2str(T)  ,' Mean=',num2str(case_origin.ext_mean(i)), ' std=',num2str(case_origin.ext_std(i))  ] );
%     disp( ['Origin ext, log10 ext = ', num2str( case_origin.ext(i) )] );
    
    [f, xi] = ksdensity( EXT0 );
    case_origin.ext_pdf{i} = [f; xi];
    plot(xi, f, '-' ,'DisplayName', ['Origin samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
    
    [f, xi] = ksdensity( EXT1 );
    plot(xi, f, '--' ,'DisplayName', ['ANN samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
    
    fig_ext = get(gcf, 'Number');
    % Get full width of half height
    loc = find(f > max(f)./2);
    case_origin.ext_hbw(:,i) = [ xi(loc(1)), xi(loc(end)), xi(f==max(f)) ];
end
xlabel(' extintion time (s) ');
ylabel(' pdf ');
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', COMP]);
end

%% Optimization
disp( 'Optimizing the surrogate subspace' );
rng default % for reproducibility
r = generate_sample(1, 1000, 2);
%% Multi-objects
if bool_multiobject
    ConstraintFunction = [];
    nvars = nrxn;
    LB = -sup*ones(1,nvars); 
    UB =  -LB;
    options = gaoptimset('Display','iter' ,'UseParallel',true);
%     options = gaoptimset('Display','iter', 'PlotFcn',@gaplotpareto,'UseParallel',true);
    ObjectiveFunction = @(x)multiobjective(x, case_origin, sigma, r, bool_psr);
    options.FitnessLimit = 1.e-5;
    options.TolFun = 1.e-5;
    options.Generations = 300;
%     options.PopulationSize = 400;
%     options.Generations = 400;
    
    pool1 = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(pool1)
        distcomp.feature( 'LocalUseMpiexec', true )
        ncpus = feature('numCores');
        disp( ['ncpus= ', num2str(ncpus)] );
        parpool('local', ncpus);
    end
    
    [x,fval,flag,output,population] = gamultiobj( ObjectiveFunction, nvars, [], [], [], [], LB, UB, ConstraintFunction, options );
%     [x,fval,flag,output,population] = ga( ObjectiveFunction, nvars, [], [], [], [], LB, UB, ConstraintFunction, options );
    
    case_surrogate.x = x;
    case_surrogate.fval = fval;
    
    case_name = [ 'case_p_', num2str(p), '_test.mat' ];
    save(case_name);
    
end