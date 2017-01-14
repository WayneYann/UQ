% unsetenv DISPLAY
% matlab >&! matlab.out << EOF
%% uq.m
clear;clc;

if ispc
    % Code to run on Windows platform
    close all;
end

addpath( '../uq_library/ ');

mech = 'grimech30.xml';
uq_input = 'GRI_uncertainty.txt';
COMP = 'CH4:1,O2:2,N2:7.52';
nrxn = 217;
N_train = 6000;
ann_layers = 12;

%% ----------------------------------------------------------------
bool_new_sample_idt = 0;
bool_new_sample_sl = 0;
bool_new_sample_ext = 0;
bool_new_sample = bool_new_sample_idt+ bool_new_sample_ext+ bool_new_sample_sl;
sup = 1; % sup bound for the optimization
bool_multiobject = 1;
bool_optimization = 1;
bool_psr = 1;

p = 1;
% T_list_idt = 1050:100:1850;
T_list_idt = 1000:100:1800;
T_list_SL = 300;
T_list_ext = 300;

index = 1:nrxn;
m = length(index);

N_estimate = N_train*5; % Number for traing
sigma = readrateuq(index, uq_input); % UF = 3 sigma

case_origin.T_list_idt = T_list_idt;
case_origin.T_list_SL = T_list_SL;
case_origin.T_list_ext = T_list_ext;

case_origin.idt_net_list = cell(1, length(T_list_idt));
case_origin.sl_net_list = cell(1, length(T_list_SL));
case_origin.ext_net_list = cell(1, length(T_list_ext));
case_origin.ann_layers_list = ones(1, length(T_list_idt)).*ann_layers;

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
%     dlmwrite(source_sample, X0(:,1:m),'delimiter','\t','precision','%.10f');
    dlmwrite(source_index, index(1:m),'delimiter','\t');
end
X1 = generate_sample(m, N_estimate, sigma);

if bool_new_sample > 0 || bool_optimization
    ncpus = feature('numCores');
    pool1 = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(pool1)
        distcomp.feature( 'LocalUseMpiexec', true );
        disp( ['ncpus= ', num2str(ncpus)] );
        parpool('local', ncpus);
    end
end

%% 
disp(['Evaluate ignition delay samples']);
%% Ignition delay time
figure;
source_out = fullfile('data', 'samples_out_idt.txt');   %ignition delay time
for i=1:length(T_list_idt)
    T = T_list_idt(i);
    destination_sample = fullfile( 'data', ['samples_idt_', num2str(p),'_', num2str(T),'_K.txt' ]);
    destination_out = fullfile( 'data', ['samples_out_idt_', num2str(p),'_', num2str(T),'_K.txt' ]);
    if bool_new_sample_idt
        % The py script will run idt  simulation
        parfor i=1:ncpus
            pause(i*10);
            disp( ['Starting cpu # ', num2str(i)] );
            system( ['python run_samples_idt.py '  num2str(p), ' ', num2str(T), ' ', COMP, ' ',mech,' ',...
                num2str(1+ ceil( N_train/ncpus )*(i-1)),' ',num2str(min( ceil( N_train/ncpus )*i, N_train)) ]);
        end
        DATA = [];
        for i=1:ncpus
            subdata = dlmread( [source_out, '_', ...
                num2str(1+ ceil( N_train/ncpus )*(i-1)) ,'_',num2str(min( ceil( N_train/ncpus )*i, N_train))] );
            DATA = [DATA; subdata];
        end
        copyfile( source_sample, destination_sample);
        dlmwrite(destination_out, DATA,'delimiter','\t','precision','%.20f');
    else
        X0 = dlmread(destination_sample);
end
    
    % The ignition delay is always handled in log scale
    Y0 = log10( dlmread( destination_out ) );
    net = train_ann(log(X0), Y0, case_origin.ann_layers_list(i) );
    case_origin.idt_net_list{i} = net;
    
    Y1 = net(log(X1'))';
    case_origin.idt_mean(i) = mean(Y0);
    case_origin.idt_std(i) = std(Y0);
     
    disp( ['p= ', num2str(p), ' T= ', num2str(T)  ,...
        ' Mean=',num2str(case_origin.idt_mean(i)), ' std=',num2str(case_origin.idt_std(i))  ] );
    
    [f, xi] = ksdensity( Y0 );
    case_origin.idt_pdf{i} = [f; xi];
    plot(xi, f, '-' ,'DisplayName', ['Origin samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
    
    [f, xi] = ksdensity( Y1 );
    plot(xi, f, '--' ,'DisplayName', ['ANN samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;

    %     fig_idt = get(gcf, 'Number');
    
end
xlabel(' log10 IDT (s) ');
ylabel(' pdf ');
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', COMP]);

%% flame speed
disp(['Evaluate flame speed samples']);
figure;
source_out = fullfile('data','samples_out_sl.txt');      %laminar flame speed
for i=1:length(T_list_SL)
    T = T_list_SL(i);
    destination_sample = fullfile( 'data', ['samples_sl_', num2str(p),'_', num2str(T),'_K.txt' ]);
    destination_out = fullfile( 'data', ['samples_out_sl_', num2str(p),'_', num2str(T),'_K.txt' ]);
    if bool_new_sample_sl
        % The py script will run flame speed simulation
        parfor i=1:ncpus
            pause(i*10);
            disp( ['Starting cpu # ', num2str(i)] );
            system( ['python run_samples_SL.py '  num2str(p), ' ', num2str(T), ' ', COMP, ' ',mech,' ',...
                num2str(1+ ceil( N_train/ncpus )*(i-1)),' ',num2str(min( ceil( N_train/ncpus )*i, N_train)) ]);
        end
        DATA = [];
        for i=1:ncpus
            subdata = dlmread( [source_out, '_', ...
                num2str(1+ ceil( N_train/ncpus )*(i-1)) ,'_',num2str(min( ceil( N_train/ncpus )*i, N_train))] );
            DATA = [DATA; subdata];
        end
        copyfile( source_sample, destination_sample);
        dlmwrite(destination_out, DATA,'delimiter','\t','precision','%.20f');
    else
        X0 = dlmread(destination_sample);
    end
    
    Y0 = dlmread( destination_out );
    net = train_ann(log(X0), Y0,ann_layers);
    case_origin.sl_net_list{i} = net;
    
    Y1 = net(log(X1'))';
    case_origin.sl_mean(i) = mean(Y1);
    case_origin.sl_std(i) = std(Y1);
    
    disp( 'Origin ANN datasets SL' );
    disp( ['T= ', num2str(T)  ,' Mean=',num2str(case_origin.sl_mean(i)), ' std=',num2str(case_origin.sl_std(i))  ] );
    
    [f, xi] = ksdensity( Y0 );
    case_origin.sl_pdf{i} = [f; xi];
    plot(xi, f, '-' ,'DisplayName', ['Origin samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
    
    [f, xi] = ksdensity( Y1 );
    plot(xi, f, '--' ,'DisplayName', ['ANN samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
    
    fig_sl = get(gcf, 'Number');
end
xlabel(' SL (m/s) ');
ylabel(' pdf ');
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', COMP]);

if bool_psr
    %% psr turning point
    disp(['Evaluate psr extinction samples']);
    figure;
    for i=1:length(T_list_ext)
        T = T_list_ext(i);
        destination_sample = fullfile( 'data', ['samples_ext_', num2str(p),'_', num2str(T),'_K.txt' ]);
        destination_out = fullfile( 'data', ['samples_out_ext_', num2str(p),'_', num2str(T),'_K.txt' ]);
        if bool_new_sample_ext
            parfor i=1:ncpus
                pause(i*10);
                disp( ['Starting cpu # ', num2str(i)] );
                [~, ~, ext0] = uq_ext(mech, p, T, COMP, index, X0, [destination_out, ...
                    '_', num2str(1+ ceil( N_train/ncpus )*(i-1)) ,...
                    '_',num2str(min( ceil( N_train/ncpus )*i, N_train))], ...
                    (1+ ceil( N_train/ncpus )*(i-1)),  (min( ceil( N_train/ncpus )*i, N_train)) );
            end
            
            DATA = [];
            for i=1:ncpus
                subdata = dlmread( [source_out, '_', ...
                    num2str(1+ ceil( N_train/ncpus )*(i-1)) ,'_',num2str(min( ceil( N_train/ncpus )*i, N_train))] );
                DATA = [DATA; subdata];
            end
            copyfile( source_sample, destination_sample);
            dlmwrite(destination_out, DATA,'delimiter','\t','precision','%.20f');
        else
            X0 = dlmread(destination_sample);
        end

        Y = log10(dlmread( destination_out )); %extinction time at turning point
        Y0 = Y(:,1); %extinction time at turning point
        net_ext = train_ann(log(X0), Y0,ann_layers); % use the same ANN with IDT
        case_origin.ext_net_list{i} = net_ext;

        Y1 = net_ext(log(X1'))';
        case_origin.ext_mean(i) = mean(Y1);
        case_origin.ext_std(i) = std(Y1);

        disp( 'Origin ANN datasets EXT' );
        disp( ['T= ', num2str(T)  , ...
            ' Mean=',num2str(case_origin.ext_mean(i)), ' std=',num2str(case_origin.ext_std(i))  ] );

        [f, xi] = ksdensity( Y0 );
        case_origin.ext_pdf{i} = [f; xi];
        plot(xi, f, '-' ,'DisplayName', ['Origin samples ', num2str(T), ' K' ], 'LineWidth',2 );
        hold all;

        [f, xi] = ksdensity( Y1 );
        plot(xi, f, '--' ,'DisplayName', ['ANN samples ', num2str(T), ' K' ], 'LineWidth',2 );
        hold all;

        %fig_ext = get(gcf, 'Number');
    end
    xlabel(' extintion time (s) ');
    ylabel(' pdf ');
    legend show;
    legend boxoff;
    title( [ num2str(p), ' atm ', COMP]);
end

if bool_optimization
    %% Optimization
    disp( 'Optimizing the surrogate subspace' );
    rng default % for reproducibility
    r = generate_sample(1, 1000, 2);
    %% Multi-objects
    ConstraintFunction = [];
    nvars = nrxn;
    LB = -sup*ones(1,nvars); 
    UB =  -LB;
    options = gaoptimset('Display','iter' ,'UseParallel',true);
    % options = gaoptimset('Display','iter', 'PlotFcn',@gaplotpareto,'UseParallel',true);
    metric = 'Total variation distance';
    % 'Hellinger distance' 'Total variation distance' 
    % 'KL-divergence' 'Jensen–Shannon divergence'
    ObjectiveFunction = @(x)multiobjective(x, case_origin, sigma, r, bool_psr, metric);
    options.FitnessLimit = 1.e-5;
    options.TolFun = 1.e-5;
    options.Generations = 600;
    %     options.PopulationSize = 400;

    if bool_multiobject
        [x,fval,flag,output,population] = gamultiobj( ObjectiveFunction, nvars, [], [], [], [], LB, UB, ConstraintFunction, options );
    else
        [x,fval,flag,output,population] = ga( ObjectiveFunction, nvars, [], [], [], [], LB, UB, ConstraintFunction, options );
    end

    case_surrogate.x = x;
    case_surrogate.fval = fval;

    case_name = [ 'case_p_', num2str(p), '_TVD.mat' ];
    save(case_name);
end