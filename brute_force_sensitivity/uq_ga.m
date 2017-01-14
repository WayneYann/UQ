% unsetenv DISPLAY
% matlab >&! matlab.out << EOF
%% uq.m
clear;clc;

if ispc
    % Code to run on Windows platform
    close all;
end

addpath( '../uq_library/ ');

[mech, uq_input, COMP] = deal( 'grimech30.xml', 'GRI_uncertainty.txt', 'CH4:1,O2:2,N2:7.52' );
[nrxn, N_train, ann_layers] = deal( 217, 6000, 12 );

%% ----------------------------------------------------------------
[bool.new_sample.idt, bool.new_sample.sl, bool.new_sample.ext] = deal( 1, 0, 0 );
bool.generate_new_sample = bool.new_sample.idt+ bool.new_sample.ext+ bool.new_sample.sl;
bool.psr = 1;
[bool.optimization, bool.multiobject] = deal( 1 , 1 );

p = 20;
% T_list_idt = 1050:100:1850;
T_list_idt = 1000:100:1800;
T_list_sl = 300;
T_list_ext = 300;

index = 1:nrxn;
m = length(index);

N_estimate = N_train*5; % Number for traing
sigma = readrateuq(index, uq_input); % UF = 3 sigma

case_origin.idt = Create_Object( T_list_idt, ann_layers, bool.new_sample.idt );
case_origin.sl = Create_Object( T_list_sl, ann_layers, bool.new_sample.sl );
case_origin.ext = Create_Object( T_list_ext, ann_layers, bool.new_sample.ext );
case_surrogate = case_origin;

case_origin.para.mech = mech;
case_origin.para.COMP = COMP;
case_origin.para.N_train = N_train;

sup = 1; % sup bound for the optimization
%% disp('Generate samples ...');
X0 = generate_sample(m, N_train, sup*sigma);
X1 = generate_sample(m, N_estimate, sigma);
if bool.generate_new_sample
    disp(['Generate samples', char(10)]);
    source_index = fullfile('data', 'samples_index.txt');    
    % dlmwrite(source_sample, X0(:,1:m),'delimiter','\t','precision','%.10f');
    dlmwrite(source_index, index(1:m),'delimiter','\t');
end

if bool.generate_new_sample > 0 || bool.optimization
    ncpus = feature('numCores');
    pool1 = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(pool1)
        distcomp.feature( 'LocalUseMpiexec', true );
        disp( ['ncpus= ', num2str(ncpus)] );
        parpool('local', ncpus);
    end
end

[case_origin.idt.net_list, case_origin.idt.mean_list,case_origin.idt.std_list,case_origin.idt.pdf_list] ...
    = run_sample( p, case_origin, X0, X1, 'idt' );
[case_origin.sl.net_list, case_origin.sl.mean_list,case_origin.sl.std_list,case_origin.sl.pdf_list] ...
    = run_sample( p, case_origin, X0, X1, 'sl' );
if case_origin.ext.bool.eval
    [case_origin.ext.net_list, case_origin.ext.mean_list,case_origin.ext.std_list,case_origin.ext.pdf_list] ...
        = run_sample( p, case_origin, X0, X1, 'ext' );
end

if bool.optimization
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
    ObjectiveFunction = @(x)multiobjective(x, case_origin, sigma, r, metric);
    options.FitnessLimit = 1.e-5;
    options.TolFun = 1.e-5;
    options.Generations = 600;
    %     options.PopulationSize = 400;

    if bool.multiobject
        [x,fval,flag,output,population] = gamultiobj( ObjectiveFunction, nvars, [], [], [], [], LB, UB, ConstraintFunction, options );
    else
        [x,fval,flag,output,population] = ga( ObjectiveFunction, nvars, [], [], [], [], LB, UB, ConstraintFunction, options );
    end

    case_surrogate.x = x;
    case_surrogate.fval = fval;

    case_name = [ 'case_p_', num2str(p), '_TVD.mat' ];
    save(case_name);
end