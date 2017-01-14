% unsetenv DISPLAY
% matlab >&! matlab.out << EOF
%% uq.m
clear;clc;
if ispc
    % Code to run on Windows platform
    close all;
end

bool.annretrain = 1;
p = 1;

addpath( '../uq_library/ ');

[mech, uq_input, COMP] = deal( 'H2Reaction_Konnov.xml', 'H2_Uncertainty_Format_LH.txt',...
    'H2:2,O2:1,N2:3.76' );
[nrxn, N_train, ann_layers] = deal( 33, 2000, 13 );
ActiveId = [6,7,9,12,15,16,17,20,22,24];
ActiveId = 1:33;

%% ----------------------------------------------------------------
[bool.new_sample.idt, bool.new_sample.sl, bool.new_sample.ext] = deal( 0, 0, 0 );
bool.generate_new_sample = bool.new_sample.idt+ bool.new_sample.ext+ bool.new_sample.sl;
[bool.optimization, bool.multiobject] = deal( 1 , 1 );
if bool.generate_new_sample > 0 || bool.optimization
    ncpus = feature('numCores');
    pool1 = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(pool1)
        distcomp.feature( 'LocalUseMpiexec', true );
        disp( ['ncpus= ', num2str(ncpus)] );
        pool1 = parpool('local', ncpus);
        pool1.IdleTimeout = 120;
    end
end
 
case_name = [ 'case_p_', num2str(p), '_noopt.mat' ];
if bool.annretrain || exist(case_name, 'file') == 0
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
    case_origin.ext.bool.eval = 1;
    case_surrogate = case_origin;
    
    case_origin.para.p = p;
    case_origin.para.mech = mech;
    case_origin.para.COMP = COMP;
    case_origin.para.N_train = N_train;
    case_origin.para.ActiveId = ActiveId;
    case_origin.para.nrxn = nrxn;
    case_origin.para.sigma = sigma;
    case_origin.para.metric = 'KL-divergence';
    case_origin.para.metric = 'Total variation distance';
    case_origin.para.metric = 'Jensen–Shannon divergence';
    % 'Hellinger distance' 'Total variation distance' 
    % 'KL-divergence' 'Jensen–Shannon divergence'

    %% disp('Generate samples ...');
    X0 = generate_sample(m, N_train, sigma);
    X1 = generate_sample(m, N_estimate, sigma);
    if bool.generate_new_sample
        disp(['Generate samples', char(10)]);
        source_index = fullfile('data', 'samples_index.txt');    
        % dlmwrite(source_sample, X0(:,1:m),'delimiter','\t','precision','%.10f');
        dlmwrite(source_index, index(1:m),'delimiter','\t');
    end

    [case_origin.idt.net_list, case_origin.idt.mean_list,case_origin.idt.std_list,case_origin.idt.pdf_list] ...
        = run_sample( p, case_origin, X0, X1, 'idt' );
    [case_origin.sl.net_list, case_origin.sl.mean_list,case_origin.sl.std_list,case_origin.sl.pdf_list] ...
        = run_sample( p, case_origin, X0, X1, 'sl' );
    if case_origin.ext.bool.eval
        [case_origin.ext.net_list, case_origin.ext.mean_list,case_origin.ext.std_list,case_origin.ext.pdf_list] ...
            = run_sample( p, case_origin, X0, X1, 'ext' );
    end
    save(case_name, 'case_origin', 'X0', 'X1');
else
    load(case_name);
end
case_origin.bool = bool;

if case_origin.bool.optimization
    
    [x_list, fval_list] = surrogate_opt( case_origin );
    case_surrogate.x = x_list{1};
    case_surrogate.fval = fval_list{1};
    
    switch case_origin.para.metric
        case 'KL-divergence'
            case_name = [ 'case_p_', num2str(p), '_KL.mat' ];
        case 'Total variation distance'
            case_name = [ 'case_p_', num2str(p), '_TVD.mat' ];
        case 'Jensen–Shannon divergence'
            case_name = [ 'case_p_', num2str(p), '_JS.mat' ];
    end
    save(case_name);
end