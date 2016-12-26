%%uq.m
% post_ignout;
close all;
clear;clc;
% global case_origin sigma


%% ----------------------------------------------------------------
bool_new_sample = 1;
N_train = 10000;
sup = 1; % sup bound for the optimization

p = 1;
T = 1000;
% COMP = 'H2:2,O2:1,N2:3.76';
COMP = 'CH4:1,O2:2,N2:7.52';
%% ----------------------------------------------------------------
nrxn = 217;
% Number for traing
N_estimate = N_train*10;
% Number for evaluating samples
index = [1:nrxn];
m = length(index); 
% sigma = readrateuq(index);
sigma = readrateuq_gri30(index);

figure();
for N_train = 1000:2000:10000
    X0 = generate_sample(m, N_train, sup*sigma);
    dlmwrite('data/samples.txt', X0(:,1:m),'delimiter','\t','precision','%.10f');
    dlmwrite('data/samples_index.txt', index(1:m),'delimiter','\t');
    source = 'data/samples_out_ign.txt';
    source_sample = 'data/samples.txt';
    
    destination = [ 'data/samples_out_ign_', num2str(p),'_' , num2str(T),'_K.txt' ];
    destination_sample = [ 'data/samples_', num2str(p),'_', num2str(T),'_K.txt' ];
    
    %The py script will run idt simulation
%     system( ['python run_samples.py '  num2str(p), ' ', num2str(T), ' ', COMP]);  
    system( ['python run_samples_gri30.py '  num2str(p), ' ', num2str(T), ' ', COMP]);  
    copyfile( source, destination);
    copyfile( source_sample, destination_sample);
    
    IDT0 = log10( dlmread( destination ) );
%     net = train_ann_idt(X0, IDT0);
    [f, xi] = ksdensity( IDT0 );
    h = plot(xi, f, '-' ,'DisplayName', [num2str(N_train), ' samples' ], 'LineWidth',2 );
    hold all;
    legend show;
end

xlabel(' log10 IDT ');
ylabel(' pdf ');
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', COMP, ' phi=1']);