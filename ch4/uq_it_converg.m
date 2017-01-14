%%uq.m
% post_ignout;
close all;
clear;clc;


%% ----------------------------------------------------------------
bool_new_sample = 1;
N_train = 10000;
sup = 1; % sup bound for the optimization

p = 1;
T = 1200;
COMP = 'CH4:1,O2:2,N2:7.52';
%% ----------------------------------------------------------------
nrxn = 217;
% Number for traing
N_estimate = N_train*10;
% Number for evaluating samples
index = [1:nrxn];
m = length(index); 
sigma = readrateuq(index);

figure();
for N_train = 10:2000:10000
    X0 = generate_sample(m, N_train, sup*sigma);
    source_sample = 'data/samples.txt';
    source = 'data/samples_out_idt.txt';
    
    dlmwrite(source_sample, X0(:,1:m),'delimiter','\t','precision','%.10f');
    dlmwrite('data/samples_index.txt', index(1:m),'delimiter','\t');
    
    destination = [ 'data_test/samples_out_idt_', num2str(p),'_' , num2str(T),'_K.txt' ];
    destination_sample = [ 'data_test/samples_', num2str(p),'_', num2str(T),'_K.txt' ];
    disp( ['No of samples = ', num2str(N_train)] );
    %The py script will run idt simulation
    system( ['python run_samples_idt_ch4.py '  num2str(p), ' ', num2str(T), ' ', COMP]);  
    copyfile( source, destination);
    copyfile( source_sample, destination_sample);
    
    IDT0 = log10( dlmread( destination ) );
%     net = train_ann_idt(X0, IDT0);
    [f, xi] = ksdensity( IDT0 );
    legend off;
    h = plot(xi, f, '-' ,'DisplayName', [num2str(N_train), ' samples' ], 'LineWidth',2 );
    hold all;
    legend show;
end

xlabel(' log10 IDT ');
ylabel(' pdf ');
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', COMP, ' phi=1']);