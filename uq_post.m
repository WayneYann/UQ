close all;
clear;clc;
close all;

p = 1;
case_name = [ 'case_p_', num2str(p), '_test.mat' ];
load(case_name);

% Pre selection
% disp('candidate solution');
if bool_psr == 1
    OptIdx = fval(:,1) < 1 & fval(:,2) < 1 &fval(:,3) < 0.2;
else
    OptIdx = fval(:,1) < 0.35 & fval(:,2) < 0.3;
end

if max(OptIdx) > 0
    disp('Opt fval');
    OptList = find(OptIdx > 0);
    for i =1:length(OptList)
        disp( ['i= ', num2str(OptList(i)), '  fval= ',num2str(fval(OptList(i), :))] );
    end
    i_opt = 12;
else
    i_opt = 1;
end

%  Samples for surrogate subspace
% 30atm~63 /50 atm~63 /100 atm~ 28
% i_opt = 28;

% switch p
%     case 30
%         i_opt = 63;
%     case 50
%         i_opt = 63;
%     case 100
%         i_opt = 28;
%     otherwise
%         i_opt = 1;
% end

x_opt = x(i_opt, :);
X = zeros( length(r), length(x_opt) );
for i=1:length(x_opt)
    X(:,i) = r.^ (x_opt(i).*sigma(i)./2);
end
    
% p = 1;
% T_list_idt = 1050:100:1850;
%% Ignition delay time
figure()
IDT_ALL = [];
IDT_surrogate = [];
for i=1:length(T_list_idt)
    T = T_list_idt(i);
    destination_sample = fullfile( 'data', ['samples_idt_', num2str(p),'_', num2str(T),'_K.txt' ]);
    destination_idt = fullfile( 'data', ['samples_out_idt_', num2str(p),'_', num2str(T),'_K.txt' ]);
    % The ignition delay is always handled in log scale
    X0 = dlmread(destination_sample);
    IDT0 = log10( dlmread( destination_idt ) );
    IDT_ALL = [IDT_ALL, IDT0];
%     net_idt = train_ann_idt(log(X0), IDT0);
%     case_origin.idt_net_list{i} = net_idt;
    net_idt = case_origin.idt_net_list{i};
    IDT1 = net_idt(log(X1'))';
    
    IDT_S = net_idt(log(X'))';
    case_surrogate.idt_mean(i) = mean(IDT_S);
    case_surrogate.idt_std(i) = std(IDT_S);
%     IDT_surrogate = [IDT_surrogate, IDT_S];
    
    case_origin.idt_mean(i) = mean(IDT0);
    case_origin.idt_std(i) = std(IDT0);
        
    disp( ['P= ', num2str(p) ,' T= ', num2str(T)  ,' Mean=',num2str(case_origin.idt_mean(i)), ' std=',num2str(case_origin.idt_std(i))  ] );
%     disp( ['Origin IDT, log10 idt = ', num2str( log10(idt))] );
   [f, xi] = ksdensity( IDT0 );
    plot(xi, f, '-.' ,'DisplayName', ['MC ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
    
    [f, xi] = ksdensity( IDT1, 'npoints', 1000 );
    plot(xi, f, '--' ,'DisplayName', ['ANN ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
    
    [f, xi] = ksdensity( IDT_S );
    plot(xi, f, '-' ,'DisplayName', ['SUB ', num2str(T), ' K' ], 'LineWidth',1 );
    hold all;
    
    fig_idt = get(gcf, 'Number');
    
end
xlabel(' log10 IDT (s) ');
ylabel(' pdf ');
xlim( [case_origin.idt_mean(end) - 3*case_origin.idt_std(end),  case_origin.idt_mean(1) + 3*case_origin.idt_std(1)] );
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', COMP]);
hold all;

%% flame speed
figure()
for i=1:length(T_list_SL)
    T = T_list_SL(i);
    destination_sample = fullfile( 'data', ['samples_sl_', num2str(p),'_', num2str(T),'_K.txt' ]);
    destination_sl = fullfile( 'data', ['samples_out_sl_', num2str(p),'_', num2str(T),'_K.txt' ]);
    
    X0 = dlmread(destination_sample);
    SL0 = dlmread( destination_sl );
    net_sl = case_origin.sl_net_list{i};
    SL1 = net_sl(log(X1'))';
    SL_surrogate = net_sl( log(X') )';
    case_surrogate.sl_mean(i) = mean(SL_surrogate);
    case_surrogate.sl_std(i) = std(SL_surrogate);
     
    case_origin.sl_mean(i) = mean(SL1);
    case_origin.sl_std(i) = std(SL1);
    
%     sl = sl( p, T, COMP );
%     case_origin.sl(i) = sl;
    
    disp( ['P= ', num2str(p) ,' T= ', num2str(T)  ,' Mean=',num2str(case_origin.sl_mean(i)), ' std=',num2str(case_origin.sl_std(i))  ] );
%     disp( ['Origin IDT, log10 idt = ', num2str( sl)] );
    [f, xi] = ksdensity( SL0 );
    plot(xi, f, '-.' ,'DisplayName', ['MC ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
    
    [f, xi] = ksdensity( SL1, 'npoints', 1000 );
    plot(xi, f, '--' ,'DisplayName', ['ANN ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
    
    [f, xi] = ksdensity( SL_surrogate );
    plot(xi, f, '-' ,'DisplayName', ['SUB ', num2str(T), ' K' ], 'LineWidth',1 );
    hold all;
    
    fig_sl = get(gcf, 'Number');
    
end
xlabel(' SL (m/s) ');
xlabel('laminar flame speed (m/s)')
ylabel(' pdf ');
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', COMP]);
hold all;

if bool_psr;
%% psr turning point
figure()
for i=1:length(T_list_ext)
    T = T_list_ext(i);
    destination_sample = fullfile( 'data', ['samples_ext_', num2str(p),'_', num2str(T),'_K.txt' ]);
    destination_ext = fullfile( 'data', ['samples_out_ext_', num2str(p),'_', num2str(T),'_K.txt' ]);
    
    X0 = dlmread(destination_sample);
    EXT = dlmread( destination_ext ); %extinction time at turning point
    EXT0 = log10(EXT(:,1)); %extinction time at turning point
    TEMP_init = EXT(:,3);
    net_ext =case_origin.ext_net_list{i}; % use the same ANN with IDT
    
    EXT1 = net_ext(log(X1'))';
    EXT_surrogate = net_ext(log(X'))';
    case_surrogate.ext_mean(i) = mean(EXT_surrogate);
    case_surrogate.ext_std(i) = std(EXT_surrogate);
     
    case_origin.ext_mean(i) = mean(EXT1);
    case_origin.ext_std(i) = std(EXT1);
    
%     case_origin.ext(i) = log10(ext0);
    
    disp( 'Origin ANN datasets EXT' );
    disp( ['P= ', num2str(p),' T= ', num2str(T)  ,' Mean=',num2str(case_origin.ext_mean(i)), ' std=',num2str(case_origin.ext_std(i))  ] );
%     disp( ['Origin ext, log10 ext = ', num2str( case_origin.ext(i) )] );
    [f, xi] = ksdensity( EXT0 );
    plot(xi, f, '-.' ,'DisplayName', ['MC ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
    
    [f, xi] = ksdensity( EXT1, 'npoints', 1000 );
    plot(xi, f, '--' ,'DisplayName', ['ANN ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
    
    [f, xi] = ksdensity( EXT_surrogate );
    plot(xi, f, '-' ,'DisplayName', ['MC ', num2str(T), ' K' ], 'LineWidth',1 );
    hold all;

    fig_ext = get(gcf, 'Number');
end
xlabel(' extintion time (s) ');
ylabel(' pdf ');
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', COMP]);
hold all;
end

%% plot ignition delay time
figure();
% h1 = plot( 1000./case_origin.T_list_idt, case_origin.idt, 's', 'DisplayName', 'Origin IDT', 'LineWidth', 2 );
% hold all;
e1=errorbar(1000./case_origin.T_list_idt, case_origin.idt_mean, case_origin.idt_std, '-o', ...
    'DisplayName', 'Detailed Mech' , 'LineWidth', 2 );
hold all;
e2=errorbar(1000./case_origin.T_list_idt, case_surrogate.idt_mean, case_surrogate.idt_std, '--^', ...
    'DisplayName', 'Surrogate', 'LineWidth', 2);
hold all;
xlabel( '1000/T (1/K)' );
ylabel( 'log10 IDT (s)' );
legend( 'Origin Mean & Std', 'Surrogate Mean & Std' );
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

if bool_psr;
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
end