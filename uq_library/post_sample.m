function [ mean_list, std_list, pdf_list ] =  post_sample( p, case_origin, X, X1, type )

%% 
disp(['Evaluate ', type  ,' samples']);
%% Ignition delay time
figure();
switch type
    case 'idt'
        obj = case_origin.idt;
    case 'sl'
        obj = case_origin.sl;
    case 'ext'
        obj = case_origin.ext;
    otherwise
        disp( 'warning: unsupported type of problem' )
end

obj_sub = obj;

for i=1:length(obj.T_list)
    T = obj.T_list(i);
    destination_out = fullfile( 'data', ['samples_out_',type,'_', num2str(p),'_', num2str(T),'_K.txt' ]);
    switch type
        case 'idt'
            Y0 = log10( dlmread( destination_out ) );
        case 'sl'
            Y0 = dlmread( destination_out );
        case 'ext'
            Y = log10(dlmread( destination_out ));
            Y0 = Y(:,1); %Extinction time at turning point
    end

    net = obj.net_list{i};
    Y1 = net(log(X1'))';
    
    YS = net(log(X'))';
    obj_sub.mean_list(i) = mean(YS);
    obj_sub.std_list(i) = std(YS);

    disp( ['p= ', num2str(p), ' T= ', num2str(T)  ,...
        ' Mean=',num2str(obj_sub.mean_list(i)), ' std=',num2str(obj_sub.std_list(i))  ] );

    [f, xi] = ksdensity( Y0 );
    plot(xi, f, '-' ,'DisplayName', ['MC, Origin ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;

    [f, xi] = ksdensity( Y1 );
    plot(xi, f, '--' ,'DisplayName', ['ANN, Origin ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;
    
    [f, xi] = ksdensity( YS, obj.pdf_list{i}(2,:) );
    obj_sub.pdf_list{i} = [f; xi];
    plot(xi, f, '-.' ,'DisplayName', ['ANN, Surrogate Subspace ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;

end

ylabel(' pdf ');
xlim( [obj.mean_list(end) - 4*obj.std_list(end),  obj.mean_list(1) + 4*obj.std_list(1)] );
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', case_origin.para.COMP]);

figure();
e1=errorbar(1000./obj.T_list, obj.mean_list, obj.std_list, '-o', ...
    'DisplayName', 'ANN for Origin Mean & Std' , 'LineWidth', 2 );
hold all;
e2=errorbar(1000./obj.T_list, obj_sub.mean_list, obj_sub.std_list, '--^', ...
    'DisplayName', 'ANN for Surrogate Mean & Std', 'LineWidth', 2);
hold all;
xlabel( '1000/T (1/K)' );
label(type);
% legend( 'Origin Mean & Std', 'Surrogate Mean & Std' );
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', case_origin.para.COMP]);

[ mean_list, std_list, pdf_list ] = deal( obj_sub.mean_list, obj_sub.std_list, obj_sub.pdf_list );

end

function label( type )
    switch type
        case 'idt'
            xlabel(' log10 (IDT) (s) ');
        case 'sl'
            xlabel('sl (m/s) ');
        case 'ext'
            xlabel(' Extintion time (s) ');
    end
end