function [ net_list, mean_list, std_list, pdf_list ] =  run_sample( p, case_origin, X0, X1, type )

mech = case_origin.para.mech;
COMP = case_origin.para.COMP;
N_train = case_origin.para.N_train;
ncpus = feature('numCores');
%% 
source_sample = fullfile('data', 'samples.txt');
disp(['Evaluate ', type  ,' samples']);
%% Ignition delay time
figure;
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

source_out = fullfile('data', ['samples_out_',type,'.txt']);   %ignition delay time
for i=1:length(obj.T_list)
    T = obj.T_list(i);
    destination_sample = fullfile( 'data', ['samples_',type,'_', num2str(p),'_', num2str(T),'_K.txt' ]);
    destination_out = fullfile( 'data', ['samples_out_',type,'_', num2str(p),'_', num2str(T),'_K.txt' ]);
    if obj.bool.new_sample
        parfor i=1:feature('numCores')
            pause(i*10);
            disp( ['Starting cpu # ', num2str(i)] );
            switch type
                case 'idt'
                        system( ['python ../uq_library/run_samples_idt.py '  num2str(p), ' ', num2str(T), ' ', COMP, ' ',mech,' ',...
                    num2str(1+ ceil( N_train/ncpus )*(i-1)),' ',num2str(min( ceil( N_train/ncpus )*i, N_train)) ]);
                case 'sl'
                        system( ['python ../uq_library/run_samples_SL.py '  num2str(p), ' ', num2str(T), ' ', COMP, ' ',mech,' ',...
                    num2str(1+ ceil( N_train/ncpus )*(i-1)),' ',num2str(min( ceil( N_train/ncpus )*i, N_train)) ]);
                case 'ext'
                        [~, ~, ~] = uq_ext(mech, p, T, COMP, index, X0, [destination_out, ...
                        '_', num2str(1+ ceil( N_train/ncpus )*(i-1)) ,...
                        '_',num2str(min( ceil( N_train/ncpus )*i, N_train))], ...
                        (1+ ceil( N_train/ncpus )*(i-1)),  (min( ceil( N_train/ncpus )*i, N_train)) );
            end
          end
        DATA = [];
        for i=1:ncpus
            subdata = dlmread( [source_out, '_', ...
                num2str(1+ ceil( N_train/ncpus )*(i-1)) ,'_',num2str(min( ceil( N_train/ncpus )*i, N_train))] );
            DATA = [DATA; subdata];
        end
        copyfile( source_sample, destination_sample);
        dlmwrite(destination_out, DATA,'delimiter','\t','precision','%.20f');
    end

    switch type
        case 'idt'
            Y0 = log10( dlmread( destination_out ) );
        case 'sl'
            Y0 = dlmread( destination_out );
        case 'ext'
            Y = log10(dlmread( destination_out ));
            Y0 = Y(:,1); %extinction time at turning point
    end

    net = train_ann( log(X0), Y0, obj.ann_layers_list(i) );
    obj.net_list{i} = net;

    Y1 = net(log(X1'))';
    % The mean and std can be assigned with either the original MC and the
    % ANN response surface, for validation of the surrogate subspace, we
    % shall use the response surface as our objective funtion.
    obj.mean_list(i) = mean(Y1);
    obj.std_list(i) = std(Y1);

    disp( ['p= ', num2str(p), ' T= ', num2str(T)  ,...
        ' Mean=',num2str(obj.mean_list(i)), ' std=',num2str(obj.std_list(i))  ] );

    [f, xi] = ksdensity( Y0 );
    obj.pdf_list{i} = [f; xi];
    plot(xi, f, '-' ,'DisplayName', ['Origin samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;

    [f, xi] = ksdensity( Y1 );
    plot(xi, f, '--' ,'DisplayName', ['ANN samples ', num2str(T), ' K' ], 'LineWidth',2 );
    hold all;

end
switch type
    case 'idt'
        xlabel(' log10 IDT (s) ');
    case 'sl'
        xlabel(' sl (m/s) ');
    case 'ext'
        xlabel(' extintion time (s) ');
end
ylabel(' pdf ');
legend show;
legend boxoff;
title( [ num2str(p), ' atm ', case_origin.para.COMP]);

[ net_list, mean_list, std_list, pdf_list ] = deal( obj.net_list, obj.mean_list, obj.std_list, obj.pdf_list );

end