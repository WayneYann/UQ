function [surrogate_mean, surrogate_std] = surrogate_evaluate_single(x, case_index, case_origin, sigma, r)
%     global case_origin sigma
    
%      Generate one-d samples
    X = zeros( length(r), length(x) );
    for i=1:length(x)
        X(:,i) = r.^ (x(i).*sigma(i)./2);
    end
    
    IDT1 = case_origin.idt_net_list{case_index}(log(X'))';
    surrogate_mean = mean(IDT1);
    surrogate_std = std(IDT1);
    
    %      Debug
    if 1
        figure();
        [f, xi] = ksdensity( IDT1 );
        h = plot(xi, f, '-' , 'LineWidth',2 );
        title('Optimized Surrogate Subpsace');        
        hold all;
        p = 1;
        T = 1045;
        COMP = 'H2:2,O2:1,N2:3.76';
%         COMP = 'CH4:1,O2:2,N2:7.52';
        
        %     Check ANN outside the training range
        source = 'data/samples_out_ign.txt';
        source_sample = 'data/samples.txt';
        dlmwrite( source_sample, X,'delimiter','\t','precision','%.10f');
        system( ['python run_samples.py '  num2str(p), ' ', num2str(T), ' ', COMP]);  
        IDT1 = log10( dlmread(source) );
        [f, xi] = ksdensity( IDT1 );
        h = plot(xi, f, '-' , 'LineWidth',2 );
        
        legend( 'ANN Surrogate', ' re-evaluate ' );
        legend boxoff;
        hold all;
    end
    
end