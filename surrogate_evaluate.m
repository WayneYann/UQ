function [surrogate_mean, surrogate_std] = surrogate_evaluate(x, netlist, sigma, r, fig_num)

%      Generate one-d samples
    X = zeros( length(r), length(x) );
    for i=1:length(x)
        X(:,i) = r.^ (x(i).*sigma(i)./2);
    end

    surrogate_mean = zeros(1, length( netlist ) );
    surrogate_std = zeros(1, length( netlist ) );
    
    for case_index=1:length( netlist )
        IDT1 = netlist{case_index}(log(X'))';
        surrogate_mean(case_index) = mean(IDT1);
        surrogate_std(case_index) = std(IDT1);
        % Debug
        if 1
            figure(fig_num);
            [f, xi] = ksdensity( IDT1 );
            h = plot(xi, f, '-.' , 'LineWidth',2, 'DisplayName', 'Surrogate Subspace');
            title('Optimized Surrogate Subpsace');        
            hold all;
        end
    end

end