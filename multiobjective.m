function y = multiobjective(x, case_origin, sigma, r, bool_psr)
    
    metric = 'KL-divergence'; 
    % 'Hellinger distance' 'Total variation distance' 
    % 'KL-divergence' 'Jensen–Shannon divergence'
    
    y_idt = zeros(1,length(case_origin.idt_net_list));
    for case_index =1:length(case_origin.idt_net_list)
        f = @(x)surrogate_it(x, case_index, case_origin, sigma, r, metric);
        y_idt(case_index) = f(x);
    end
    y(1) = mean(y_idt);
    
%     y = mean(y_idt);
    
    y_sl = zeros(1,length(case_origin.sl_net_list));
    for case_index =1:length(case_origin.sl_net_list)
        f = @(x)surrogate_sl(x, case_index, case_origin, sigma, r, metric);
        y_sl(case_index) = f(x);
    end
    y(2) = mean(y_sl);
    
    if bool_psr;
        y_ext = zeros(1,length(case_origin.ext_net_list));
        for case_index =1:length(case_origin.ext_net_list)
            f = @(x)surrogate_ext(x, case_index, case_origin, sigma, r, metric);
            y_ext(case_index) = f(x);
        end
        y(3) = mean(y_ext);
    end
    
end