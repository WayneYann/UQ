function y = surrogate_it(x, case_index, case_origin,sigma, r, metric)

    X = zeros( length(r), length(x) );
    for i=1:length(x)
        X(:,i) = r.^ (x(i).*sigma(i)./2);
    end
    
    IDT1 = case_origin.idt_net_list{case_index}(log(X'))';

    pdf_origin = case_origin.idt_pdf{case_index};
    [f, xi] = ksdensity( IDT1, pdf_origin(2,:) );
    pdf_surrogate = [f;xi];
    
    y = func_metric( pdf_surrogate, pdf_origin, metric );
    
%     mean_IDT = mean(IDT1);
%     std_IDT = std(IDT1);    
%     loc = find(f > max(f)./2);
%     
%     y = (std_IDT - case_origin.idt_std(case_index) ).^2 + ...
%         (mean_IDT - case_origin.idt_mean(case_index)).^2 + ...
%         (xi(loc(1)) - case_origin.idt_hbw(1,case_index)).^2 + ...
%         (xi(loc(end)) - case_origin.idt_hbw(2,case_index)).^2 + ...
%         (xi(f == max(f)) - case_origin.idt_hbw(3,case_index)).^2;
%     
%     y = y./(case_origin.idt_mean(case_index).^2);
    
end