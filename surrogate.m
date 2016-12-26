function y = surrogate(x, case_index, case_origin,sigma, r)

    X = zeros( length(r), length(x) );
    for i=1:length(x)
        X(:,i) = r.^ (x(i).*sigma(i)./2);
    end
    IDT1 = case_origin.idt_net_list{case_index}(log(X'))';
    mean_IDT = mean(IDT1);
    std_IDT = std(IDT1);
    
    y = (std_IDT - case_origin.std(case_index) ).^2 + ...
        (mean_IDT - case_origin.mean(case_index)).^2;

end