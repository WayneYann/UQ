function y = multiobjective(x, case_origin)
    
    y = [];
    for obj =[case_origin.idt, case_origin.sl, case_origin.ext]
        if obj.bool.eval
            y_obj = zeros(1,length(obj.T_list));
            for case_index =1:length(obj.T_list)
                f = @(x)surrogate(x, obj.net_list{case_index}, obj.pdf_list{case_index}, ...
                    obj.std_list(case_index),case_origin.para.sigma, case_origin.para.metric);
                y_obj(case_index) = f(x);
            end
            y(end+1) = mean(y_obj);
        end
    end
    if case_origin.bool.multiobject == 0
        y = sum(y);
    end
end

function y = surrogate(x, net, pdf_origin, std_origin, x_sigma, metric)
    r= log(generate_sample(1, 1000, exp(norm(x))));
    X = zeros( length(r), length(x) );
    for i=1:length(x)
        X(:,i) = exp( r* x(i)/norm(x)*log(x_sigma(i)) );
    end
    
    debug_flag = 0;
    if (debug_flag)
        figure(10);
         j = 1;
        [f, xi] = ksdensity( log(X(:,j)) );
        plot(xi, f, '-');
        hold all;
        r= generate_sample(1, 1000, 2);
        X = zeros( length(r), length(x) );
        for i=1:length(x)
            X(:,i) = r.^ (x(i).*log(x_sigma(i))/log(2));
        end
        [f, xi] = ksdensity( log(X(:,j)) );
        plot(xi, f, '-.');     
        hold all;
    end

    Y = net(log(X'))';
    [f, xi] = ksdensity( Y, pdf_origin(2,:) );
    pdf_surrogate = [f;xi];
    std_surrogate = std(Y);
    y = func_metric( pdf_surrogate, pdf_origin, metric );
    
    if y == 0 || std_origin > 10*std_surrogate
        y = Inf;
    end
        
end

function y = func_metric(pdf_surrogate, pdf_origin, metric)
    % by setting threshhold for calculating the distance, we shall focus on 
    % caparing the important part
    epsilong = 1.e-2;
    goodIdx = pdf_surrogate(1,:)>epsilong & pdf_origin(1,:)>epsilong;
    diff_x = [diff(pdf_origin(2,:)), pdf_origin(2,end)-pdf_origin(2,end-1)] ;
    switch metric
        case 'KL-divergence'
            d2 = ( (pdf_origin(1,goodIdx)- pdf_surrogate(1,goodIdx)).* log(pdf_origin(1,goodIdx) ./ pdf_surrogate(1,goodIdx)) ) * diff_x(goodIdx)';
            y = d2;
        case 'Hellinger distance'
            d2 = ( pdf_origin(1,goodIdx).*(sqrt(pdf_origin(1,goodIdx) ./ pdf_surrogate(1,goodIdx)) -1 ).^2 .* pdf_surrogate(1,goodIdx) )* diff_x(goodIdx)';
            y = d2;
        case 'Total variation distance'
            d2 = (pdf_origin(1,goodIdx).*abs(pdf_origin(1,goodIdx) ./ pdf_surrogate(1,goodIdx) -1 )./2) * diff_x(goodIdx)';
            y = d2;
        case 'Jensen-Shannon divergence'
            pdf_M = 0.5 .*(pdf_surrogate + pdf_origin);
            d2_1 = ( pdf_origin(1,goodIdx) .* log(pdf_origin(1,goodIdx) ./ pdf_M(1,goodIdx)) )* diff_x(goodIdx)';
            d2_2 = ( pdf_surrogate(1, goodIdx) .* log(pdf_surrogate(1,goodIdx) ./ pdf_M(1,goodIdx)) )* diff_x(goodIdx)';
            y = 0.5*(d2_1 + d2_2);
        otherwise %default KL-div
            d2 = ( (pdf_origin(1,goodIdx)- pdf_surrogate(1,goodIdx)).* log(pdf_origin(1,goodIdx) ./ pdf_surrogate(1,goodIdx)) ) * diff_x(goodIdx)';
            y = d2;
    end
    
end