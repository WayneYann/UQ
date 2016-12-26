function y = func_metric(pdf_surrogate, pdf_origin, metric)

    goodIdx = pdf_surrogate(1,:)>0 & pdf_origin(1,:)>0;
    switch metric
        case 'KL-divergence'
            d1 = sum( pdf_surrogate(1, goodIdx) .* log(pdf_surrogate(1,goodIdx) ./ pdf_origin(1,goodIdx)) );
            d2 = sum( pdf_origin(1,goodIdx) .* log(pdf_origin(1,goodIdx) ./ pdf_surrogate(1,goodIdx)) );
            y = abs(d2);
        case 'Hellinger distance'
            d1 = sum( (sqrt(pdf_surrogate(1,goodIdx) ./ pdf_origin(1,goodIdx)) -1 ).^2 .* pdf_origin(1,goodIdx) );
            d2 = sum( (sqrt(pdf_origin(1,goodIdx) ./ pdf_surrogate(1,goodIdx)) -1 ).^2 .* pdf_surrogate(1,goodIdx) );
            y = d2;
        case 'Total variation distance'
            d1 = sum( (abs(pdf_surrogate(1,goodIdx) ./ pdf_origin(1,goodIdx)) -1 )./2 .* pdf_origin(1,goodIdx) );
            d2 = sum( (abs(pdf_origin(1,goodIdx) ./ pdf_surrogate(1,goodIdx)) -1 )./2 .* pdf_surrogate(1,goodIdx) );
            y = d2;
        case 'Jensen–Shannon divergence'
            pdf_M = 0.5 .*(pdf_surrogate + pdf_origin);
            d2_1 = sum( pdf_origin(1,goodIdx) .* log(pdf_origin(1,goodIdx) ./ pdf_M(1,goodIdx)) );
            d2_2 = sum( pdf_surrogate(1, goodIdx) .* log(pdf_surrogate(1,goodIdx) ./ pdf_M(1,goodIdx)) );
            d2 = 0.5*(d2_1 + d2_2);
            y = d2;
        otherwise %default KL-div
            d1 = sum( pdf_surrogate(1, goodIdx) .* log(pdf_surrogate(1,goodIdx) ./ pdf_origin(1,goodIdx)) );
            d2 = sum( pdf_origin(1,goodIdx) .* log(pdf_origin(1,goodIdx) ./ pdf_surrogate(1,goodIdx)) );
            y = d2;
    end
    
end