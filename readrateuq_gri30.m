function sigma = readrateuq_gri30(index)
    
    if nargin == 0
        index = 1:3;
    end
    
    data = dlmread('GRI_uncertainty.txt');
    sigma = data( index, 2 );

end