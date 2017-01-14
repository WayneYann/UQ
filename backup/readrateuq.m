function sigma = readrateuq(index)
    
    if nargin == 0
        index = 1:3;
    end
    
    data = dlmread('H2_Uncertainty_Format_LH.txt');
    sigma = data( index, 2 );

end