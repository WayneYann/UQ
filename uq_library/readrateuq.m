function sigma = readrateuq(index, uq_input)
    
    if nargin == 0
        index = 1:3;
    end
    
    data = dlmread(uq_input);
    sigma = data( index, 2 );

end