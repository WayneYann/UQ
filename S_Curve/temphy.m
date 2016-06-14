function temp = temphy(gas, y)
% return temperatrue given H and massfraction
%y(1:nsp+1)
    setMassFractions(gas, y(1:end-1), 'nonorm');
    set(gas,'H',y(end),'P',pressure(gas));
    temp = temperature(gas);
end