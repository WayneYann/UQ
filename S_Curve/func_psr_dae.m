function dydt = func_psr_dae( t, y, gas, gslope) %#ok<INUSL>
global yref  Tref  Omega_ref; 
global yin nsp mw;

% Set the state of the gas, based on the current solution vector.
nsp = nSpecies(gas);
setMassFractions(gas, y(2:nsp+1), 'norm');
set(gas,'H',y(nsp+2),'P',pressure(gas));
TEMP = temperature(gas);

xx = [y(2:nsp+1)/yref; TEMP/Tref; log( y(nsp+3)/Omega_ref )];

% Energy equation
rrho = 1.0/density(gas);
wdot = netProdRates(gas);

dydt = zeros(nsp+3, 1);
dydt(1) = 1.;
%Yi
dydt(2:nsp+1) = rrho*mw.*wdot - y(nsp+3)*( y(2:nsp+1)-yin(1:nsp) );
%h
dydt(nsp+2) = - y(nsp+3)*( y(nsp+2)-yin(nsp+1) );
dydt(nsp+3) = y(1) - gslope'*xx;

end