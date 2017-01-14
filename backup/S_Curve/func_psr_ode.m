function dydt = func_psr_ode( t, y, Omega, gas) %#ok<INUSL>
global yin nsp mw;

% Set the state of the gas, based on the current solution vector.
setMassFractions(gas, y(1:nsp), 'norm');
set(gas,'H',y(nsp+1),'P',pressure(gas));

% energy equation
rrho = 1.0/density(gas);
wdot = netProdRates(gas);

dydt = zeros(nsp+1, 1);

dydt(1:nsp) = rrho*mw.*wdot - Omega*( y(1:nsp)-yin(1:nsp) );
dydt(nsp+1) = - Omega*( y(nsp+1)-yin(nsp+1) );

end