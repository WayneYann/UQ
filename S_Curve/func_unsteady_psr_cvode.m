function y = func_unsteady_psr_cvode(gas, y0, Omega)
global nsp;

options = odeset('RelTol',1.e-12,'AbsTol',1.e-20,'Stats','off');
tel = [0 20/Omega];
out = ode15s(@func_psr_ode,tel,y0,options, Omega, gas); %pass extra paras
y = out.y(:,end);
y = [Omega; y; temphy(gas, y(1:nsp+1))];

end