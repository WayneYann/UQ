function IDT = idt_hp_uq(P, T, COMP, X, t_end)
%  IGNITE_HP   Solves the same ignition problem as 'ignite', but uses function
%  conhp instead of reactor. 
%
% tic;
if nargin == 0
    COMP = 'H2:2,O2:1,N2:3.76';
    P = 1;
    T = 1234.03;
    index = [1:nrxn];
    m = length(index); 
    sigma = readrateuq(index);
    %      Generate one-d samples
    r = generate_sample(1, 1000, 2);
    sigma_r = sigma' .* x;
    X = r*sigma_r/2;
    t_end = 1;
end

[rows, cols] = size(X);
IDT = zeros(rows,1);

for i = 1:rows
    gas = IdealGasMix('H2Reaction_Konnov.xml');
    mw = molecularWeights(gas);
    disp(i);
    for j = 1:2
        setMultiplier(gas, j, 2 );
    end
    set(gas,'T',T,'P', P*oneatm,'X',COMP);

    y0 = [temperature(gas); massFractions(gas)];
    tel = [0 t_end];
    options = odeset('RelTol',1.e-8,'AbsTol',1.e-12,'Stats','off');
    out = ode15s(@conhp,tel,y0,options,gas,mw);
    [~, idt_index] = max( diff(out.y(1,:))./diff(out.x) );
    IDT(i) = log10( out.x(idt_index) );
    
    % disp( ['ignition @ ', num2str(idt), ' s'] );

    if out.y(1,idt_index) < T+100
        disp( ['no ignition @Tin=', num2str(T)] );
    end

    if 0
       % plot the temperature and OH mole fractions.
       figure(1);
       plot(out.x,out.y(1,:));
       xlabel('time');
       ylabel('Temperature');
       title(['Final T = ' num2str(out.y(1,end)) ' K']);

       figure(2);
       ioh = speciesIndex(gas,'OH');
       plot(out.x,out.y(1+ioh,:));
       xlabel('time');
       ylabel('Mass Fraction');
       title('OH Mass Fraction');   
    end
    
end
%     toc;
end

function dydt = conhp(t, y, gas, mw) %#ok<INUSL>
% CONHP ODE system for a constant-pressure, adiabatic reactor.
%
%    Function CONHP evaluates the system of ordinary differential
%    equations for an adiabatic, constant-pressure,
%    zero-dimensional reactor. It assumes that the 'gas' object
%    represents a reacting ideal gas mixture.


% Set the state of the gas, based on the current solution vector.
setMassFractions(gas, y(2:end), 'nonorm');
set(gas, 'T', y(1), 'P', pressure(gas));
nsp = nSpecies(gas);

% energy equation
wdot = netProdRates(gas);
tdot = - temperature(gas) * gasconstant * enthalpies_RT(gas)' ...
    * wdot / (density(gas)*cp_mass(gas));

% set up column vector for dydt
dydt = [ tdot
    zeros(nsp, 1) ];

% species equations
rrho = 1.0/density(gas);
for i = 1:nsp
    dydt(i+1) = rrho*mw(i)*wdot(i);
end
end