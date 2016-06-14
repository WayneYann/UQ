function [t_ext, TEMP_ext] = psr_extinction(P_in, T_in, COMP, Omega_init, index, factor) 
    % ///////////////////////////////////////////////////////////////////
    % PSR: PERFECTLY STIRRED REACTOR
    % 
    %  WRITTEN BY:
    %      ZHUYIN REN (Laniu), July 24, 2013  
    %  Modified by WEIQI JI, Jan 12, 2016
    %  Coding in Cantera-based by WEIQI JI, June 12th, 2016
    % C COMMENTS: 
    % C The current available PSR code developed at Sandia is unable to realiably predict the
    % C unstable branches of the solution. Here a new Fortran code based on arc-length continuation
    % C in the composition(mass fractions)-frequence-temperature space is developed and 
    % C demonstrated. The new method is capable of predicting the combustion curve accurately and    
    % C is able to treat cases with heat loss. 
    % C
    % C The steady-state PSR equations solved are:
    % C 
    % C  0=-\omega ( Y- Y^{in})+ S(Y),    where Y is the vector of mass fractions (KK)
    % C  0=-\omega (h-h^{in})-Q_{LOSS}/(\rho V),
    % C
    % C where: omega is the inverse residence time; Y is the vector of mass fractions; 
    % C        S its rate of change due to chemical reactions; h is specific enthalpy;
    % C        Q_{LOSS} is the heat transfer rate from the reactor; rho is density;
    % C        and V is the volume of the reactor.
    % C
    % C*********************************************************
    % C Method:
    % C  Arc-length continuation.  The variables Y, h, T, omega are considered as functions
    % C  of the arclength, s.  From a consistent initial condition, the equations solved are:
    % C 
    % C  0=-\omega ( Y- Y^{in})+ S(Y),    where Y is the vector of mass fractions (KK)
    % C  0=-\omega (h-h^{in})-Q_{LOSS}/(\rho V),
    % C  dZ/ds=1,  
    % C  0=Z- g'*X,
    % C  where X=[ Y/Y^{ref}, T/T^{ref}, \ln(\omega/\omega^{ref})] and  g=(X^k- X^{k-1})/(s_k- s_{k-1}),
    % C  where k and k-1 denote the previous two computed points.
    % C The default value for the reference values are 1. And we employ ddasac to solve the above DAEs.
    % /////////////////////////////////////////////////////////////////////////////////////

    gas = IdealGasMix('H2Reaction_Konnov.xml');   
    if nargin == 0
        COMP = 'H2:2,O2:1,N2:3.76';
        P_in = 1.0;
        T_in = 700.0;
        Omega_init = 1./(1E-5);
        index = [13,  29];
        factor = [1 1];
    end
    %% PARA
%     Omega_init = 1./(2E-5);
    Omega_max = 1e+8;
    reso = 2e-1;
    nmax = 30000;
    
%     T_in =300.0;
%     P_in = 1.0;
%     COMP = 'H2:2,O2:1,N2:3.76';
    set(gas,'T',T_in,'P', P_in*oneatm,'X',COMP);

    global yin nsp mw;
    mw = molecularWeights(gas);
    nsp = nSpecies(gas);
    yin = [massFractions(gas); enthalpy_mass(gas)];
    
    setMultiplier(gas, index, factor)
    
    %% EQUILIBRIUM
    equilibrate(gas,'HP');
    % tad = temperature(gas);
    y0 = [massFractions(gas); enthalpy_mass(gas)];

    %% Solve initial points
    y1 = func_unsteady_psr(gas, y0, Omega_init);
    y2 = func_unsteady_psr(gas, y0, Omega_init*1.01);
    y3 = func_unsteady_psr(gas, y0, Omega_init*1.02);

    %% Arc-length Continuation
    global yref  Tref  Omega_ref; 
    yref = 1.; Tref = 1.; Omega_ref = 1.;

    OUT = zeros(nmax, nsp+3 );
    % fprintf( 'niter\t t_res[1/s]\t TEMP[K]\n');
    niter = 0; 
    t0 = cputime;

    while( niter < nmax && y3(1) < Omega_max)
        niter = niter +1;
        xx1 = [y1(2:nsp+1)/yref; y1(nsp+3)/Tref; log(y1(1)/Omega_ref) ];
        xx2 = [y2(2:nsp+1)/yref; y2(nsp+3)/Tref; log(y2(1)/Omega_ref) ];
        xx3 = [y3(2:nsp+1)/yref; y3(nsp+3)/Tref; log(y3(1)/Omega_ref) ];

        xx_dif12 = xx2-xx1;
        xx_dif23 = xx3-xx2;

        s12 = sqrt( xx_dif12'*xx_dif12 );
        s23 = sqrt( xx_dif23'*xx_dif23 );

        gslope = xx_dif23/s23;
        curv = xx_dif23/s23 - xx_dif12/s12;
        curvature = sqrt( curv'*curv )/s23;
        ds = reso*sqrt( gslope'*gslope ) / max( curvature, 8e-2 ) ;
        ds = min(ds, 10*s23);
        ds = max(ds, 1e-6);

        ds = ds*exp(-1);
        s = 0.; sout = ds;

        y0 = y3;
        y0(nsp+3) = y3(1);
        y0(1) = gslope'*xx3;

        yprime = 0*y0;
        yprime(1) = 1.0;
        M = zeros( nsp+3 );
        M(1) = 1;
        % Disable the tolerence setting if it is hard to converge
        %options = odeset('RelTol',1.e-5,'AbsTol',1.e-12,'Stats','off','mass', M, 'InitialSlope', yprime);
        options = odeset('mass', M, 'InitialSlope', yprime);
        tel = [s sout];
        out = ode15s(@func_psr_dae, tel, y0, options, gas, gslope);
        y = out.y(:,end);
        y = [y(nsp+3,end);  y(2:nsp+2,end);  temphy(gas, y(2:nsp+2))];
        OUT(niter, :) = y';
        y1 = y2;
        y2 = y3;
        y3 = y;
    %     fprintf( '%i\t%.6e\t%.2f\n', niter, 1./y3(1), y3(end));

        if( OUT(niter, 1) < OUT(max(niter-1,1), 1))
            break;
        end
    end
%     disp(['CPU time = ' num2str(cputime - t0) ' s']);

    OUT(niter+1:end, :) = [];
    t_ext = 1./OUT(max(niter-1,1), 1);
    TEMP_ext = OUT(max(niter-1,1), end);
    fprintf( '%i\t %.6e\t %.2f\n', niter-1, t_ext, TEMP_ext);

    if 0
        %plot temperature verus residence time
        figure(1);
        h = semilogx( 1./OUT(:,1) , OUT(:,end), '-'  );
        set(h, 'LineWidth', 1);
        xlabel('\it t_{res} (s)');
        ylabel( '\it T_{max} (K)' );
        title([COMP, ' T=', int2str(T_in), 'K',' P=',int2str(P_in),'atm' ])
        hold all;
        txt1 = ['t_{ext} = ', num2str(t_ext), ' [s] @', num2str( OUT(end, end) ), ' K'];
        x1_txt = exp( -(log(max(OUT(:,1))) + log(min(OUT(:,1))))./2 );
        y1_txt = (max(OUT(:,end)) + min(OUT(:,end)))/2;
        text(x1_txt, y1_txt, txt1);
        hold all;
        %plot temperature verus Omega (Mixing frequency)
        %     figure;
        %     h = semilogx( OUT(:,1) , OUT(:,end), 'o-'  );
        %     xlabel('\it Omega (1/s)');
        %     ylabel( '\it T_{max} (K)' );
        %     title([COMP, ' T=', int2str(T_in), 'K',' P=',int2str(P_in),'atm' ])
        %     set(h, 'LineWidth', 1);
        %     hold all;
    end
end

function temp = temphy(gas, y)
    % return temperatrue given H and massfraction
    %y(1:nsp+1)
    setMassFractions(gas, y(1:end-1), 'nonorm');
    set(gas,'H',y(end),'P',pressure(gas));
    temp = temperature(gas);
end

function y = func_unsteady_psr(gas, y0, Omega)
    global nsp;

    options = odeset('RelTol',1.e-5,'AbsTol',1.e-12,'Stats','off');
    tel = [0 20/Omega];
    out = ode15s(@func_psr_ode,tel,y0,options, Omega, gas); %pass extra paras
    y = out.y(:,end);
    y = [Omega; y; temphy(gas, y(1:nsp+1))];
end

function dydt = func_psr_ode( t, y, Omega, gas) %#ok<INUSL>
    global yin nsp mw;
    % Set the state of the gas, based on the current solution vector.
    setMassFractions(gas, y(1:nsp), 'nonorm');
    set(gas,'H',y(nsp+1),'P',pressure(gas));
    % energy equation
    rrho = 1.0/density(gas);
    wdot = netProdRates(gas);
    dydt = zeros(nsp+1, 1);
    dydt(1:nsp) = rrho*mw.*wdot - Omega*( y(1:nsp)-yin(1:nsp) );
    dydt(nsp+1) = - Omega*( y(nsp+1)-yin(nsp+1) );
end

function dydt = func_psr_dae( t, y, gas, gslope) %#ok<INUSL>
    global yref  Tref  Omega_ref; 
    global yin nsp mw;
    % Set the state of the gas, based on the current solution vector.
    nsp = nSpecies(gas);
    setMassFractions(gas, y(2:nsp+1), 'nonorm');
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
