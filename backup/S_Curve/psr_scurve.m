
clear;
digits(64);
gas = IdealGasMix('H2_Li.cti');

%% PARA
Omega_init = 1./(1.0E-5);
Omega_max = 1e+6;
reso = 2e-1;
nmax = 30000;
T_in =1000.0;
P_in = 1.0;
COMP = 'H2:2,O2:1,N2:3.76';
set(gas,'T',T_in,'P', P_in*oneatm,'X',COMP);

setMultiplier(gas,9,0.5)

global yin nsp mw;
mw = molecularWeights(gas);
nsp = nSpecies(gas);
yin = [massFractions(gas); enthalpy_mass(gas)];

%% EQUILIBRIUM
equilibrate(gas,'HP');
% tad = temperature(gas);
y0 = [massFractions(gas); enthalpy_mass(gas)];

%% Solve initial points
y1 = func_unsteady_psr(gas, y0, Omega_init);
y2 = func_unsteady_psr(gas, y0, Omega_init*1.01);
y3 = func_unsteady_psr(gas, y0, Omega_init*1.02);
% disp([y1(end), y2(end), y3(end)]);

%% Arc-length Continuation
global yref  Tref  Omega_ref; 
yref = 1.; Tref = 1.; Omega_ref = 1.;

OUT = zeros(nmax, nsp+3 );
fprintf( 'niter\t t_res[1/s]\t TEMP[K]\n');
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
%     options = odeset('RelTol',1.e-5,'AbsTol',1.e-9,'Stats','off','mass', M, 'InitialSlope', yprime, 'MaxOrder',2);
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
disp(['CPU time = ' num2str(cputime - t0) ' s']);

if( OUT(niter, end) == 0)
    OUT(niter:end, :) = [];
else
    OUT(niter+1:end, :) = [];
end
t_ext = 1./OUT(max(niter-1,1), 1);
TEMP_ext = OUT(max(niter-1,1), end);
fprintf( '%i\t %.6e\t %.2f\n', niter-1, t_ext, TEMP_ext);
s_curve = [1./OUT(:,1) , OUT(:,end)];
save s_curve;

if 1
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
%     text(x1_txt, y1_txt, txt1);
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