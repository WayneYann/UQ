function log_normal_sensitivity =  run_bf_sensitivity( bf, T, type )

gas = IdealGasMix(bf.mech);

%% 
disp(['Evaluate ', type  ,' samples']);

source_out = fullfile('data', ['samples_out_',type,'.txt']);   %ignition delay time
destination_out = fullfile( 'data', [bf.fuel,'_','samples_out_',type,'_', num2str(bf.p),'_', num2str(T),'_K.txt' ]);

switch type
    case 'idt'
        system( ['python bf_sensitivity_idt.py '  num2str(bf.p), ' ', num2str(T),...
            ' ', bf.COMP, ' ',bf.mech, ' ', num2str(bf.perturbation),' ', num2str(length(bf.index)) ]);
    case 'sl'
        system( ['python bf_sensitivity_SL.py '  num2str(bf.p), ' ', num2str(T),...
            ' ', bf.COMP, ' ',bf.mech, ' ', num2str(bf.perturbation),' ', num2str(length(bf.index)) ]);
    case 'ext'
        [~,~,~] = bf_ext(bf, T, source_out);
end
copyfile( source_out, destination_out);

normal_ratio = dlmread(destination_out);
log_normal_sensitivity = log(normal_ratio)./log( bf.perturbation );

abs_sens = abs( log_normal_sensitivity );

[~, I_sort] = sort( abs_sens, 'descend' );
cutoff = 10;

rxn_equations = reactionEqn(gas, I_sort(1:cutoff));

%% Sensitivity
figure();
barh( log_normal_sensitivity( I_sort(1:cutoff) ));
% Set the axis limits
axis([-max( abs_sens ), max(abs_sens), 0, cutoff+0.5]);

set(gca, 'YTick', 1:cutoff);
set(gca, 'YTickLabel', strcat(rxn_equations,' (', num2str(I_sort(1:cutoff)), ' )' ));
xlabel( 'log normal sensitivity' );

title( [ num2str(bf.p), ' atm ', num2str(T), ' K ', bf.COMP, ' (', type, ' )']);
fname = [bf.fuel,'_',num2str(bf.p), '_atm_', num2str(T), '_K_',type, '_log_normal_sens.fig'];
savefig( fname );
fname = [bf.fuel,'_',num2str(bf.p), '_atm_', num2str(T), '_K_',type, '_log_normal_sens.tiff'];
saveas( gcf, fname );

%% UQ estimated by local sensitivity
UQ = log(normal_ratio)./(bf.perturbation-1).*(bf.sigma-1);
abs_UQ = abs(UQ);
[~, I_sort] = sort( abs_UQ, 'descend' );
rxn_equations = reactionEqn(gas, I_sort(1:cutoff));
figure();
barh( UQ( I_sort(1:cutoff)) );
% Set the axis limits
% axis([-max( abs_UQ ), max( abs_UQ), 0, cutoff+0.5]);

set(gca, 'YTick', 1:cutoff);
set(gca, 'YTickLabel', strcat(rxn_equations,' (', num2str(I_sort(1:cutoff)), ' )' ));
xlabel( 'log normalized UQ' );

title( [num2str(bf.p), ' atm ', num2str(T), ' K ', bf.COMP, ' (', type, ' )']);
fname = [bf.fuel,'_',num2str(bf.p), '_atm_', num2str(T), '_K_',type, '_normal_UQ.fig'];
savefig( fname );
fname = [bf.fuel,'_',num2str(bf.p), '_atm_', num2str(T), '_K_',type, '_normal_UQ.tiff'];
saveas( gcf, fname );

end