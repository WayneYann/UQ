function [EXT, TEMP_EXT, ext0] = uq_ext(P_in, T_in, COMP,index, X0, destination_ext)
    % ext0 : nominal extinction resident time
    gas = IdealGasMix('H2Reaction_Konnov.xml');   
    Omega_init = 1./(1.E-4);
    
    [ext0, ~, ~] = psr_extinction(gas, P_in, T_in, COMP, Omega_init);
    
    EXT = zeros(length(X0),1);
    TEMP_EXT = zeros(length(X0),1);
    TEMP_init = zeros(length(X0),1);
    Omega_init = 1./(10*ext0); 
    for i = 1:length(X0)
        setMultiplier(gas, index, X0(i,:))
        
        [EXT(i), TEMP_EXT(i), TEMP_init(i), niter] = psr_extinction(gas, P_in, T_in, COMP, Omega_init);
        
        fprintf( '%i\t %i\t %.6e\t %.2f %.2f\n', i, niter-1, EXT(i), TEMP_EXT(i), TEMP_init(i));
        
        setMultiplier(gas, 1.0)
    end
    
    dlmwrite( destination_ext, [EXT, TEMP_EXT, TEMP_init],'precision','%.20f' );
    
end