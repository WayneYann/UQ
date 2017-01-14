% unsetenv DISPLAY
% matlab >&! matlab.out << EOF
%% uq.m
clear;clc;close all;

if ispc
    % Code to run on Windows platform
    close all;
end

addpath( '../uq_library/ ');
bf.fuel = 'h2';

switch bf.fuel
    case 'h2'
        [bf.mech, uq_input, bf.COMP] = deal( 'H2Reaction_Konnov.xml', 'H2_Uncertainty_Format_LH.txt',...
            'H2:2,O2:1,N2:3.76' );
        bf.nrxn = 33;
    case 'ch4'
        [bf.mech, uq_input, bf.COMP] = deal( 'grimech30.xml', 'GRI_uncertainty.txt',...
            'CH4:1,O2:2,N2:7.52' );
        bf.nrxn = 217;
end

bf.p = 1;
bf.idt.T_list = 1000:400:1800;
bf.sl.T_list = 300;
bf.ext.T_list = 300;

bf.index = 1:bf.nrxn;
bf.sigma = readrateuq(bf.index, uq_input); % UF = 3 sigma
bf.perturbation = 1.1;

for T = bf.idt.T_list
%     bf.idt.log_normal_sensitivity = run_bf_sensitivity( bf, T, 'idt' );
end

T = 300;
% bf.sl.log_normal_sensitivity = run_bf_sensitivity( bf, T, 'sl' );

bf.ext.log_normal_sensitivity = run_bf_sensitivity( bf, T, 'ext' );

movefile( '*.fig', 'fig/.' );
movefile( '*.tiff', 'fig/.' );