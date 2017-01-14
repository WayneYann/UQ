
clear;clc;
%generate extinction time versus extinction temperature curve
P_in = 1;
T_in_list = 300:100:1000;
COMP = 'H2:2,O2:1,N2:3.76';
Npoints = length(T_in_list);
OUT = zeros( Npoints, 4 );

Omega_init = 1./(2E-5);
   
for i=1:length(T_in_list)
    T_in = T_in_list(i);
    [t_ext, TEMP_ext] = psr_extinction(P_in, T_in, COMP, Omega_init);
    idt = idt_hp(P_in, TEMP_ext, COMP);
    OUT(i,1:4) = [T_in, TEMP_ext, t_ext, idt];
    Omega_init = (1./t_ext)./2;
    fprintf( '%i\t %.2f\t %.6e\t %.6e\n', T_in, TEMP_ext,  t_ext, idt);
end

figure(1);
h1 = semilogy( 1000./OUT(:,2) , OUT(:,3), 'o-' ,'DisplayName','t_{ext}' );
hold all;
h2 = semilogy( 1000./OUT(:,2) , OUT(:,4), '^-', 'DisplayName','idt_{ext}' );
legend('show');
xlabel('1000/T_{ext} (K)');
ylabel('Time (s)');
hold all;

% save OUT;

