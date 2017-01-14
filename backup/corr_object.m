
x1 = T_list_idt;
y1 = T_list_idt;

for i=1:length(x1)
    y1(i) = corr( IDT_surrogate(:,1),  IDT_surrogate(:,i));
end

figure;
plot(x1, y1, '-o', 'LineWidth',2);
xlabel(' T (K) ');
ylabel(' correlation coeficient ');

%  Scatter plot
s1 = EXT0;
s2 = TEMP_init;
figure;
scatter(s1, s2, 'o');
xlabel('extintion time at 300 K');
ylabel('temperature at 10 x origional extinction time at 300 K');
box on;
title( ['origional para space, corr = ', num2str(corr(s1,s2))] );

% ksdensity plot
figure;
[f, xi] = ksdensity( TEMP_init );
plot( xi, f, '--', 'LineWidth', 2 );
hold all;
xlabel('temperature at 10 x origional extinction time at 300 K');
ylabel('pdf')

% ksdensity plot
figure;
[f, xi] = ksdensity( log10(TEMP_init) );
plot( xi, f, '--', 'LineWidth', 2 );
hold all;
xlabel('log10(temperature) at 10 x origional extinction time at 300 K');
ylabel('pdf')