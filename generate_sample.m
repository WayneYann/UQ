function [X0] = generate_sample(m, N, sigma, debug_flag)
% m:dimensions
%N: no of samples

p = sobolset(m,'Skip',1e3,'Leap',1e2);
% p = scramble(p,'MatousekAffineOwen');
X0 = net(p,N);

for i =1:m
    X0(:,i) = exp( norminv( X0(:,i), 0, log(sigma(i)) ) ); 
end

if ~isnan(debug_flag)
    debug_flag = 0;
end

if(debug_flag)
    [f, xi] = ksdensity(X0(:,1));
    figure;
    scatter(log(X0(:,1)) , log(X0(:,2)), 'o');
    figure;
    plot(log(xi), f);
end