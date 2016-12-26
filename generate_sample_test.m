% m:dimensions
%N: no of samples
clear;clc;
close all;

if nargin == 0
    m = 2;
    N = 1000;
    sigma = [2 2];
end

rng(1);
p = sobolset(m,'Skip',1e2,'Leap',1e1);
% p = scramble(p,'MatousekAffineOwen');
X0 = net(p,N);

for i =1:m
    X0(:,i) = exp(norminv( X0(:,i), 0, log(sigma(i)) )) ; 
end

debug_flag = 1;

if(debug_flag)
    [f, xi] = ksdensity( log(X0(:,1)) );
    figure;
    scatter(log(X0(:,1)) , log(X0(:,2)), 'o');
    figure;
    plot(xi, f);
end