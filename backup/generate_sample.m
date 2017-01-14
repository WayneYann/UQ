function [X0] = generate_sample(m, N, sigma)
% m:dimensions
%N: no of samples

if nargin == 0
    m = 2;
    N = 300;
    sigma = [exp(1), exp(1)];
end

rng(1);
p = sobolset(m,'Skip',1e2,'Leap',1e1);
% p = scramble(p,'MatousekAffineOwen');
X0 = [];
while length(X0) < N,
    X0 = [X0; net( p, int32(N/(0.9973.^m)*1.2) )];
    X0 = exp(norminv( X0, 0, 1 )) ; 
    X0(any(abs(log(X0'))>3),:) = [];
end
% min(log(X0))
% max(log(X0))
X0 = X0(1:N,:);

for i =1:m
%     X0(:,i) = exp(norminv( X0(:,i), 0, log(sigma(i))./3 )) ; 
    X0(:,i) = exp(log(X0(:,i)).*log(sigma(i))./3 );
end

debug_flag = 0;
if (debug_flag)
    close all;
    [f, xi] = ksdensity( log(X0(:,1)) );
    figure(1);
    scatter(log(X0(:,1)) , log(X0(:,2)), 'o');
    hold all;
    figure(2);
    plot(xi, f);
    hold all;
    plot([log(sigma(i)), log(sigma(i))], [0, max(f)], '--');
    plot([log(sigma(i)), log(sigma(i))]./3.*2, [0, max(f)], '-.');
    plot([log(sigma(i)), log(sigma(i))]./3.*1, [0, max(f)], '-');
    hold all;
end

end