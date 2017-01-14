clear;clc;close all;addpath( '../uq_library/ ');

p = 10; case_name = [ 'case_p_', num2str(p), '_JS.mat' ]; load(case_name);

if ~case_origin.bool.multiobject
    x_opt = x_list{1}; fval = fval_list{1};
else
    for opt_step = 1:3;
        x = x_list{opt_step}; fval = fval_list{opt_step}; fval_min = min(fval);
        disp(['min fval = ',num2str(min(fval))]);
    end
    if case_origin.ext.bool.eval == 1
        scale = 2;
        OptIdx = fval(:,1) < scale*fval_min(1) & fval(:,2) < scale*fval_min(2) &fval(:,3) < scale*fval_min(3);
%         OptIdx = fval(:,1) < 05 & fval(:,2) < 1.e-0 &fval(:,3) < 1.e-8;
    else
        OptIdx = fval(:,1) < 0.35 & fval(:,2) < 0.3;
    end

    if max(OptIdx) > 0
        disp('Opt fval');
        OptList = find(OptIdx > 0);
        for i =1:length(OptList)
            disp( ['i= ', num2str(OptList(i)), '  fval= ',num2str(fval(OptList(i), :))] );
        end
    end

    switch p
        case 1
            i_opt = 64;
        case 10
            i_opt = 14;
        otherwise
            i_opt = 1;
    end

    x_opt = x(i_opt, :);
end

r= log(generate_sample(1, 1000, exp(norm(x_opt))));
X = zeros( length(r), length(x_opt) );
for i=1:length(x_opt)
    X(:,i) = exp( r* x_opt(i)/norm(x_opt)*log(case_origin.para.sigma(i)) );
end

[case_surrogate.idt.mean_list,case_surrogate.idt.std_list,case_surrogate.idt.pdf_list] ...
    = post_sample( p, case_origin, X, X1, 'idt' );
[case_surrogate.sl.mean_list,case_surrogate.sl.std_list,case_surrogate.sl.pdf_list] ...
    = post_sample( p, case_origin, X, X1,'sl' );
if case_origin.ext.bool.eval
    [case_surrogate.ext.mean_list,case_surrogate.ext.std_list,case_surrogate.ext.pdf_list] ...
        = post_sample( p, case_origin, X, X1,'ext' );
end