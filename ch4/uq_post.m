clear;clc;close all;
p = 1;
case_name = [ 'case_p_', num2str(p), '_TVD.mat' ];  load(case_name);

if ~case_origin.bool.multiobject
    x = x_list{1};
    fval = fval_list{1};
else
    for opt_step = 1:3;
        x = x_list{opt_step};
        fval = fval_list{opt_step};
        fval_min = min(fval);
    disp(['min fval = ',num2str(min(fval))]);
    end
    if case_origin.ext.bool.eval == 1
        scale = 4;
        OptIdx = fval(:,1) < scale*fval_min(1) & fval(:,2) < scale*fval_min(2) &fval(:,3) < scale*fval_min(3);
        OptIdx = fval(:,1) < 0.1 & fval(:,2) < 0.03 &fval(:,3) < 0.05;
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
            i_opt = 90;
        case 10
            i_opt = 83;
        case 30
            i_opt = 63;
        otherwise
            i_opt = 1;
    end
    x_opt = x(i_opt, :);
end
rng default % for reproducibility

XTS = X0;
for i = 1:length(x_opt)
    XTS(:,i) = log(X0(:,i))./log(case_origin.para.sigma(i));
end
r = XTS * x_opt';
X = zeros( length(r), length(x_opt) );
for i=1:length(x_opt)
    X(:,i) = exp( r* x_opt(i)/norm(x_opt)*log(case_origin.para.sigma(i)) );
end

%case_origin.idt.T_list = 1000:400:1800;
[case_surrogate.idt.mean_list,case_surrogate.idt.std_list,case_surrogate.idt.pdf_list] ...
    = post_sample( p, case_origin, X, X1, 'idt' );
[case_surrogate.sl.mean_list,case_surrogate.sl.std_list,case_surrogate.sl.pdf_list] ...
    = post_sample( p, case_origin, X, X1,'sl' );
if case_origin.ext.bool.eval
    [case_surrogate.ext.mean_list,case_surrogate.ext.std_list,case_surrogate.ext.pdf_list] ...
        = post_sample( p, case_origin, X, X1,'ext' );
end