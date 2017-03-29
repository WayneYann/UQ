clear;clc;close all;addpath( '../uq_library/');

p = 1; case_name = [ 'case_p_', num2str(p), '_JS.mat' ]; load(case_name);

if ~case_origin.bool.multiobject
    x = x_list{1}; fval = fval_list{1};
else
    for opt_step = 1:3;
        xx = x_list{opt_step}; 
        fval = fval_list{opt_step}; 
        fval_min = min(fval);
        disp(['Step ',num2str(opt_step),' min fval = ',num2str(min(fval))]);
    end
    if case_origin.ext.bool.eval == 1
        scale = 1.5;
        OptIdx = fval(:,1) < scale*fval_min(1) & ...
            fval(:,2) < scale*fval_min(2) & ...
            fval(:,3) < scale*fval_min(3);
    else
        OptIdx = fval(:,1) < 0.35 & fval(:,2) < 0.3;
    end

    if max(OptIdx) > 0
        disp('Optimal fval');
        OptList = find(OptIdx > 0);
        for i =1:length(OptList)
            disp( ['i= ', num2str(OptList(i)), '  fval= ',num2str(fval(OptList(i), :))] );
        end
    end

    switch p
        case 1
            i_opt = 64;
            disp( ['i= ', num2str(i_opt), '  fval= ',num2str(fval(i_opt, :))] );
        case 10
            i_opt = 14;
        otherwise
            i_opt = 1;
    end
    x = xx(i_opt, :);
end
    
    read_active = dlmread('active_eignvector_all.txt')';   
    active_eignvec = read_active([1,2,3,4],:);
    
    A = active_eignvec';
    AAT = A*A';
    w = ones(size(active_eignvec,1),1);
%     [U,S,V] = svd(AAT);
%     r = 1;
%     for i=1:size(S,1)
%         if trace(S(1:i,1:i)) > 0.999*trace(S)
%           r = i;
%           disp(r);
%           break;
%         end
%     end
    
%     v = A'*(AAT\w;
%     v = ( U(:,1:r)*S(1:r,1:r)*V(1:r,:))\(A*w);
%     v = AAT_approx\(A*w);
%     disp('Basis vector v:');
%     disp(v'./norm(v));
%     disp('Lengthe of v:');
%     disp(norm(v));
%     tol=1e-15; 
%     maxit=100;
%     v = bicgstab( AAT, A*w, tol, maxit );
%     UB = ones(size(active_eignvec,2),1);
%     UB = ones(r,1);
%     LB = -UB;
%     options = optimoptions('lsqlin','Algorithm', ...
%         'trust-region-reflective','Display','iter','FunctionTolerance', 1.e-12);
%     Alpha = lsqlin( A*A'*U(:,1:r), A*w, [], [], [], [], LB, UB);
%     v = U(:,1:r)*Alpha
%     
%     A'*v
    v = A*pinv(A'*A)*(w)
    disp('Basis vector v:');
    disp(v'./norm(v));
    disp('Length of v:');
    disp(norm(v));
%     v = v./v
%     v = x';
%     v = read_active(4,:)';

%% Visulaize the basis vector for the subspace
if 0
   figure();
   index = 1:length(x_dir);
   plot(index, x_dir, 'o--','LineWidth',1);
   hold all;
   plot( [1, length(x_dir)],[0,0], '-k', 'LineWidth',1 );
   xlabel('Index');
   ylabel('Component');
   title(['Scaling factor: ', num2str(x_length)]);
   xlim([1,length(x_dir)]);
   ylim([-0.5,0.5]);
   set(gca, 'FontSize',14);
end

%% Project the whole space onto the subspace
if 0
   % Recale the input set to check if the surrogate subspace still work.
   X1 = generate_sample(m, N_estimate, sigma);
end

X1_norm = log(X1)./( ones(size(X1,1), 1)*log(case_origin.para.sigma') );
X = exp( X1_norm*v./norm(v')*v'.*...
    ( ones(size(X1,1), 1)*log(case_origin.para.sigma')));

[case_surrogate.idt.mean_list,case_surrogate.idt.std_list,case_surrogate.idt.pdf_list] ...
    = post_sample( p, case_origin, X, X0, X1, 'idt' );
[case_surrogate.sl.mean_list,case_surrogate.sl.std_list,case_surrogate.sl.pdf_list] ...
    = post_sample( p, case_origin, X, X0, X1,'sl' );
if case_origin.ext.bool.eval
    [case_surrogate.ext.mean_list,case_surrogate.ext.std_list,case_surrogate.ext.pdf_list] ...
        = post_sample( p, case_origin, X, X0, X1,'ext' );
end
post_sample( p, case_origin, X, X0, X1,'diffext' );

h = get(0,'children');
if ~isdir('fig')
    mkdir('fig');
end
for i=1:length(h)
    %stuff where i is the numbering of the figure *and* the handle to use,
    fig = h(i);
    set(findall(fig,'-property','FontSize'),'FontSize',10);
    set(gca, 'FontSize',14);
    saveas(fig,['fig/figure_number_' num2str(fig.Number) '.tiff'])
    saveas(fig,['fig/figure_number_' num2str(fig.Number) '.fig'])
end
% close all;
