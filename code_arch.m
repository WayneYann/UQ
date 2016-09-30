
%% Sensitivity
% N_train = nrxn;
% EXT = zeros(N_train, 2);
% for i =1:N_train
%     disp(i);    
%     index = i;
%     factor = 2;
%     EXT(i,:) = psr_extinction( P_in, T_in, COMP, Omega_init, index, factor);
% end
% barh(EXT(:,1) - mean(EXT(:,1)));
% set(gca,'YTick',[1:N_train]) ;
% set(gca, 'YTickLabel', reactionEqn(gas));
% xlabel('Sensitivity');
% title('Exdtinction residence time sensitivity');
% return;
