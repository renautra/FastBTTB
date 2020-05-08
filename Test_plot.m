% Plot results for finding files in the paper
% A Tutorial and Open Source Software for the Efficient Evaluation of Gravity and Magnetic Kernels (2019)
% Jarom Hogue, Rosemary Renaut and Saeed Vatankhah
% Trademarks: 
% Rosemary Renaut and Jarom Hogue (TM)
% Test_plot Version 1: December 13, 2019
%%
% Padding needs to be set to determine which file to load
padding=1;
load(['Results',int2str(padding)]')
ll = {'b*','bh','bx','bp','s','d'}; % linestyles
elem = ((lowerscale:upperscale).*25).^2.*((lowerscale:upperscale).*15).^2.*((lowerscale:upperscale).*2); % number of elements for each size: m*n
upperind = upperscale-lowerscale+1; indmax = scalemax-lowerscale+1; % upper indecies 
%% Gravity and Magnetic Kernels timings
figure(5)
A = zeros(totalruns,upperind); B = zeros(totalruns,indmax); 
A(:,:) = timeG_bttb(1,:,:);                     % T_gravity
B(:,:) = timeG(1,:,1:indmax);                   % G_gravity
loglog(elem(1:indmax),mean(B),ll{1},elem,mean(A),ll{2}) % plot gravity results on log scale
ylabel('time (s)','Interpreter','latex')
A(:,:) = timeG_bttb(2,:,:);                     % T_magnetic
B(:,:) = timeG(2,:,1:indmax);                   % G_magnetic
hold on
loglog(elem(1:indmax),mean(B),ll{3},elem,mean(A),ll{4}) % plot magnetic results on log scale
set(gca, 'YScale', 'log')                       % enforce log scale for y axis
axis tight
legend('$\mathbf G_{gravity}$','$\mathbf{\hat{T}}_{gravity}$','$\mathbf G_{magnetic}$','$\mathbf{\hat{T}}_{magnetic}$','interpreter','latex','location','NW')
figure_properties
if padding>1 print('-depsc',['figure5b','.eps']);else print('-depsc',['figure5a','.eps']);end
%% Gravity and Magnetic forward multiplication timings with errors
figure(6)
A = zeros(totalruns,upperind); B = zeros(totalruns,indmax); 
A(:,:) = timeGx_bttb(1,:,:);                    % T_gravity
B(:,:) = timeGx(1,:,1:indmax);                  % G_gravity
yyaxis left, loglog(elem(1:indmax),mean(B),ll{1},elem,mean(A),ll{2}), ylabel('time (s)') % plot using left y axis
A = zeros(totalruns,upperind); B = zeros(totalruns,indmax); 
A(:,:) = timeGx_bttb(2,:,:);                    % T_magnetic
B(:,:) = timeGx(2,:,1:indmax);                  % G_magnetic
hold on, loglog(elem(1:indmax),mean(B),ll{3},elem,mean(A),ll{4}) % add magnetic timing results to plot
A = zeros(totalruns,indmax); B = zeros(totalruns,indmax); 
A(:,:) = error_bttb(2,:,1:indmax);              % E_magnetic
B(:,:) = error_bttb(1,:,1:indmax);              % E_gravity
yyaxis right, loglog(elem(1:indmax),mean(B),ll{5},elem(1:indmax),mean(A),ll{6}) 
set(gca, 'YScale', 'log') % plot errors on right y axis
legend('$\mathbf G_{gravity}$','$\mathbf{\hat{T}}_{gravity}$','$\mathbf G_{magnetic}$','$\mathbf{\hat{T}}_{magnetic}$','$\mathbf E_{gravity}$','$\mathbf E_{magnetic}$','interpreter','latex','location','NW')
ylabel('mean error')
figure_properties
if padding>1 print('-depsc',['figure6b','.eps']);else print('-depsc',['figure6a','.eps']);end
%% Gravity and Magnetic transpose multiplication with errors
figure(7)
A = zeros(totalruns,upperind); B = zeros(totalruns,indmax); 
A(:,:) = timeGTx_bttb(1,:,:);                   % T_gravity
B(:,:) = timeGTx(1,:,1:indmax);                 % G_gravity
yyaxis left, loglog(elem(1:indmax),mean(B),ll{1},elem,mean(A),ll{2}), ylabel('time (s)') % plot using left y axis
A(:,:) = timeGTx_bttb(2,:,:);                   % T_magnetic
B(:,:) = timeGTx(2,:,1:indmax);                 % G_magnetic
hold on, loglog(elem(1:indmax),mean(B),ll{3},elem,mean(A),ll{4}) % add magnetic timing results to plot
A = zeros(totalruns,indmax); B = zeros(totalruns,indmax); 
A(:,:) = error_transpose_bttb(2,:,1:indmax);    % E_magnetic
B(:,:) = error_transpose_bttb(1,:,1:indmax);    % E_gravity
yyaxis right, loglog(elem(1:indmax),mean(B),ll{5},elem(1:indmax),mean(A),ll{6})
set(gca, 'YScale', 'log'),                      % plot errors on right y axis
legend('$\mathbf G_{gravity}$','$\mathbf{\hat{T}}_{gravity}$','$\mathbf G_{magnetic}$','$\mathbf{\hat{T}}_{magnetic}$','$\mathbf E_{gravity}$','$\mathbf E_{magnetic}$','interpreter','latex','location','NW')
ylabel('mean error')
figure_properties
if padding>1 print('-depsc',['figure7b','.eps']);else print('-depsc',['figure7a','.eps']);end


