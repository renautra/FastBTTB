%% Short script used for generating figures consistently with figures 5,6,7 in
% A Tutorial and Open Source Software for the Efficient Evaluation of Gravity and Magnetic Kernels (2019)
% Jarom Hogue, Rosemary Renaut and Saeed Vatankhah
% Trademarks: 
% Rosemary Renaut and Jarom Hogue (TM)
% figure_properties Version 1: December 13, 2019. Rosemary Renaut and Jarom Hogue
%%
xticks([25^2*15^2*2 50^2*30^2*4 100^2*60^2*8 175^2*105^2*14 300^2*180^2*24]) % locations for x ticks
xticklabels({'(25,15,2)','(50,30,4)','(100,60,8)','(175,105,14)','(300,180,24)'}) % display for x ticks
xlabel('$(s_{x},s_{y},n_{z})$','Interpreter','latex')
axislims=axis;
axis([25^2*15^2*2  300^2*180^2*24 axislims(3) axislims(4)])
set(gca,'TickLabelInterpreter', 'latex');
grid on
set(findobj('Type','line'),'MarkerSize',12) % Update formatting of the figure for readability
set(findobj('Type','line'),'LineWidth',2) 
set(findobj('Type','axes'),'FontSize',16)
set(findobj('Type','text'),'FontWeight','bold')
set(findobj('Type','axes'),'FontWeight','bold')
set(findobj('Type','title'),'FontWeight','bold')
set(findobj('Type','axes'),'FontName','Courier')