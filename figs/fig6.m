clear; clc;

%%

hf = figure('Unit', 'Centimeters', 'Position', [1 1 8 7]);
mkr_size = 2; lw = 0.75;
font_sz = 7;

% ax2 =  axes('OuterPosition', [0.01, 0.01, 0.35, 0.97]); 
ax4 =  axes('OuterPosition', [-0.01, 0.01, 1.05, 0.98]); 
% ax4 =  axes('OuterPosition', [0.67, 0.01, 0.35, 0.97]);

clist = colororder();

% strong pump limit
% w_L = linspace(0.1,3,1000);
w_L = logspace(-1,2,1000);
cfs_strong = get_characteristic_functions_strongPump(w_L);


loglog(ax4, w_L, cfs_strong.F5, '-', 'DisplayName', '$\mathcal{F}_5(w/L)$'); hold(ax4, 'on'); grid(ax4, 'on');
loglog(ax4, w_L, cfs_strong.F6, '-', 'DisplayName', '$\mathcal{F}_6(w/L)$');
loglog(ax4, w_L, w_L.^(-2)/90, '--', 'DisplayName', '$(w/L)^{-2}/90$', 'Color', 'k');
loglog(ax4, w_L, w_L.^(-6)/6.2e5, '--', 'DisplayName', '$(w/L)^{-6}$/6.2e5', 'Color', clist(5,:));

xlabel(ax4, '$w/L$', 'interpreter', 'latex');
legend(ax4, {},'location','best','interpreter','latex');
xlim(ax4, [0.15, 20]);


% set font size
allAxes = hf.findobj('Type','axes');
allLegend = hf.findobj('Type','legend');
allLine = hf.findobj('Type','line');
allErrorbar = hf.findobj('Type','errorbar');
allTextbox = hf.findobj('Type','textbox');
set([allAxes;allLegend;allTextbox],'fontsize',font_sz);
set([allLine;allErrorbar],'LineWidth',lw, 'MarkerSize', mkr_size);

% annotation('textbox', [0.00 0.94 0.05 0.05], 'String', '(a)', 'FontSize', 8, 'EdgeColor', 'none');
% annotation('textbox', [0.33 0.94 0.05 0.05], 'String', '(b)', 'FontSize', 8, 'EdgeColor', 'none');
% annotation('textbox', [0.66 0.94 0.05 0.05], 'String', '(c)', 'FontSize', 8, 'EdgeColor', 'none');

% Export figure
print(hf, ['Fig6'], '-dpng', '-r600');
hf.PaperSize = hf.Position(3:4);
print(hf, ['Fig6'], '-dpdf', '-painters');
