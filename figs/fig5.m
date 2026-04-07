clear; clc;

%%

hf = figure('Unit', 'Centimeters', 'Position', [1 1 18 6]);
mkr_size = 2; lw = 0.75;
font_sz = 7;

ax2 =  axes('OuterPosition', [0.01, 0.01, 0.35, 0.97]); 
ax3 =  axes('OuterPosition', [0.34, 0.01, 0.35, 0.97]); 
ax5 =  axes('OuterPosition', [0.67, 0.01, 0.35, 0.97]);

clist = colororder();

% weak pump limit
w_L = linspace(0,2.0,1000);
OD = [1, 2.25, 3.75, 6, 10];
dOD = 0.01;
cfs_weak = get_characteristic_functions_weakPump(w_L, OD);
cfs_weak_dOD = get_characteristic_functions_weakPump([1], [OD - dOD;  OD + dOD]);    % 计算 dG / dOD
nameList = {'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8'};
for ii = 1:length(nameList)
    nstr = nameList{ii};
    cfs_weak_dOD.(nstr) = ( cfs_weak_dOD.(nstr)(2,:) - cfs_weak_dOD.(nstr)(1,:) )/(2*dOD);    % 中心差分公式求导数
end

% 1st order correction
line_spec = {'-','-','-','--','--','--'};
for ii = 1:length(OD)
    plot(ax2, w_L, cfs_weak.F1*cfs_weak.G1(ii) + cfs_weak.F2*cfs_weak.G2(ii), line_spec{ii}, 'DisplayName', sprintf('%.3g', OD(ii)), 'color', clist(ii,:));  
    hold(ax2, 'on'); grid(ax2, 'on');
    C1 = cfs_weak.F1*cfs_weak.G1(ii) + cfs_weak.F2*cfs_weak.G2(ii);
    dC1_dtau = cfs_weak.F1*cfs_weak_dOD.G1(ii) + cfs_weak.F2*cfs_weak_dOD.G2(ii);
%     plot(ax2, w_L, C1 + OD(ii)*dC1_dtau, '-.', 'DisplayName', sprintf('%g', OD(ii)), 'color', clist(ii,:));  
end
xlabel(ax2, '$w/L$', 'interpreter', 'latex');
ylabel(ax2, '$\mathcal{C}_1$', 'interpreter', 'latex');
lg = legend(ax2, {},'location','best','interpreter','latex');
title(lg, 'OD');

% 1st order correction  温度系数
for ii = 1:length(OD)
    C1 = cfs_weak.F1*cfs_weak.G1(ii) + cfs_weak.F2*cfs_weak.G2(ii);
    dC1_dtau = cfs_weak.F1*cfs_weak_dOD.G1(ii) + cfs_weak.F2*cfs_weak_dOD.G2(ii);
    plot(ax3, w_L, C1 + OD(ii)*dC1_dtau, line_spec{ii}, 'DisplayName', sprintf('%g', OD(ii)), 'color', clist(ii,:));  
    hold(ax3, 'on'); grid(ax3, 'on');
end
xlabel(ax3, '$w/L$', 'interpreter', 'latex');
ylabel(ax3, '$\mathcal{C}_1 + \tau \cdot \partial \mathcal{C}_1 / \partial \tau$', 'interpreter', 'latex');
% lg = legend(ax3, {},'location','best','interpreter','latex');
% title(lg, 'OD');

% 3rd order correction
for ii = 1:length(OD)
    plot(ax5, w_L, cfs_weak.F3*cfs_weak.G7(ii) + cfs_weak.F4*cfs_weak.G8(ii), line_spec{ii}, 'DisplayName', sprintf('%g', OD(ii)));  
    hold(ax5, 'on'); grid(ax5, 'on');
end
xlabel(ax5, '$w/L$', 'interpreter', 'latex');
ylabel(ax5, '$\mathcal{F}_3\mathcal{G}_7 + \mathcal{F}_4\mathcal{G}_8$', 'interpreter', 'latex');
% lg = legend(ax5, {},'location','best','interpreter','latex');
% title(lg, 'OD');



% set font size
allAxes = hf.findobj('Type','axes');
allLegend = hf.findobj('Type','legend');
allLine = hf.findobj('Type','line');
allErrorbar = hf.findobj('Type','errorbar');
allTextbox = hf.findobj('Type','textbox');
set([allAxes;allLegend;allTextbox],'fontsize',font_sz);
set([allLine;allErrorbar],'LineWidth',lw, 'MarkerSize', mkr_size);

annotation('textbox', [0.00 0.94 0.05 0.05], 'String', '(a)', 'FontSize', 8, 'EdgeColor', 'none');
annotation('textbox', [0.33 0.94 0.05 0.05], 'String', '(b)', 'FontSize', 8, 'EdgeColor', 'none');
annotation('textbox', [0.66 0.94 0.05 0.05], 'String', '(c)', 'FontSize', 8, 'EdgeColor', 'none');

% Export figure
print(hf, ['Fig5'], '-dpng', '-r600');
hf.PaperSize = hf.Position(3:4);
print(hf, ['Fig5'], '-dpdf', '-painters');
