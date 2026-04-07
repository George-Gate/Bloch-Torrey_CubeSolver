clear; clc;

calcSzDist;   % get sol_alpha1D and sol_BA_2D




%% plot
% % generate colormap for intensity plot
% clims = [0 20];
% crt1 = 0.5;  crt2 = 2;
% cval = linspace(clims(1),clims(2),2000)';
% cmp = [(cval-clims(1))/(crt1-clims(1)), (cval-crt1)/(crt2-crt1), (cval-crt2)/(clims(2)-crt2)];
% cmp = min(1,max(0,cmp));

hf = figure('Unit', 'Centimeters', 'Position', [1 1 18 7]);
mkr_size = 2; lw = 0.75;
font_sz = 7;
x1 = [0.08, 0.52]; y1 = [0.01, 0.50]; 
ax(1) =  axes('OuterPosition', [-0.01, 0.01, 0.55, 1.00]); 
% ax(2) =  axes('OuterPosition', [0.51, 0.50, 0.45, 0.47]); 
% ax(3) =  axes('OuterPosition', [0.51, 0.04, 0.45, 0.47]); 
ax(2) =  axes('InnerPosition', [0.55, 0.55, 0.37, 0.39]); 
ax(3) =  axes('InnerPosition', [0.55, 0.12, 0.37, 0.39]); 

% 1D plot of photon flux
lineID = [length(sol_alpha1D.alpha0_list):-1:1];
plot(ax(1), sol_alpha1D.z_L_list, sol_alpha1D.alpha_z(lineID,:)./sol_alpha1D.alpha0_list(lineID));  hold(ax(1), 'on');
grid(ax(1), 'on');
plot(ax(1), sol_alpha1D.z_L_list, exp(-OD*(0.5+sol_alpha1D.z_L_list)), '--k')
xlabel(ax(1), '$z/L$', 'interpreter', 'latex');
ylabel(ax(1), '$\Phi(z)/\Phi_0$', 'interpreter', 'latex');
lg = legend(ax(1), num2str([sol_alpha1D.alpha0_list(lineID); 0]), 'Position',[0.371118573698493 0.293798521101819 0.0941176464890732 0.253529405225726]);
title(lg, '$\Phi_0 / \Phi^*$', 'interpreter', 'latex');
% alpha0 = sol_alpha1D.alpha0_list(lineID);
% KK = OD.*(0.5+sol_alpha1D.z_L_list);
% plot(ax(1), sol_alpha1D.z_L_list, 1 - 1./(1+alpha0).*KK, '--');  % 强光极限线性近似
% plot(ax(1), sol_alpha1D.z_L_list, 1 - 1./(1+alpha0).*KK + 1./2./(1+alpha0).^3.*KK.*KK, '-.');  % 强光极限二阶近似

% alpha0 --> Phi0/Phi*

Pid = [1, 4];
for iP = 1:length(Pid)
    cax = ax(iP+1);
    da = sol_BA_2D(Pid(iP));
    sObj = pcolor(cax, da.z_L, da.x_L, abs(da.B_A));   hold(cax, 'on');
    set(sObj,  'FaceColor', 'interp', 'EdgeColor', 'none');
    contour(cax, da.z_L, da.x_L, abs(da.B_A), 'LineColor', 0.4*[1 1 1]);
    axes(cax);    % make ax become gca
    crameri('-acton');
    if iP==2
        xlabel(cax, '$z/L$', 'interpreter', 'latex');
    end
    ylabel(cax, '$\sqrt{x^2+y^2}/L$', 'interpreter', 'latex');
    box(cax, 'on'); grid(cax, 'on');
%     title(cax, ['\alpha_0(0) = ' num2str(da.max_alpha0)]);
    fprintf('max(alpha_0): %g\n', da.max_alpha0);
    set(cax,'Layer','top');
    caxis(cax, [0,45]);
end

set(ax(2), 'XTickLabel', {});

% show color bar
colorbar(ax(3),'Position', [0.94 0.10 0.02 0.80]);
annotation(hf,'textbox',  [0.88 0.96 0.16 0.03], 'String',{'$B_\mathrm{A}$ (nT)'}, 'HorizontalAlignment','center', 'FontSize', font_sz, 'EdgeColor','none','interpreter','latex');

% set font size
allAxes = hf.findobj('Type','axes');
allLegend = hf.findobj('Type','legend');
allLine = hf.findobj('Type','line');
allErrorbar = hf.findobj('Type','errorbar');
allTextbox = hf.findobj('Type','textbox');
set([allAxes;allLegend;allTextbox],'fontsize',font_sz);
set([allLine;allErrorbar],'LineWidth',lw, 'MarkerSize', mkr_size);

annotation('textbox', [0.01 0.93 0.05 0.05], 'String', '(a)', 'FontSize', 8, 'EdgeColor', 'none');
annotation('textbox', [0.50 0.93 0.05 0.05], 'String', '(b)', 'FontSize', 8, 'EdgeColor', 'none');
annotation('textbox', [0.50 0.49 0.05 0.05], 'String', '(c)', 'FontSize', 8, 'EdgeColor', 'none');


% Export figure
print(hf, ['Fig1'], '-dpng', '-r600');
hf.PaperSize = hf.Position(3:4);
print(hf, ['Fig1'], '-dpdf', '-painters');

