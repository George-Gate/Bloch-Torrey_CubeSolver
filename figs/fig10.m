clear; clc;

spp = load(['dataSet' filesep 'sweepPumpPower' filesep 'sweepPumpPower 20250626_170731.mat']);   % numBasis = 100, Rrel=1/(100us), w/L=0.4
sqg = load(['dataSet' filesep 'sweepQuadraticGradient' filesep 'sweepQuadraticGradient 20230417_012656.mat']);   % numBasis = 200
sct = load(['dataSet' filesep 'sweepCellTemperature' filesep 'sweepCellTemperature_v2 20250626_171113.mat']);   % numBasis = 100, Rrel=1/(100us), w/L=0.4, Imax=50mW/cm^2


% disp. pars. used in calculation
% spp.cellPars_129
% spp.cellPars_131
% spp.result.pumpBeamProfile
% spp.result.I_max_list
% spp.result.numBasis
% spp.aperture

%% 数值结果 与 二阶梯度解析结果 比较
% figure();
% G2 = sqg.result.G2_list';
% g129 = sqg.cellPars_129.gXe;
% g131 = spp.cellPars_131.gXe;
% L = sqg.cellPars_129.L(1);
% D = sqg.cellPars_129.D;
% lambda129 = sqg.cellPars_129.lambda(1);
% lambda131 = sqg.cellPars_131.lambda(1);
% 
% theory_1st_order_129 = 1i*g129*G2*L^2/12*(1-2/15*lambda129);
% theory_1st_order_131 = 1i*g131*G2*L^2/12*(1-2/15*lambda131);
% theory_2nd_order_129 = g129^2*G2.^2*L^6/7560/D*(1-2/15*lambda129);
% theory_2nd_order_131 = g131^2*G2.^2*L^6/7560/D*(1-2/15*lambda131);
% chi1 = 6.68056e-8;
% chi2 = 0.646886;
% theory_3rd_order_129 = -1i*g129^3*G2.^3*L^10/D^2*chi1*(1+chi2*lambda129);
% theory_3rd_order_131 = -1i*g131^3*G2.^3*L^10/D^2*chi1*(1+chi2*lambda131);
% 
% 
% subplot(1,2,1);
% plot(sqg.result.G2_list , abs(-sqg.result.ptb_129.order1st - theory_1st_order_129)./abs(-sqg.result.ptb_129.order1st), 'DisplayName', 's1 error');  hold on;
% plot(sqg.result.G2_list , abs(-sqg.result.ptb_129.order2nd - theory_2nd_order_129)./abs(-sqg.result.ptb_129.order2nd), 'DisplayName', 's2 error');
% % bA_theory = -1i*(theory_1st_order_129/g129 - theory_1st_order_131/g131);   % real bA
% % bA_numeric = +1i*(sqg.result.ptb_129.order1st/g129 - sqg.result.ptb_131.order1st/g131);
% bA_theory = -1i*(theory_3rd_order_129/g129^3 - theory_3rd_order_131/g131^3);    % 3
% bA_numeric = +1i*((sqg.result.ptb_129.order3rd_a+sqg.result.ptb_129.order3rd_b)/g129^3 - (sqg.result.ptb_131.order3rd_a+sqg.result.ptb_131.order3rd_b)/g131^3);
% plot(sqg.result.G2_list , abs(bA_numeric - bA_theory)./abs(bA_numeric), 'DisplayName', 'bA error');
% 
% subplot(1,2,2);
% plot(sqg.result.G2_list , abs(-sqg.result.ptb_129.order3rd_a -sqg.result.ptb_129.order3rd_b - theory_3rd_order_129)./abs(theory_3rd_order_129), 'DisplayName', 's3 error');  hold on;
% % plot(sqg.result.G2_list , -1i*sqg.result.ptb_129.order3rd_a -1i*sqg.result.ptb_129.order3rd_b, 'DisplayName', 's3 numeric');  hold on;
% % plot(sqg.result.G2_list , 1i*theory_3rd_order_129, 'DisplayName', 's3 theory');  hold on;

%% plot
hf = figure('Unit', 'Centimeters', 'Position', [1 1 9 15]);
mkr_size = 2; lw = 0.75;
font_sz = 8;
ax1 =  axes('OuterPosition', [0, 0.67, 1.00, 0.33]); 
ax2 =  axes('OuterPosition', [0, 0.34, 1.00, 0.33]); 
ax3 =  axes('OuterPosition', [0, 0.01, 1.00, 0.33]); 

colorList = colororder();

g129 = spp.cellPars_129.gXe;
g131 = spp.cellPars_131.gXe;
R0 = g129/g131;
re = sqg.result;

plot(ax1, re.G2_list, -sqg.cellPars_129.B0*( abs(re.R) - abs(R0))/(R0), '-', 'Color', colorList(1,:), 'DisplayName', 'Numerical Diagonalization');  hold(ax1, 'on'); grid(ax1, 'on');
plot(ax1, re.G2_list, -sqg.cellPars_129.B0*( abs(re.ptb_R) - abs(R0) )/(R0), '--', 'Color', colorList(2,:), 'DisplayName', 'Perturbation Approx. $b_{\rm A}$');
plot(ax1, re.G2_list, -1i*(R0*re.ptb_131.order1st - re.ptb_129.order1st)/g131/(R0), '-.', 'Color', colorList(5,:), 'DisplayName', '$b_{\rm A}^{(1)}$'); 
plot(ax1, re.G2_list, -1i*(R0*re.ptb_131.order3rd_a - re.ptb_129.order3rd_a)/g131/(R0)-1i*(R0*re.ptb_131.order3rd_b - re.ptb_129.order3rd_b)/g131/(R0), 'k-.', 'DisplayName', '$b_{\rm A}^{(3)}$'); 
% plot(ax2, re.G2_list, (sqg.cellPars_131.lambda(1) - sqg.cellPars_129.lambda(1))*re.G2_list*sqg.cellPars_129.L(1)^2/90, '--g' , 'DisplayName', '$b_{\rm A}^{(1)}$ Theory' );

xlabel(ax1, 'Quadratic Gradient $G_2$ $\rm (nT/cm^2)$', 'interpreter', 'latex');
ylabel(ax1, '$b_{\rm A}$ (nT)', 'interpreter', 'latex');

re = spp.result;

plot(ax2, 1e3*re.I_max_list, -spp.cellPars_129.B0*( abs(re.R) - abs(R0))/(R0), '-', 'Color', colorList(1,:),'DisplayName', 'Numerical Diagonalization');  hold(ax2, 'on'); grid(ax2, 'on');
plot(ax2, 1e3*re.I_max_list, -spp.cellPars_129.B0*( abs(re.ptb_R) - abs(R0) )/(R0), '--', 'Color', colorList(2,:), 'DisplayName', 'Perturbation Approx. $b_{\rm A}$');
plot(ax2, 1e3*re.I_max_list, -1i*(R0*re.ptb_131.order1st - re.ptb_129.order1st)/g131/(R0), '-.', 'Color', colorList(5,:), 'DisplayName', '$b_{\rm A}^{(1)}$'); 
plot(ax2, 1e3*re.I_max_list, -1i*(R0*re.ptb_131.order3rd_a - re.ptb_129.order3rd_a)/g131/(R0)-1i*(R0*re.ptb_131.order3rd_b - re.ptb_129.order3rd_b)/g131/(R0), 'k-.', 'DisplayName', '$b_{\rm A}^{(3)}$'); 

xlabel(ax2, 'Center Intensity $h \nu \Phi_{\rm m}$ $\rm (mW/cm^2)$', 'interpreter', 'latex');
ylabel(ax2, '$b_{\rm A}$ (nT)', 'interpreter', 'latex');
% lg = legend(ax1, {}, 'location', 'best', 'interpreter', 'latex');
% xlim(ax2, [0,1]);

re = sct.result;

plot(ax3, re.cellTemp_list, -sct.cellPars_129.B0*( abs(re.R) - abs(R0))/(R0), '-', 'Color', colorList(1,:),'DisplayName', 'Numerical Diagonalization');  hold(ax3, 'on'); grid(ax3, 'on');
plot(ax3, re.cellTemp_list, -sct.cellPars_129.B0*( abs(re.ptb_R) - abs(R0) )/(R0), '--', 'Color', colorList(2,:), 'DisplayName', 'Perturbation Approx. $b_{\rm A}$');
plot(ax3, re.cellTemp_list, -1i*(R0*re.ptb_131.order1st - re.ptb_129.order1st)/g131/(R0), '-.', 'Color', colorList(5,:), 'DisplayName', '$b_{\rm A}^{(1)}$'); 
plot(ax3, re.cellTemp_list, -1i*(R0*re.ptb_131.order3rd_a - re.ptb_129.order3rd_a)/g131/(R0)-1i*(R0*re.ptb_131.order3rd_b - re.ptb_129.order3rd_b)/g131/(R0), 'k-.', 'DisplayName', '$b_{\rm A}^{(3)}$'); 

xlabel(ax3, 'Cell Temperature $\rm (^\circ C)$', 'interpreter', 'latex');
ylabel(ax3, '$b_{\rm A}$ (nT)', 'interpreter', 'latex');
lg = legend(ax3, {}, 'location', 'best', 'interpreter', 'latex');
% set(lg, 'Position',[0.417230206433547 0.600407764162535 0.514793878791094 0.152703114585446]);

annotation(hf, 'textbox', [0.00 0.94 0.05 0.05], 'String', '(a)', 'FontSize', 8, 'EdgeColor', 'none');
annotation(hf, 'textbox', [0.00 0.61 0.05 0.05], 'String', '(b)', 'FontSize', 8, 'EdgeColor', 'none');
annotation(hf, 'textbox', [0.00 0.28 0.05 0.05], 'String', '(c)', 'FontSize', 8, 'EdgeColor', 'none');

% set font size
allAxes = hf.findobj('Type','axes');
allLegend = hf.findobj('Type','legend');
allLine = hf.findobj('Type','line');
allErrorbar = hf.findobj('Type','errorbar');
allTextbox = hf.findobj('Type','textbox');
set([allAxes;allLegend;allTextbox],'fontsize',font_sz);
set([allLine;allErrorbar],'LineWidth',lw, 'MarkerSize', mkr_size);

%%
% Export figure
print(hf, ['Fig10'], '-dpng', '-r600');
hf.PaperSize = hf.Position(3:4);
print(hf, ['Fig10'], '-dpdf', '-painters');
    
