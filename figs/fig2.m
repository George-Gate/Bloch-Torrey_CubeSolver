clear; clc;

% 这里的数据文件应该跟 Fig4(b) 用的数据文件一致
warning off;
load(['dataSet' filesep 'sweepPumpPower' filesep 'sweepPumpPower 20250626_170731.mat']);   % numBasis = 100, Rrel=1/(100us), w/L=0.4
warning on;

%%
gS = 2;
mu0 = 4*pi*1e-7;       % N/A^2     % 1 T = 1 N/(m*A)
muB = 927.400968e-26;  % J/T
h = 6.62606957e-34;    % Plank's constant, J*s
g129 = cellPars_129.gXe;
g131 = cellPars_131.gXe;
R0 = g129/g131;
re = result;

% characteristic cell length
try
    re.characteristic_info.nRb;
catch
    solver = Table0.CubeSolver2();
    re.characteristic_info.nRb = solver.getRbNumberDensity(cellPars_129.Temp);
end
Bm = -cellPars_129.pumpPolarization/3*cellPars_129.kappa_RbXe*(1e6*re.characteristic_info.nRb)*gS*mu0*muB * 1e9;  % in unit of nT
L0 = pi*(cellPars_131.lambda(1) - cellPars_129.lambda(1))^(1/4)*sqrt( cellPars_129.D/g129/Bm );                   % in unit of cm

hf = figure('Unit', 'Centimeters', 'Position', [1 1 18 5]);
mkr_size = 2; lw = 0.75;
font_sz = 7;
ax1 =  axes('OuterPosition', [0.01, 0.01, 0.33, 0.99]); 
ax2 =  axes('OuterPosition', [0.34, 0.01, 0.33, 0.99]); 
ax3 =  axes('OuterPosition', [0.67, 0.01, 0.33, 0.99]); 

clist = colororder();

% subplot(1,3,1);
plot(ax1, 1e3*re.I_max_list, -imag( re.ptb_129.order1st )/g129, '-', 'DisplayName', '$\delta B^{\rm (1)}_{\rm 129}$'); hold(ax1, 'on'); grid(ax1, 'on');
plot(ax1, 1e3*re.I_max_list, -imag( re.ptb_131.order1st )/g131, '--', 'DisplayName', '$\delta B^{\rm (1)}_{\rm 131}$');   
plot(ax1, 1e3*re.I_max_list, -1000*(  imag( re.ptb_131.order1st )/g131 - imag( re.ptb_129.order1st )/g129  ), '-.', 'DisplayName', '$10^3 \times (\delta B^{\rm (1)}_{\rm 131} - \delta B^{\rm (1)}_{\rm 129})$');
xlabel(ax1, 'Center Intensity $h \nu \Phi_{\rm m}$ $\rm (mW/cm^2)$', 'interpreter', 'latex');
ylabel(ax1, 'First Order Correction (nT)', 'interpreter', 'latex');
lg = legend(ax1, {}, 'interpreter', 'latex', 'location', 'best');
set(lg, 'position', [0.075662641420837,0.753617353341379,0.202699332077909,0.23452116204877]);


% subplot(1,3,2);   
plot(ax2, 1e3*re.I_max_list, abs( re.ptb_129.order3rd_a / g129^3 )*g129^2, '-', 'color', clist(1,:), 'DisplayName', '$\delta B^{\rm (3a)}_{\rm 129}$');  hold(ax2, 'on'); grid(ax2, 'on');
plot(ax2, 1e3*re.I_max_list, abs( re.ptb_131.order3rd_a / g131^3 )*g129^2, '--', 'color', clist(2,:), 'DisplayName', '$\delta B^{\rm (3a)}_{\rm 131}$');

plot(ax2, 1e3*re.I_max_list, abs( re.ptb_129.order3rd_b / g129^3 )*g129^2, '-', 'color', clist(5,:), 'DisplayName', '$\delta B^{\rm (3b)}_{\rm 129}$');  
plot(ax2, 1e3*re.I_max_list, abs( re.ptb_131.order3rd_b / g131^3 )*g129^2, '--', 'color', clist(4,:), 'DisplayName', '$\delta B^{\rm (3b)}_{\rm 131}$');

xlabel(ax2, 'Center Intensity $h \nu \Phi_{\rm m}$ $\rm (mW/cm^2)$', 'interpreter', 'latex');
ylabel(ax2, '$\left| \delta B^{(k)}_{u} / \gamma_{u}^2 \right| \times \gamma_{129}^2 $ (nT)', 'interpreter', 'latex');
legend(ax2, {}, 'interpreter', 'latex', 'location', 'best');


% subplot(1,3,3);
plot(ax3, 1e3*re.I_max_list, -imag( re.ptb_129.order3rd_a + re.ptb_129.order3rd_b )/g129, '-', 'DisplayName', '$\delta B^{\rm (3)}_{\rm 129}$');  hold(ax3, 'on'); grid(ax3, 'on');
plot(ax3, 1e3*re.I_max_list, -imag( re.ptb_131.order3rd_a + re.ptb_131.order3rd_b )/g131, '--', 'DisplayName', '$\delta B^{\rm (3)}_{\rm 131}$');  
plot(ax3, 1e3*re.I_max_list, -imag( re.ptb_129.order3rd_a + re.ptb_129.order3rd_b )/g129 + imag( re.ptb_131.order3rd_a + re.ptb_131.order3rd_b )/g131, 'k-.', 'DisplayName', '$\delta B^{\rm (3)}_{\rm 129} - \delta B^{\rm (3)}_{\rm 131}$'); 

xlabel(ax3, 'Center Intensity $h \nu \Phi_{\rm m}$ $\rm (mW/cm^2)$', 'interpreter', 'latex');
ylabel(ax3, 'Third Order Correction (nT)', 'interpreter', 'latex');
legend(ax3, {}, 'interpreter', 'latex', 'location', 'best');



% set font size
allAxes = hf.findobj('Type','axes');
allLegend = hf.findobj('Type','legend');
allLine = hf.findobj('Type','line');
allErrorbar = hf.findobj('Type','errorbar');
allTextbox = hf.findobj('Type','textbox');
set([allAxes;allLegend;allTextbox],'fontsize',font_sz);
set([allLine;allErrorbar],'LineWidth',lw, 'MarkerSize', mkr_size);

annotation('textbox', [0.00 0.95 0.05 0.05], 'String', '(a)', 'FontSize', 8, 'EdgeColor', 'none');
annotation('textbox', [0.33 0.95 0.05 0.05], 'String', '(b)', 'FontSize', 8, 'EdgeColor', 'none');
annotation('textbox', [0.66 0.95 0.05 0.05], 'String', '(c)', 'FontSize', 8, 'EdgeColor', 'none');


% Export figure
print(hf, ['Fig2'], '-dpng', '-r600');
hf.PaperSize = hf.Position(3:4);
print(hf, ['Fig2'], '-dpdf', '-painters');
