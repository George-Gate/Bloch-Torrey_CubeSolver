clear; clc;

%% load data
load(['dataSet' filesep 'compare3_sweepPumpPower' filesep 'sweepPumpPower_v2 20250626_191458.mat']);

%% calc. 近似解析式
gS = 2;
mu0 = 4*pi*1e-7;       % N/A^2     % 1 T = 1 N/(m*A)
muB = 927.400968e-26;  % J/T
h = 6.62606957e-34;    % Plank's constant, J*s
g129 = cellPars_129.gXe;
g131 = cellPars_131.gXe;
R0 = g129/g131;

re = result;
L = cellPars_129.L(1);
D = cellPars_129.D;
OD = cellPars_129.OD;
Bm = -cellPars_129.pumpPolarization/3*cellPars_129.kappa_RbXe*(1e6*re.characteristic_info.nRb)*gS*mu0*muB * 1e9;  % in unit of nT
I_star = re.characteristic_info.I_star;
base_ampl_weak = Bm*re.I_max_list/I_star;   % Bm*Phi_m/Phi_star,  in unit of nT
base_ampl_strong = Bm*I_star./re.I_max_list;   % Bm*Phi_star/Phi_m,  in unit of nT

cfs_weak = get_characteristic_functions_weakPump(re.pumpBeamProfile.w/L, OD);
cfs_strong_1 = get_characteristic_functions_strongPump(re.pumpBeamProfile.w/L);


bA1st_approx_weak = base_ampl_weak*(cellPars_129.lambda(1) - cellPars_131.lambda(1))*(cfs_weak.F1*cfs_weak.G1+cfs_weak.F2*cfs_weak.G2);
bA1st_approx_strong = base_ampl_strong*(cellPars_129.lambda(1) - cellPars_131.lambda(1))*(cfs_strong_1.F5);    % 强光极限近似

bA1st = -1i*(R0*re.ptb_131.order1st - re.ptb_129.order1st)/g131/(R0);

fprintf('L=%g cm, w/L=%g, I*=%.3fmW, OD=%g, Bm=%.3fnT\n', L, result.pumpBeamProfile.w/L, 1e3*I_star, OD, Bm);

%% plot
hf = figure('Unit', 'Centimeters', 'Position', [1 1 9 7]);
mkr_size = 2; lw = 0.75;
font_sz = 8;
ax =  axes('OuterPosition', [0, 0, 1.00, 1.00]); 

colorList = colororder();

x_data = (re.I_max_list/I_star);

loglog(ax, x_data, 1e3*bA1st, 'k-', 'DisplayName', '$b_{\rm A}^{(1)}$'); hold(ax, 'on'); grid(ax, 'on');
loglog(ax, x_data, 1e3*bA1st_approx_weak, '--', 'Color', colorList(1,:), 'DisplayName', '$b_{\rm A}^{(1)}$ weak approx.');
loglog(ax, x_data, 1e3*bA1st_approx_strong, '--', 'Color', colorList(2,:), 'DisplayName', '$b_{\rm A}^{(1)}$ strong approx.');

xlabel(ax, 'Normalized Center Intensity $\Phi_{\rm m}/\Phi^*$', 'interpreter', 'latex');
ylabel(ax, '$\left| b_{\rm A}^{(1)} \right|$ (pT)', 'interpreter', 'latex');
xlim(ax, [1e-3, 1e3]);
ylim(ax, [1e-3, 1e+2]);
xticks(ax, [1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3]);
lg = legend(ax, {}, 'location', 'best', 'interpreter', 'latex');


% 
% semilogx(ax, lambda, kappaL(:,end:-1:1)/pi);  hold on;
% semilogx(ax, lambda, sqrt(2*lambda)/pi, '--', 'color', colorList(5,:));
% % semilogx(ax, lambda, sqrt(2*lambda+lambda.*lambda)/pi, '-.', 'color', colorList(5,:));
% semilogx(ax, lambda, (pi + 2*lambda/pi )/pi, '--', 'color', colorList(4,:));
% semilogx(ax, lambda, (2*pi + 2*lambda/(2*pi) )/pi, '--', 'color', colorList(3,:));
% semilogx(ax, lambda, (3*pi + 2*lambda/(3*pi) )/pi, '--', 'color', colorList(2,:));
% semilogx(ax, lambda, (4*pi + 2*lambda/(4*pi) )/pi, '--', 'color', colorList(1,:));
% grid on;
% xlabel(ax, '$\lambda$', 'interpreter', 'latex');
% ylabel(ax, '$\kappa_p L / \pi$', 'interpreter', 'latex');
% legend(ax, num2str([num_eigs-1:-1:0]', '$p=%d$'), 'Position',[0.217034420518732 0.63397709253039 0.214137804700447 0.253529405225726], 'interpreter', 'latex');
% xlim(ax, minmax(lambda(:)'));
% ylim(ax, [0 5]);
% xticks(ax, [1e-3,1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3])

% set font size
allAxes = hf.findobj('Type','axes');
allLegend = hf.findobj('Type','legend');
allLine = hf.findobj('Type','line');
allErrorbar = hf.findobj('Type','errorbar');
allTextbox = hf.findobj('Type','textbox');
set([allAxes;allLegend;allTextbox],'fontsize',font_sz);
set([allLine;allErrorbar],'LineWidth',lw, 'MarkerSize', mkr_size);


% Export figure
print(hf, ['Fig7'], '-dpng', '-r600');
hf.PaperSize = hf.Position(3:4);
print(hf, ['Fig7'], '-dpdf', '-painters');

   
