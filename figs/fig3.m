clear; clc;

warning off


% (xc, yc)/L = (0.02,0.03)
dataPath = ['dataSet' filesep 'sweepPumpPower_varyBeamWidth2'];
vBW{1} = load([dataPath filesep 'sweepPumpPower_varyBeamWidth_v2 L=0.30cm 20250626_173539.mat']);   % numBasis = 100, L=0.3cm;
vBW{2} = load([dataPath filesep 'sweepPumpPower_varyBeamWidth_v2 L=0.50cm 20250626_174823.mat']);   % numBasis = 100, L=0.5cm
vBW{3} = load([dataPath filesep 'sweepPumpPower_varyBeamWidth_v2 L=0.80cm 20250626_182816.mat']);   % numBasis = 100, L=0.8cm


warning on;

%%
xAxis = 'power';   % power / intensity
if strcmp(xAxis, 'intensity')
    is_x_pow = false;
else
    is_x_pow = true;
end

gS = 2;
mu0 = 4*pi*1e-7;       % N/A^2     % 1 T = 1 N/(m*A)
muB = 927.400968e-26;  % J/T
h = 6.62606957e-34;    % Plank's constant, J*s

g129 = vBW{1}.cellPars_129.gXe;
g131 = vBW{1}.cellPars_131.gXe;
R0 = g129/g131;

hf = figure('Unit', 'Centimeters', 'Position', [1 1 18 6]);
mkr_size = 2; lw = 0.75;
font_sz = 7;

ax1 =  axes('OuterPosition', [0.01, 0.01, 0.35, 0.97]); 
ax2 =  axes('OuterPosition', [0.34, 0.01, 0.35, 0.97]); 
ax3 =  axes('OuterPosition', [0.67, 0.01, 0.35, 0.97]);
ax = [ax1, ax2, ax3];

clist = colororder();


numBasis_list = cellfun(@(x)x.result{1, 1}.numBasis, vBW);
fprintf('numBasis:[%s]\n', num2str(numBasis_list(:)','%d '));
re = vBW{1}.result{1};
diff_lambda = vBW{1}.cellPars_131.lambda(1) - vBW{1}.cellPars_129.lambda(1);
D = vBW{1}.cellPars_131.D;
Bm = -vBW{1}.cellPars_129.pumpPolarization/3*vBW{1}.cellPars_129.kappa_RbXe*(1e6*re.characteristic_info.nRb)*gS*mu0*muB * 1e9;  % in unit of nT
L0 = pi*abs(diff_lambda)^0.25*sqrt( D/abs(g129*Bm) );

wanted_w_L = [0.25, 0.3, 0.5, 0.8];
xlims = {[0, 0.015], [0, 0.030], [0, 0.080]};

for iL = 1:length(vBW)
    beamProfiles = cellfun(@(re)re.pumpBeamProfile, vBW{iL}.result);
    L = vBW{iL}.cellPars_129.L(1);
    fprintf('(%d)OD: %g, (xc,yc)/L: [\n', iL, vBW{iL}.cellPars_129.OD);
    disp(num2str([vertcat(beamProfiles.xc)/L, vertcat(beamProfiles.yc)/L],'%g %g'));
    fprintf(']\n');
    
    w_L_list = vBW{iL}.w_L_list;
    w_L_sel = arrayfun(@(x)find(abs(w_L_list - x)<3*eps), wanted_w_L);
    % 确保所有子图中的 w_L 一致
    tmp = w_L_list(w_L_sel);
    assert(all( abs(tmp(:)- wanted_w_L(:)) <2*eps ));
    
    iC = 1;
    
    for ii = w_L_sel
        re = vBW{iL}.result{ii};

        % calc bA curves
        bA = -vBW{iL}.cellPars_129.B0*( abs(re.R) - abs(R0))/(R0);
        bA_1st = -1i*(R0*re.ptb_131.order1st - re.ptb_129.order1st)/g131/(R0);
        

        tmp_Ilist = re.I_max_list;
        tmp_P_beam = re.P_beam_list(end)/re.I_max_list(end)  *tmp_Ilist;   % P_beam 正比于 I_max
        if is_x_pow
            tmp_x = 1e3*tmp_P_beam;    % total power of the whole laser beam
        else
            tmp_x = tmp_Ilist;
        end
        plot(ax(iL), tmp_x, 1e3*bA          , 'Color', clist(iC,:), 'DisplayName', sprintf('%.2f', w_L_list(ii)));  
        hold(ax(iL), 'on'); grid(ax(iL), 'on');
        plot(ax(iL), tmp_x, 1e3*bA_1st, ':', 'Color', 0.9*clist(iC,:), 'DisplayName', '$b_A^{(1)}$', 'Tag','dashed'); 
        iC = iC+1;
    end
    
    if is_x_pow
        xlabel(ax(iL), 'Pump Beam Power $P_{\rm in}$ $\rm (mW)$', 'interpreter', 'latex');
        xlim(ax(iL),1e3*xlims{iL});
    else
        xlabel(ax(iL), 'Center Intensity $I_0$ $\rm (W/cm^2)$', 'interpreter', 'latex');
        xlim(ax(iL),xlims{iL});
    end
    ylabel(ax(iL), '$b_{\rm A}$ (pT)', 'interpreter', 'latex');
    title(ax(iL), sprintf('$L=%.1f~{\\rm cm}$', L), 'interpreter', 'latex');    
    
end

lg_ax = ax1;
% 重排线的顺序，只显示实线的legend
lineList = lg_ax.Children;
N = length(lineList);
swapIdx = [1:2:N-1, 2:2:N];
lg_ax.Children = lineList(swapIdx);
lg = legend(ax1, {}, 'location', 'best', 'interpreter', 'latex');
lg.String = lg.String(1:N/2);
% set(lg, 'Position', [0.185528022648871,0.575894329394458,0.102941175637876,0.379735673025316]);
title(lg, '$w/L$', 'interpreter', 'latex');

% set font size
allAxes = hf.findobj('Type','axes');
allLegend = hf.findobj('Type','legend');
allLine = hf.findobj('Type','line');
allErrorbar = hf.findobj('Type','errorbar');
allTextbox = hf.findobj('Type','textbox');
set([allAxes;allLegend;allTextbox],'fontsize',font_sz);
set([allLine;allErrorbar],'LineWidth',lw, 'MarkerSize', mkr_size);

% 特殊处理标记的线
allLine2 = hf.findobj('Type','line', 'Tag', 'dashed');
set(allLine2,'LineWidth',lw+0.5);

annotation('textbox', [0.00 0.94 0.05 0.05], 'String', '(a)', 'FontSize', 8, 'EdgeColor', 'none');
annotation('textbox', [0.33 0.94 0.05 0.05], 'String', '(b)', 'FontSize', 8, 'EdgeColor', 'none');
annotation('textbox', [0.66 0.94 0.05 0.05], 'String', '(c)', 'FontSize', 8, 'EdgeColor', 'none');

% Export figure
print(hf, ['Fig3_' xAxis], '-dpng', '-r600');
hf.PaperSize = hf.Position(3:4);
print(hf, ['Fig3_' xAxis], '-dpdf', '-painters');
