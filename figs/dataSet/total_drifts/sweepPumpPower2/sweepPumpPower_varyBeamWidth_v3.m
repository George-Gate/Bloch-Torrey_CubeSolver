% 计算不同 w/L 下的 bA vs. Imax, dbA/dT vs. Imax 曲线
% 只选取三个温度点计算导数
% D, lambda, nRb都随温度变化
% lambda 正比于 L

% cellLen = 0.3;
solver = Table0.CubeSolver2();

g129 = -2*pi*11.86015e-3;  % rad/s/nT
g131 = +2*pi*3.515769e-3;  % rad/s/nT

% define cell parameters
cellPars = solver.cellPars_example;    % an example of cellPars structure is defined in the constructor of @CubeSolver2
cellPars.L = cellLen*[1, 1, 1];
cellPars.D = 0.2;
cellPars.D_alpha = 1;    % 扩散系数的温度依赖，D \propto T^alpha
cellPars.OD = 6/0.8*cellPars.L(1); 
cellPars.Temp = 110;
cellPars.Rrel = 1/(100e-6);

if abs( cellPars.L(1) - 0.3 ) < 1e-6
    P_max = 15e-3;
elseif abs( cellPars.L(1) - 0.5 ) < 1e-6
    P_max = 50e-3;
elseif abs( cellPars.L(1) - 0.8 ) < 1e-6
    P_max = 100e-3;
else
    error();
end

cellPars_129 = cellPars;
cellPars_129.lambda = 1e-3*[1,1,1]/0.8*cellPars.L(1);   % lambda 正比于 L
cellPars_129.T_lambda = 0.10 * 1.602e-19/1.38e-23;      % lambda的特征温度（对应Ebar）
cellPars_129.gXe = g129;
cellPars_129.G2c = 1/20;      % 1/s

cellPars_131 = cellPars;
cellPars_131.lambda = 2e-2*[1,1,1]/0.8*cellPars.L(1);   % lambda 正比于 L
cellPars_131.T_lambda = 0.15 * 1.602e-19/1.38e-23;      % lambda的特征温度（对应Ebar）
cellPars_131.gXe = g131;
cellPars_131.G2c = 1/20;      % 1/s

% 计算特征气室尺寸
[~, auxInfo] = solver.calc_polarization_field(@(x,y)ones(size(x)), cellPars, @(x,y)ones(size(x)), false);   % 获取Bm
L0 = (   (abs(cellPars_131.lambda(1)-cellPars_129.lambda(1))*pi^4*cellPars.D^2)/( cellPars.L(1)*abs(g129*auxInfo.Bm)^2 )    )^(1/3)

pumpBeamProfile = struct();
% pumpBeamProfile.w = 0.3*cellPars.L(1);   % cm
pumpBeamProfile.xc = 0.02*cellPars.L(1);
pumpBeamProfile.yc = 0.03*cellPars.L(1);
% pumpBeamProfile.xc = 0;
% pumpBeamProfile.yc = 0;
pumpBeamProfile.aperture = 10.0;
circleMask = @(x,y,r) 1./((1+((x.*x+y.*y)/r^2).^20));  % 圆形通光孔
% pumpBeamProfile.aperture_mask = @(x,y)circleMask(x,y,pumpBeamProfile.aperture/2);           % 通光孔外形，通光区域为1，不通光区域为0
pumpBeamProfile.aperture_mask = @(x,y)ones(size(x));

% 定义温度列表，并计算不同温度下的各种参数
Twork = cellPars.Temp;   % degC
dT = 0.05;
cellTemp_list = Twork + [-dT, 0, +dT];
Tlist = cellTemp_list + 273.15;   % K
T0 = Twork + 273.15;   % K
% get sigma_abs from OD and cellTemp    
[OD_list, sigma_abs] = solver.get_OD_list(cellTemp_list, cellPars);
D_list = cellPars.D ./ (T0^cellPars.D_alpha).*(Tlist.^cellPars.D_alpha);
lambda129_list = cellPars_129.lambda(1)* exp(  cellPars_129.T_lambda./Tlist -  cellPars_129.T_lambda/T0   );
lambda131_list = cellPars_131.lambda(1)* exp(  cellPars_131.T_lambda./Tlist -  cellPars_131.T_lambda/T0   );
% 不同温度下的nRb会在 Table0.CubeSolver2 中自动计算

clear T0 Tlist

fileName = sprintf('sweepPumpPower_varyBeamWidth_v3 L=%.2fcm %s', cellPars.L(1), datestr(now(),'yyyymmdd_HHMMSS'));



% w_L_list = [0.25, 0.3, 0.35, 0.4, 0.5, 0.7, 0.8];
w_L_list = linspace(0.2, 1.0, 100);
% w_L_list = linspace(0.2, 1.0, 10);

% aperture_list = 3.0;
result = {};
tS = tic;
for ii = 1:length(w_L_list)
    pumpBeamProfile.w = w_L_list(ii)*cellPars.L(1);
    I_max_list = P_max/(2*pi*pumpBeamProfile.w^2)*unique([linspace(0, 0.25, 64), linspace(0.005, 1, 96)]);
    % I_max_list = P_max/(2*pi*pumpBeamProfile.w^2)*unique([linspace(0, 0.25, 32), linspace(0.005, 1, 32)]);
    for iT = 1:length(cellTemp_list)
        % 更新温度和OD参数，假设sigma_abs是固定不变的，根据nRb来估算OD的变化
        cP_129 = cellPars_129;            cP_131 = cellPars_131;
        cP_129.Temp = cellTemp_list(iT);  cP_131.Temp = cellTemp_list(iT);  
        cP_129.OD = OD_list(iT);          cP_131.OD = OD_list(iT);
        cP_129.D = D_list(iT);            cP_131.D = D_list(iT);
        cP_129.lambda = lambda129_list(iT)*[1,1,1];
        cP_131.lambda = lambda131_list(iT)*[1,1,1];
        
        result{ii, iT} = solver.get_pump_power_dependence(I_max_list, pumpBeamProfile, cP_129, cP_131, false);
    end
    fprintf('%d/%d finished. Time elapsed: %.3f hr\n', ii,length(w_L_list), toc(tS)/3600 );
end

save([fileName '.mat'], 'result', 'cellPars_129', 'cellPars_131', 'w_L_list', 'Twork', 'dT', 'cellTemp_list', 'OD_list', 'fileName');


%%
hf = figure('Position', [667 377 1048 522]);
ax1 = subplot(2,2,1);
ax2 = subplot(2,2,3);
ax3 = subplot(2,2,2);
ax4 = subplot(2,2,4);

g129 = cellPars_129.gXe;
g131 = cellPars_131.gXe;
R0 = g129/g131;

for ii = 1:length(w_L_list)
    re = result{ii, 2};
    plot(ax1, re.I_max_list, re.T2_129, 'DisplayName', sprintf('%.2f', w_L_list(ii)));  hold(ax1, 'on'); grid(ax1, 'on');
%     plot(ax1, re.I_max_list, -1./real(re.ptb_129.eig0), '--'); 

    % dbA / dI0
    bA = -cellPars_129.B0*( abs(re.R) - abs(R0))/(R0);
    bA_fun = griddedInterpolant(re.I_max_list, bA, 'spline');
    bA_chebfun = chebfun(@(x)bA_fun(x), minmax(re.I_max_list(:)'));
    dbA_dI0_chebfun = diff(bA_chebfun,1);
    
    tmp_Ilist = linspace(min(re.I_max_list), max(re.I_max_list), 5000);
    plot(ax2, tmp_Ilist, dbA_dI0_chebfun(tmp_Ilist), 'DisplayName', sprintf('%.2f', w_L_list(ii))); hold(ax2, 'on'); grid(ax2, 'on');

    plot(ax3, re.I_max_list, re.f_129, 'DisplayName', sprintf('%.2f', w_L_list(ii))); hold(ax3, 'on'); grid(ax3, 'on');
%     plot(ax3, re.I_max_list, abs(imag(re.ptb_129.eig0))/2/pi, '--');  
    plot(ax4, re.I_max_list, bA, 'DisplayName', sprintf('%.2f', w_L_list(ii)));  hold(ax4, 'on'); grid(ax4, 'on');
%     plot(ax4, re.I_max_list, cellPars_129.B0*( re.ptb_R - re.ptb_R(1) ), '--'); 
end

xlabel(ax1, 'Center intensity $I_0$ $\rm (W/cm^2)$', 'interpreter', 'latex');
ylabel(ax1, '$T_{2,129}$ (s)', 'interpreter', 'latex');
title(ax1, re.type, 'interpreter', 'none');

xlabel(ax2, 'Center intensity $I_0$ $\rm (W/cm^2)$', 'interpreter', 'latex');
ylabel(ax2, '${\rm d}b_{\rm A}/{\rm d}I_0$ $\rm (nT \cdot cm^2 / W)$', 'interpreter', 'latex');
title(ax2, sprintf( 'Num. Basis: %d, $\\lambda_{129}=%g$, $\\lambda_{131}=%g$' ,  re.numBasis, cellPars_129.lambda(1), cellPars_131.lambda(1)), 'interpreter', 'latex');   % basis info
xlim(ax2,[0,0.5])

xlabel(ax3, 'Center intensity $I_0$ $\rm (W/cm^2)$', 'interpreter', 'latex');
ylabel(ax3, '$f_{129}$ (Hz)', 'interpreter', 'latex');
title(ax3, sprintf( 'Pump Beam: Aperture=%g~cm, center=$(%g,%g)$~cm' ,  re.pumpBeamProfile.aperture, re.pumpBeamProfile.xc, re.pumpBeamProfile.yc ), 'interpreter', 'latex');   % pump beam info

xlabel(ax4, 'Center intensity $I_0$ $\rm (W/cm^2)$', 'interpreter', 'latex');
ylabel(ax4, '$b_{\rm A}$ (nT)', 'interpreter', 'latex');
title(ax4, sprintf( '$L=%.2f~{\\rm cm}$, OD($%g~{\\rm ^\\circ C}$) = %g, $D=%g~{\\rm cm^2/s}$' ,  cellPars_129.L(1), cellPars_129.Temp, cellPars_129.OD, cellPars_129.D ), 'interpreter', 'latex');   % pump beam info
xlim(ax4,[0,0.5])

lg = legend(ax3, {}, 'location', 'best');
title(lg, '$w/L$', 'interpreter','latex');

print(hf, [fileName '.png'], '-dpng', '-r330');

%% 3D plot
% construct data
P_beam_mat = {};
w_L_mat = {};
R_mat = {};
bA_mat = {};

g129 = cellPars_129.gXe;
g131 = cellPars_131.gXe;
R0 = g129/g131;

for iT = 1:length(cellTemp_list)
    tmp = cellfun(@(re)re.P_beam_list, result(:,iT), 'UniformOutput', false);
    P_beam_mat{iT} = cell2mat(tmp');
    tmp = cellfun(@(re)re.R, result(:,iT), 'UniformOutput', false);
    R_mat{iT} = cell2mat(tmp');
    bA_mat{iT} = -cellPars_129.B0*( abs(R_mat{iT}) - abs(R0))/(R0);
    w_L_mat{iT} = repmat(w_L_list,  size(bA_mat{iT},1), 1);
end

% 温度导数
dbA_dT_mat = ( bA_mat{3} - bA_mat{1} ) / 2/dT;   % nT/K
% 功率导数
iT=2;
dbA_dPeam_mat = nan(size(bA_mat{iT}));
for ii = 1:size(bA_mat{iT}, 2)
    bA_fun = griddedInterpolant(P_beam_mat{iT}(:,ii), bA_mat{iT}(:,ii), 'spline');
    bA_chebfun = chebfun(@(x)bA_fun(x), minmax(P_beam_mat{iT}(:,ii)'));
    dbA_dPeam_chebfun = diff(bA_chebfun,1);
    dbA_dPeam_mat(:,ii) = dbA_dPeam_chebfun(P_beam_mat{iT}(:,ii));    % nT/W
end

uP = 1e-3; % 功率涨落
deltaT = 10e-3;  % 温度涨落

xData = 1e3*P_beam_mat{2};
yData = w_L_mat{2};
zData = {1e6*deltaT*dbA_dT_mat; ...
         1e3*bA_mat{2};  ...
         1e6*uP*P_beam_mat{2}.*dbA_dPeam_mat};
zlabel_str = {['$\delta T \times {\rm d}b_{\rm A}/{\rm d}T$, ', sprintf('$\\delta T = %g~{\\rm mK}$',1e3*deltaT)];
              '$b_{\rm A}$';
              ['$\delta P \times {\rm d}b_{\rm A}/{\rm d}P_{\rm in} $, ', sprintf('$\\delta P = 10^{%g} P_{\\rm in}$',log10(uP))]};
unit_str = {'$(\rm fT)$';
            '$(\rm pT)$';
            '$(\rm fT)$'};
          
hf2 = figure('Position', [102 426 1603 402]);
for ii = 1:3
    ax = subplot(1,3,ii);
    sObj = pcolor(ax, xData, yData, zData{ii});   hold(ax, 'on');  box(ax, 'on');  grid(ax, 'on');
    set(sObj,  'FaceColor', 'interp', 'EdgeColor', 'none');
    contour(ax, xData, yData, zData{ii}, 'LineColor', 0.2*[1 1 1]);
    
    caxis(ax, max(abs(caxis(ax)))*[-1,+1]);   % colorbar 上下对称
    crameri('vik', 'pivot', 0);
    set(ax,'Layer','top');
    xlabel(ax, 'Laser Beam Power $P_{\rm in}$ $\rm (mW)$', 'interpreter', 'latex');
    ylabel(ax, '$w/ L$', 'interpreter', 'latex');
    title(ax, zlabel_str{ii}, 'interpreter', 'latex');
    cObj = colorbar(ax);
    title(cObj, unit_str{ii}, 'interpreter', 'latex');
end

mkr_size = 2; lw = 0.75;
font_sz = 14;
% set font size
allAxes = hf2.findobj('Type','axes');
allLegend = hf2.findobj('Type','legend');
allLine = hf2.findobj('Type','line');
allErrorbar = hf2.findobj('Type','errorbar');
allTextbox = hf2.findobj('Type','textbox');
set([allAxes;allLegend;allTextbox],'fontsize',font_sz);
set([allLine;allErrorbar],'LineWidth',lw, 'MarkerSize', mkr_size);

print(hf2, [fileName '_2.png'], '-dpng', '-r330');