% 计算不同 w/L 下的 bA vs. Imax, dbA/dxc, dbA/dsigma_abs vs. Imax 曲线
% 只选取三个温度点计算导数
% D, lambda, nRb都随温度变化
% lambda 正比于 L

cellLen = 0.8;
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


fileName = sprintf('sweepPumpPower_varyBeamWidth_v4 L=%.2fcm %s', cellPars.L(1), datestr(now(),'yyyymmdd_HHMMSS'));


% w_L_list = [0.25, 0.3, 0.35, 0.4, 0.5, 0.7, 0.8];
w_L_list = linspace(0.2, 1.0, 100);
% w_L_list = linspace(0.2, 1.0, 10);

rel_dxc = 1e-4;             % 计算dbA/dxc时使用的步长为  rel_dxc*L 或 rel_dxc*4w，取较小值
rel_dsigma = 1e-3;          % 计算dbA/dsigma_abs时使用的步长为 rel_dsigma*sigma_abs

% aperture_list = 3.0;

result_xc = {};
result_sigma_abs = {};
dxc_list = [];
tS = tic;
% calc dbA/dxc
for ii = 1:length(w_L_list)
    pumpBeamProfile.w = w_L_list(ii)*cellPars.L(1);
    I_max_list = P_max/(2*pi*pumpBeamProfile.w^2)*unique([linspace(0, 0.25, 64), linspace(0.005, 1, 96)]);
    % I_max_list = P_max/(2*pi*pumpBeamProfile.w^2)*unique([linspace(0, 0.25, 32), linspace(0.005, 1, 32)]);
    dxc = rel_dxc*min(cellPars_129.L(1), 4*pumpBeamProfile.w);
    dxc_list(1,ii) = dxc;
    for ixc = 1:2
        % 更新xc
        pBP = pumpBeamProfile;
        pBP.xc = pBP.xc + dxc*(ixc-1.5);  % +-0.5*dxc
        
        result_xc{ii, ixc} = solver.get_pump_power_dependence(I_max_list, pBP, cellPars_129, cellPars_131, false);
    end
    fprintf('(Step 1) %d/%d finished. Time elapsed: %.3f hr\n', ii,length(w_L_list), toc(tS)/3600 );
end
% calc dbA/dsigma_abs
for ii = 1:length(w_L_list)
    pumpBeamProfile.w = w_L_list(ii)*cellPars.L(1);
    I_max_list = P_max/(2*pi*pumpBeamProfile.w^2)*unique([linspace(0, 0.25, 64), linspace(0.005, 1, 96)]);
    % I_max_list = P_max/(2*pi*pumpBeamProfile.w^2)*unique([linspace(0, 0.25, 32), linspace(0.005, 1, 32)]);
    
    for iS = 1:2
        % 更新OD参数，由于L和温度固定不变，更改OD相当于更改了sigma_abs
        % 程序中计算Phi^*时，也是根据OD先反推得到sigma_abs，因此这里只需要更新OD的值
        cP_129 = cellPars_129;            cP_131 = cellPars_131;
        cP_129.OD = cP_129.OD*(1+(iS-1.5)*rel_dsigma);    % +-0.5*dsigma      
        cP_131.OD = cP_131.OD*(1+(iS-1.5)*rel_dsigma);
        
        result_sigma_abs{ii, iS} = solver.get_pump_power_dependence(I_max_list, pumpBeamProfile, cP_129, cP_131, false);
    end
    fprintf('(Step 2) %d/%d finished. Time elapsed: %.3f hr\n', ii,length(w_L_list), toc(tS)/3600 );
end

save([fileName '.mat'], 'result_xc', 'result_sigma_abs', 'cellPars_129', 'cellPars_131', 'w_L_list', 'rel_dxc', 'rel_dsigma', 'dxc_list', 'fileName');


%%
hf = figure('Position', [100 377 1048 522]);
ax1 = subplot(2,2,1);
ax2 = subplot(2,2,3);
ax3 = subplot(2,2,2);
ax4 = subplot(2,2,4);

g129 = cellPars_129.gXe;
g131 = cellPars_131.gXe;
R0 = g129/g131;

for ii = 1:length(w_L_list)
    re = result_xc{ii, 2};
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
bA_mat_xc = {};
bA_mat_sigma = {};

g129 = cellPars_129.gXe;
g131 = cellPars_131.gXe;
R0 = g129/g131;

for ii = 1:2
    tmp = cellfun(@(re)re.R, result_xc(:,ii), 'UniformOutput', false);
    R_mat{ii} = cell2mat(tmp');
    bA_mat_xc{ii} = -cellPars_129.B0*( abs(R_mat{ii}) - abs(R0))/(R0);

    tmp = cellfun(@(re)re.R, result_sigma_abs(:,ii), 'UniformOutput', false);
    R_mat{ii} = cell2mat(tmp');
    bA_mat_sigma{ii} = -cellPars_129.B0*( abs(R_mat{ii}) - abs(R0))/(R0);

    tmp = cellfun(@(re)re.P_beam_list, result_xc(:,ii), 'UniformOutput', false);
    P_beam_mat{ii} = cell2mat(tmp');
    w_L_mat{ii} = repmat(w_L_list,  size(bA_mat_xc{ii},1), 1);
end

% xc导数
dbA_dxc_mat = ( bA_mat_xc{2} - bA_mat_xc{1} ) ./dxc_list;   % nT/cm
% sigma_abs导数
dbA_dsigma_mat = ( bA_mat_sigma{2} - bA_mat_sigma{1} ) /rel_dsigma;   % sigma_abs * dbA/dsigma_abs, in unit of nT

u_sigma = 1e-3; % sigma_abs的相对涨落
delta_xc = 1e-4;  % xc涨落, cm

xData = 1e3*P_beam_mat{1};
yData = w_L_mat{1};
zData = {1e6*delta_xc*dbA_dxc_mat; ...
         1e3*(bA_mat_xc{1}+bA_mat_xc{2})/2;  ...
         1e6*u_sigma*dbA_dsigma_mat};
% zData = {zData_new{1}-zData_old{1};  zData_new{2}-zData_old{2};  zData_new{3}-zData_old{3}};
zlabel_str = {['$\delta x \times \partial b_{\rm A}/\partial x_c$, ', sprintf('$\\delta x = %g~{\\rm \\mu m}$',1e4*delta_xc)];
              '$b_{\rm A}$';
              ['$\delta \sigma \times \partial b_{\rm A}/\partial \sigma_{\rm abs} $, ', sprintf('$\\delta \\sigma = 10^{%g} \\sigma_{\\rm abs}$',log10(u_sigma))]};
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