clear; clc;

load(['dataSet' filesep 'total_drifts' filesep 'preparedData.mat']);

uP = 1e-3; % 功率相对涨落
deltaT = 10e-3;  % 温度涨落
delta_xc = 1e-4;  % xc涨落, cm
u_sigma = 1e-3; % sigma_abs的相对涨落
u_Gamma = 2e-3; % Rrel的相对涨落


%% plot
hf = figure('Unit', 'Centimeters', 'Position', [1 5 18 10]);
mkr_size = 2; lw = 0.75;
font_sz = 7;

x1 = 0.01;  x2 = 0.34;  x3 = 0.67;
y1 = 0.01;  y2 = 0.49;
HH = 0.49;  WW = 0.36;
ax1 =  axes('OuterPosition', [x1, y2, WW+0.008, HH]); 
ax2 =  axes('OuterPosition', [x2, y2, WW      , HH]); 
ax3 =  axes('OuterPosition', [x3, y2, WW      , HH]);

ax4 =  axes('OuterPosition', [x1, y1, WW+0.008, HH+0.01]); 
ax5 =  axes('OuterPosition', [x2, y1, WW      , HH+0.01]); 
ax6 =  axes('OuterPosition', [x3, y1, WW      , HH+0.01]);

axList = {ax3, ax2, ax1};
axList2 = {ax6, ax5, ax4};

colorList = colororder();
colorList(8,:) = [51 204 51]/255;


% characteristic cell length
d = dataList{1}{1};
Bm_at_Twork = d.result{1, 2}.characteristic_info.Bm;
p129 = d.cellPars_129;
p131 = d.cellPars_131;
L0 = (  abs(p131.lambda(1) - p129.lambda(1))*pi^4*p129.D^2 / (p129.L(1)*abs(p129.gXe*Bm_at_Twork)^2)  )^(1/3);

for ii = 1:3
    dA = dataList{ii}{1};
    dB = dataList{ii}{2};
    fprintf('L=%gcm\n', dA.cellPars_129.L(1));
    
    L = dA.cellPars_129.L(1);
    


    % 计算各因素导致的fluctuation
    % in unit of nT
    flu_P = uP*dA.P_beam_mat{2}.*dA.dbA_dPeam_mat;
    flu_T = deltaT*dA.dbA_dT_mat;
    flu_xc = delta_xc*dB.dbA_dxc_mat;
    flu_sigma = u_sigma*dB.rel_dbA_dsigma_mat;
    flu_Gamma = u_Gamma*dB.rel_dbA_dRrel_mat;
    % calculate polarization
    w_mat = dA.w_L_mat{2}*dA.cellPars_129.L(1);
    Phim_Phistar = dA.P_beam_mat{2}./(2*pi*w_mat.*w_mat)/dA.I_star;   % dA.I_star in unit of W/cm^2
    polarization = Phim_Phistar./(1+Phim_Phistar);
    % average polarization
    Bm = dA.result{1,2}.characteristic_info.Bm;
    tmp_mat = cellfun(@(re)re.ptb_129.order1st, dA.result(:,2), 'UniformOutput', false);
    tmp_mat = cell2mat(tmp_mat');
    avg_pol = abs(imag(tmp_mat)/dA.cellPars_129.gXe/Bm);
    
    
    % evaluation function
    delta_bA = sqrt(flu_P.*flu_P + flu_T.*flu_T + flu_xc.*flu_xc + flu_sigma.*flu_sigma + flu_Gamma.*flu_Gamma)./(avg_pol);

    if abs(L-0.3)<1e-6
        minP = 0.5e-6;   % 筛选画图数据，最小Pin
        contour_levels = [3, Inf];   % 等高线图的levels
    elseif abs(L-0.5)<1e-6
        minP = 1.5e-6;
        contour_levels = [6, Inf];
    elseif abs(L-0.8)<1e-6
        minP = 3e-6;
        contour_levels = [12, Inf];
    end
    
    % idx = find(dA.result{1,2}.P_beam_list >= minP);
    idx = 1;
    
    xData = 1e3*dA.P_beam_mat{2}(idx:end,:);
    yData = dA.w_L_mat{2}(idx:end,:);
    zData = {1e3*dA.bA_mat{2};           % bA      , pT
             1e6*delta_bA(idx:end,:)};  % delta_bA , fT
    avg_pol_sel = avg_pol(idx:end,:);

    % 重新构建xData, yData矩阵，以消除数值误差。但需要校验数值误差不大于1e-10
    [xData2, yData2] = ndgrid(xData(:,1), yData(1,:));
    assert(max(max(abs(xData2-xData))) < 1e-10);
    assert(max(max(abs(yData2-yData))) < 1e-10);
    [xMat, yMat] = meshgrid( linspace(min(xData(:)),max(xData(:)),500), linspace(min(yData(:)),max(yData(:)),500)   );
    zMat = {};
    for kk = 1:2
        zData{kk}(isinf(zData{kk})) = 0;
        % interpret data to make plot more smooth
        F = griddedInterpolant(xData2, yData2, zData{kk}, 'spline');
        zMat{kk} = F(xMat, yMat);
    end

    % find the minimum of delta_bA
    tmp_delta_bA_mat = zMat{2};
    filterIdx =  xMat/max(xMat(:)) < 0.01;   % delta_bA=0 at P_in=0, filter those points
    tmp_delta_bA_mat(filterIdx) = nan;
    [wP_target(ii), idx] = min(tmp_delta_bA_mat,[],'all')
    wP_pow(ii) = xMat(idx)
    wP_w_L(ii) = yMat(idx)
    % wP_polarization(ii) = polarization(idx)


    % =============== plot ===============================
    ax1 = axList{ii};
    ax2 = axList2{ii};
    
    % ------------------ bA -----------------------------
    sObj = pcolor(ax1, xData, yData, zData{1});   hold(ax1, 'on');  box(ax1, 'on');  grid(ax1, 'on');
    % sObj = mesh(ax, xMat, yMat, zMat);   hold(ax, 'on');  box(ax, 'on');  grid(ax, 'on');
    set(sObj,  'FaceColor', 'interp', 'EdgeColor', 'none');
    contour(ax1, xMat, yMat, zMat{1}, 'LineColor', 0.3*[1 1 1]);
    caxis(ax1, max(abs(caxis(ax1)))*[-1,+1]);   % colorbar 上下对称
    axes(ax1);    % make ax become gca
    crameri('vik', 'pivot', 0);
    cObj = colorbar(ax1);
    title(cObj, '(pT)', 'interpreter', 'latex');
    
    % Phi_m = Phi* curve
    w_L_list = [min(yData(:)), max(yData(:))];
    w_L_list = linspace(w_L_list(1), w_L_list(2), 500);
    P_star = dA.I_star*2*pi*(w_L_list*L).^2;
    plot(ax1, 10*1e3*P_star, w_L_list, '--', 'Color', colorList(7,:));
    plot(ax1, 20*1e3*P_star, w_L_list, '--', 'Color', colorList(8,:));

    title(ax1, sprintf('$L=%.1f~{\\rm cm}$ ($L/L_0 \\approx %.1f$)', L, L/L0),'interpreter','latex');
    


    % --------------------- delta_bA ---------------------------------------------------------------
    sObj = pcolor(ax2, xMat, yMat, zMat{2});   hold(ax2, 'on');  box(ax2, 'on');  grid(ax2, 'on');
    set(sObj,  'FaceColor', 'interp', 'EdgeColor', 'none');
    contour(ax2, xMat, yMat, zMat{2}, contour_levels, 'LineColor', 0.3*[1 1 1]);
    
    % clim(ax, max(abs(caxis(ax)))*[-1,+1]);
    % clim(ax, [1e-1,1e3]);
    axes(ax2);    % make ax2 become gca
    crameri('-roma', 'pivot', 0);
    cObj = colorbar(ax2);
    title(cObj, '(fT)', 'interpreter', 'latex');
    set(ax2, 'ColorScale', 'log');
    
    % average polarization curve
    contour(ax2, xData, yData, avg_pol_sel, [0.65, 0.8, 0.9, 0.95], '-.', 'LineColor', colorList(6,:));
    
    
    if abs(L-0.3)<1e-6
        xticks([ax2, ax1], [0:5:15]);
        set(ax2, 'clim', [1e-1, 1e2]);
    elseif abs(L-0.5)<1e-6
        xticks([ax2, ax1], [0:10:50]);
        set(ax2, 'clim', [1e-1, 1e2]);
    elseif abs(L-0.8)<1e-6
        xticks([ax2, ax1], [0:20:100]);
        set(ax2, 'clim', [1e-1, 2e3]);
        set(cObj, 'ticks', 10.^[-1:4]);
    end
    yticks([ax2, ax1], 0.2:0.1:1.0);
     
    xlabel(ax2, 'Pump Beam Power $P_{\rm in}$ $\rm (mW)$', 'interpreter', 'latex');
    if ii == 3
        ylabel([ax2, ax1], '$w/ L$', 'interpreter', 'latex');
    end
    set([ax2, ax1],'Layer','top'); % to make grid line visiable
    
end


% set font size
allAxes = hf.findobj('Type','axes');
allLegend = hf.findobj('Type','legend');
allLine = hf.findobj('Type','line');
allErrorbar = hf.findobj('Type','errorbar');
allTextbox = hf.findobj('Type','textbox');
set([allAxes;allLegend;allTextbox],'fontsize',font_sz);
set([allLine;allErrorbar],'LineWidth',lw, 'MarkerSize', mkr_size);

%
annotation(hf, 'textbox', [0.00  0.94 0.05 0.05], 'String', '(a)', 'FontSize', 8, 'EdgeColor', 'none');
annotation(hf, 'textbox', [0.33 0.94 0.05 0.05], 'String', '(b)', 'FontSize', 8, 'EdgeColor', 'none');
annotation(hf, 'textbox', [0.655  0.94 0.05 0.05], 'String', '(c)', 'FontSize', 8, 'EdgeColor', 'none');
annotation(hf, 'textbox', [0.00  0.46 0.05 0.05], 'String', '(d)', 'FontSize', 8, 'EdgeColor', 'none');
annotation(hf, 'textbox', [0.33 0.46 0.05 0.05], 'String', '(e)', 'FontSize', 8, 'EdgeColor', 'none');
annotation(hf, 'textbox', [0.655  0.46 0.05 0.05], 'String', '(f)', 'FontSize', 8, 'EdgeColor', 'none');


%
annotation('textbox', [0.46177355664488 0.838841931721518  0.05 0.05], 'String', '$b_{\rm A}$', 'FontSize', 8, 'EdgeColor', 'none', 'interpreter', 'latex');   % bA
annotation('textbox', [0.511770152505449 0.373057540018933,0.05,0.05], 'String', '$\delta b_{\rm A}$', 'FontSize', 8, 'EdgeColor', 'none', 'interpreter', 'latex');  % \delta bA

%% Export figure
print(hf, ['Fig8'], '-dpng', '-r600');
hf.PaperSize = hf.Position(3:4);
print(hf, ['Fig8'], '-dpdf', '-painters');

