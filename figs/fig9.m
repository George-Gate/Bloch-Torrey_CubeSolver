clear; clc;

load(['dataSet' filesep 'total_drifts' filesep 'preparedData.mat']);

uP = 1e-3; % 功率相对涨落
deltaT = 10e-3;  % 温度涨落
delta_xc = 1e-4;  % xc涨落, cm
u_sigma = 1e-3; % sigma_abs的相对涨落
u_Gamma = 2e-3; % Rrel的相对涨落


%% plot
hf = figure('Unit', 'Centimeters', 'Position', [1 1 18 22]);
mkr_size = 2; lw = 0.75;
font_sz = 7;

HH = 0.205;
x1 = 0.01;  x2 = 0.34;  x3 = 0.67;

y5 = 0.005;  y1 = 0.79;
dH = (y1-y5)/((4*HH+0.01)/HH);

y5 = y5;  y4 = y5+dH*(HH+0.01)/HH;  y3 = y4+dH; y2 = y3+dH; y1 = y1;


ax11 =  axes('OuterPosition', [x1, y1, 0.37, HH+0.002]); 
ax21 =  axes('OuterPosition', [x2, y1, 0.36, HH+0.002]); 
ax31 =  axes('OuterPosition', [x3, y1, 0.36, HH+0.002]);

ax12 =  axes('OuterPosition', [x1, y2, 0.37, HH]); 
ax22 =  axes('OuterPosition', [x2, y2, 0.36, HH]); 
ax32 =  axes('OuterPosition', [x3, y2, 0.36, HH]);

ax13 =  axes('OuterPosition', [x1, y3, 0.37, HH]); 
ax23 =  axes('OuterPosition', [x2, y3, 0.36, HH]); 
ax33 =  axes('OuterPosition', [x3, y3, 0.36, HH]);

ax14 =  axes('OuterPosition', [x1, y4, 0.37, HH]); 
ax24 =  axes('OuterPosition', [x2, y4, 0.36, HH]); 
ax34 =  axes('OuterPosition', [x3, y4, 0.36, HH]);

ax15 =  axes('OuterPosition', [x1, y5, 0.37, HH+0.01]); 
ax25 =  axes('OuterPosition', [x2, y5, 0.36, HH+0.01]); 
ax35 =  axes('OuterPosition', [x3, y5, 0.36, HH+0.01]);

axList = [ax11, ax21, ax31; 
          ax12, ax22, ax32; 
          ax13, ax23, ax33;
          ax14, ax24, ax34;
          ax15, ax25, ax35;];

colorList = colororder();
colorList(8,:) = [51 204 51]/255;


% characteristic cell length
d = dataList{1}{1};
Bm_at_Twork = d.result{1, 2}.characteristic_info.Bm;
p129 = d.cellPars_129;
p131 = d.cellPars_131;
L0 = (  abs(p131.lambda(1) - p129.lambda(1))*pi^4*p129.D^2 / (p129.L(1)*abs(p129.gXe*Bm_at_Twork)^2)  )^(1/3);


for iL = 1:length(dataList)
    dA = dataList{iL}{1};
    dB = dataList{iL}{2};


    xData = {1e3*dA.P_beam_mat{2};
             1e3*dA.P_beam_mat{2};
             1e3*dB.P_beam_mat;
             1e3*dB.P_beam_mat;
             1e3*dB.P_beam_mat;};
    yData = {dA.w_L_mat{2};
             dA.w_L_mat{2};
             dB.w_L_mat;
             dB.w_L_mat;
             dB.w_L_mat;};
    zData = {1e6*deltaT*dA.dbA_dT_mat; ...                 % dbA / dT  , in unit of fT 
             1e6*uP*dA.P_beam_mat{2}.*dA.dbA_dPeam_mat; ... % dbA / dP
             1e6*delta_xc*dB.dbA_dxc_mat; ...              % dbA / dxc
             1e6*u_sigma*dB.rel_dbA_dsigma_mat;  ...       % dbA / dsigma
             1e6*u_Gamma*dB.rel_dbA_dRrel_mat};            % dbA / dRrel
    I_star = [dA.I_star, dA.I_star, dB.I_star, dB.I_star, dB.I_star];

    
    L=dA.cellPars_129.L(1);

    for ii = 1:length(zData)
        ax = axList(ii, length(dataList) - iL + 1);
        
        
        % 重新构建xData, yData矩阵，以消除数值误差。但需要校验数值误差不大于1e-10
        [xData2, yData2] = ndgrid(xData{ii}(:,1), yData{ii}(1,:));
        assert(max(max(abs(xData2-xData{ii}))) < 1e-10);
        assert(max(max(abs(yData2-yData{ii}))) < 1e-10);
        % interpret data to make plot more smooth
        F = griddedInterpolant(xData2, yData2, zData{ii}, 'spline');
        [xMat, yMat] = meshgrid( linspace(min(xData{ii}(:)),max(xData{ii}(:)),500), linspace(min(yData{ii}(:)),max(yData{ii}(:)),500)   );
        zMat = F(xMat, yMat);

        sObj = pcolor(ax, xData{ii}, yData{ii}, zData{ii});   hold(ax, 'on');  box(ax, 'on');  grid(ax, 'on');
        % sObj = mesh(ax, xMat, yMat, zMat);   hold(ax, 'on');  box(ax, 'on');  grid(ax, 'on');
        set(sObj,  'FaceColor', 'interp', 'EdgeColor', 'none');
        contour(ax, xMat, yMat, zMat, 'LineColor', 0.3*[1 1 1]);
        caxis(ax, max(abs(caxis(ax)))*[-1,+1]);   % colorbar 上下对称
        axes(ax);    % make ax become gca
        crameri('vik', 'pivot', 0);
        
        % Phi_m = Phi* curve
        w_L_list = [min(yData{ii}(:)), max(yData{ii}(:))];
        w_L_list = linspace(w_L_list(1), w_L_list(2), 500);
        P_star = I_star(ii)*2*pi*(w_L_list*L).^2;
        plot(ax, 10*1e3*P_star, w_L_list, '--', 'Color', colorList(7,:));
        plot(ax, 20*1e3*P_star, w_L_list, '--', 'Color', colorList(8,:));
        
        if iL==3
            xticks(ax, [0:5:15]);
        elseif iL==2
            xticks(ax, [0:10:50]);
        elseif iL==1
            xticks(ax, [0:20:100]);
        end

        set(ax,'Layer','top'); % to make grid line visiable
        if ii == length(zData)
            xlabel(ax, 'Pump Beam Power $P_{\rm in}$ $\rm (mW)$', 'interpreter', 'latex');
        end
        if iL==3
            ylabel(ax, '$w/ L$', 'interpreter', 'latex');
        end
        if ii == 1
            title(ax, sprintf('$L=%.1f~{\\rm cm}$ ($L/L_0 \\approx %.1f$)', L, L/L0),'interpreter','latex');
        end
        yticks(ax, 0.2:0.1:1.0);
        cObj = colorbar(ax);
        title(cObj, '(fT)', 'interpreter', 'latex');
    end
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
%mark subplot label
xGrid = [0.00, 0.325, 0.655];
yGrid = [0.17, 0.36, 0.56, 0.75, 0.95];
legendList = {'(a)', '(b)', '(c)';
              '(d)', '(e)', '(f)';
              '(g)', '(h)', '(i)';
              '(j)', '(k)', '(l)';
              '(m)', '(n)', '(o)'};
for ix = 1:length(xGrid)
    for iy = length(yGrid):-1:1
        annotation(hf, 'textbox', [xGrid(ix) yGrid(iy) 0.05 0.05], 'String', legendList{end-iy+1, ix}, 'FontSize', 8, 'EdgeColor', 'none');
    end
end


% mark gradient name
annotation('textbox', [0.466185321350762 0.911576413959086  0.05 0.05], 'String', '$b_{\rm A}^{(T)}$', 'FontSize', 8, 'EdgeColor', 'none', 'interpreter', 'latex');   
annotation('textbox', [0.466181917211331 0.717593261131167,0.05,0.05], 'String', '$b_{\rm A}^{(P)}$', 'FontSize', 8, 'EdgeColor', 'none', 'interpreter', 'latex');  
annotation('textbox', [0.500005446623096 0.523983152827918,0.05,0.05], 'String', '$b_{\rm A}^{(x)}$', 'FontSize', 8, 'EdgeColor', 'none', 'interpreter', 'latex'); 
annotation('textbox', [0.46912309368192 0.332033694344164,0.05,0.05], 'String', '$b_{\rm A}^{(\nu)}$', 'FontSize', 8, 'EdgeColor', 'none', 'interpreter', 'latex');  
annotation('textbox', [0.498534858387802 0.0935258724428402,0.05,0.05], 'String', '$b_{\rm A}^{(\Gamma)}$', 'FontSize', 8, 'EdgeColor', 'none', 'interpreter', 'latex');  


%% Export figure
print(hf, ['Fig9'], '-dpng', '-r600');
hf.PaperSize = hf.Position(3:4);
print(hf, ['Fig9'], '-dpdf', '-painters');

