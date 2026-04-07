% plot the relation of a^(1) and a^(3) vs. OD and Phi_m

dataList = {
            'sweep_OD_and_Phim w_L=0.2 20250616_091401.mat';
            % 'sweep_OD_and_Phim w_L=0.3 20250616_164159.mat';
            'sweep_OD_and_Phim w_L=0.4 20250615_165748.mat';
            % 'sweep_OD_and_Phim w_L=0.5 20250616_195755.mat';
            % 'sweep_OD_and_Phim w_L=0.6 20250616_011325.mat';
            % 'sweep_OD_and_Phim w_L=0.8 20250616_030229.mat';
            'sweep_OD_and_Phim w_L=1.0 20250616_142155.mat';
            };
dataPath = ['dataSet' filesep 'sweep_OD_and_Phim'];

g129 = -2*pi*11.86015e-3;  % rad/s/nT
g131 = +2*pi*3.515769e-3;  % rad/s/nT
R0 = g129/g131;
colorList = colororder();

%%
hf = figure('Unit', 'Centimeters', 'Position', [5 5 18 10]);
mkr_size = 2; lw = 0.75;
font_sz = 7;


x0 = 0.05;
y0 = 0.08;
dx = 0.305; dy = 0.47;
HH = 0.41;
WW = 0.25;

x1 = x0;  x2 = x1+dx;  x3 = x2+dx;
y1 = y0+dy;  y2 = y0;

axList{1,1} =  axes('InnerPosition', [x1, y1, WW, HH]);    hold on; box on; grid on; ylabel('OD');
axList{2,1} =  axes('InnerPosition', [x2, y1, WW, HH]);    hold on; box on; grid on;   
axList{3,1} =  axes('InnerPosition', [x3, y1, WW, HH]);    hold on; box on; grid on;
axList{1,2} =  axes('InnerPosition', [x1, y2, WW, HH]);    hold on; box on; grid on; ylabel('OD');  xlabel('$\Phi_{\rm m}/\Phi^*$', 'Interpreter','latex');
axList{2,2} =  axes('InnerPosition', [x2, y2, WW, HH]);    hold on; box on; grid on; xlabel('$\Phi_{\rm m}/\Phi^*$', 'Interpreter','latex');
axList{3,2} =  axes('InnerPosition', [x3, y2, WW, HH]);    hold on; box on; grid on; xlabel('$\Phi_{\rm m}/\Phi^*$', 'Interpreter','latex');


for iD = 1:3
    
    load([dataPath filesep dataList{iD}]);

    a1 = {}; a3 = {};
    OD_list = {};
    
    lambda129 = cellPars_129.lambda(1);
    lambda131 = cellPars_131.lambda(1);
    L = cellPars_129.L(1);
    D = cellPars_129.D;
    for ii = 1:length(Im_Istar_list)
        ptb_129 = result{ii}.ptb_129;
        ptb_131 = result{ii}.ptb_131;
        Bm = result{ii}.polarizationFieldInfo.Bm;
        OD_list{ii} = result{ii}.OD_list;
    
        bA_1st = +1i*(ptb_129.order1st/g129 - ptb_131.order1st/g131);
        bA_3rd = +1i*((ptb_129.order3rd_a+ptb_129.order3rd_b)/g129 - (ptb_131.order3rd_a+ptb_131.order3rd_b)/g131);
        a1{ii} = bA_1st./Bm/(lambda131-lambda129);
        a3{ii} = bA_3rd./(R0^2-1)/g131^2/L^4*D^2*pi^4./Bm.^3;
    end
    
    % 组装数据
    [OD_mat, Phi_rel_mat] = ndgrid( OD_list{1}, Im_Istar_list);    % Phi_m / Phi_star
    a1_mat = cell2mat(a1);
    a3_mat = cell2mat(a3);
    w_L = result{1}.pumpBeamProfile.w / L;
    
    % plot 
    
    xData = Phi_rel_mat;
    yData = OD_mat;
    zData = {a1_mat, a3_mat};
    
    fprintf('min(|a1|)=%.3g, max(|a1|)=%.3g\n', min(min(abs(a1_mat(1:50,50:end)))), max(abs(a1_mat(:))));
    fprintf('min(|a3|)=%.3g, max(|a3|)=%.3g\n', min(min(abs(a3_mat(1:50,50:end)))), max(abs(a3_mat(:))));

    for ii = 1:2
        F = griddedInterpolant(xData', yData', zData{ii}', 'spline');
        [xMat, yMat] = meshgrid( linspace(min(xData(:)),max(xData(:)),300), linspace(min(yData(:)),max(yData(:)),300)   );
        zMat = F(xMat, yMat);
        ax = axList{iD,ii};
        sObj = pcolor(ax, xMat, yMat, abs(zMat));   
        set(sObj,  'FaceColor', 'interp', 'EdgeColor', 'none');
        
        % caxis(ax, max(abs(caxis(ax)))*[-1,+1]);   % colorbar 上下对称
        
        axes(ax);    % make ax become gca
        % set colormap
        % crameri('vik', 'pivot', 0);
        crameri('lapaz');
        
        set(ax, 'ColorScale', 'log');
        
        set(ax,'Layer','top');  % to make grid line visiable
        if ii==1
            title(sprintf('$w/L=%.1f$', w_L), 'Interpreter','latex');
            clim(ax, [1e-4, 1e-1]);   % colorbar upper/lower limits
            [~,hc] = contour(ax, xMat, yMat, log10(abs(zMat)), [-1,-2,-3,-4],'LineColor', 0.3*[1 1 1]);
            hc.ShowText = 'on';
            hc.TextList = [-2];
            hc.LabelFormat = '1E%+.0f';
            hc.LabelSpacing = 200;
            clabel([],hc,'FontSize',font_sz-1);
            
            % contour(ax, xMat, yMat, log10(abs(zMat)), [-1,-2],'LineColor', 0.3*[1 1 1], 'ShowText', 'on');
        else
            clim(ax, [1e-11, 1e-1]);   % colorbar upper/lower limits
            [~,hc] = contour(ax, xMat, yMat, log10(abs(zMat)), [-11:-1],'LineColor', 0.3*[1 1 1]);
            hc.ShowText = 'on';
            if iD == 1
                hc.TextList = [-3];
            else
                hc.TextList = [-2, -6];
            end
            hc.LabelFormat = '1E%+.0f';
            hc.LabelSpacing = 200;
            clabel([],hc,'FontSize',font_sz-1);
        end

        contour(ax, xMat, yMat, zMat, [-100 0 100], '--', 'LineColor', colorList(7,:));  % zero crossing line
    end

end

cObj = colorbar(axList{3,1});
set(cObj, 'Position', [0.94 0.55 0.02 0.41]);
cObj = colorbar(axList{3,2});
set(cObj, 'Position', [0.94 0.08 0.02 0.41],'Ticks',[1e-11, 1e-6, 1e-1]);

% figure mark
xGrid = [0.00, 0.31, 0.62];
yGrid = [0.47, 0.95];
legendList = {'(a)', '(b)', '(c)';
              '(d)', '(e)', '(f)'};
for ix = 1:length(xGrid)
    for iy = length(yGrid):-1:1
        annotation(hf, 'textbox', [xGrid(ix) yGrid(iy) 0.05 0.05], 'String', legendList{end-iy+1, ix}, 'FontSize', 8, 'EdgeColor', 'none');
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

%% Export figure
print(hf, 'Fig4.png', '-dpng', '-r600');
hf.PaperSize = hf.Position(3:4);
print(hf, 'Fig4.pdf', '-dpdf', '-vector');
