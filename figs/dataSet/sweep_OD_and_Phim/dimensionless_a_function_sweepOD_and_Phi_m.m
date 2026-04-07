% find the relation of a^(1) and a^(3) vs. OD and Phi_m

solver = Table0.CubeSolver2();

g129 = -2*pi*11.86015e-3;  % rad/s/nT
g131 = +2*pi*3.515769e-3;  % rad/s/nT

% define cell parameters
cellPars = solver.cellPars_example;    % an example of cellPars structure is defined in the constructor of @CubeSolver2
cellPars.L = [0.8, 0.8, 0.8];
cellPars.D = 0.2;
cellPars.OD = 6/0.8*cellPars.L(3);
cellPars.Rrel = 1/(100e-6);

cellPars_129 = cellPars;
cellPars_129.lambda = 1e-3*[1,1,1];
cellPars_129.gXe = g129;
cellPars_129.G2c = 1/20;      % 1/s

cellPars_131 = cellPars;
cellPars_131.lambda = 2e-2*[1,1,1];
cellPars_131.gXe = g131;
cellPars_131.G2c = 1/20;      % 1/s

w_L = 0.9;

pumpBeamProfile = struct();
pumpBeamProfile.w = w_L*cellPars.L(1);
pumpBeamProfile.xc = 0.02*cellPars.L(1);
pumpBeamProfile.yc = 0.03*cellPars.L(1);
% pumpBeamProfile.xc = 0;
% pumpBeamProfile.yc = 0;
pumpBeamProfile.aperture = 3.0;
circleMask = @(x,y,r) 1./((1+((x.*x+y.*y)/r^2).^20));  % 圆形通光孔
pumpBeamProfile.aperture_mask = @(x,y)circleMask(x,y,pumpBeamProfile.aperture/2);           % 通光孔外形，通光区域为1，不通光区域为0

nRb = solver.getRbNumberDensity(cellPars.Temp);  % cm^(-3)
sigma_abs = cellPars.OD/cellPars.L(3)/nRb;    % cm^2
Phi_star = cellPars.Rrel/sigma_abs;  % in unit of 1/(s*cm^2)
h = 6.62606957e-34;    % Plank's constant, J*s
I_star = Phi_star*h*cellPars.laserFreq;


% I_max_list = unique([linspace(0, 7, 2*56), linspace(0, 1, 3*56)])/5;
% I_max_list = [30, 50, 300, 1, 2, 3];

fileName = sprintf('sweep_OD_and_Phim w_L=%.1f %s', w_L, datestr(now(),'yyyymmdd_HHMMSS'));



Im_Istar_list = [0.5, linspace(1,100,100)];
cellTemp_list = - 4040./log10(linspace(10^(-4040/(80+273.15)), 10^(-4040/(140+273.15)), 128)) - 273.15;     % 80~140degC, OD linspace
% aperture_list = 3.0;
result = {};
for ii = 1:length(Im_Istar_list)
    tS = tic;
    pumpBeamProfile.I_max = Im_Istar_list(ii)*I_star;
    result{ii} = solver.get_cell_temperature_dependence(cellTemp_list, pumpBeamProfile, cellPars_129, cellPars_131, false);
    fprintf('Im/Istar=%.1f finished. Time cost: %.1fs\n', Im_Istar_list(ii), toc(tS));
end

save([fileName '.mat'], 'result', 'cellPars_129', 'cellPars_131', 'Im_Istar_list', 'I_star', 'cellTemp_list', 'fileName');

%% plot different Im/I* line in one figure
g129 = -2*pi*11.86015e-3;  % rad/s/nT
g131 = +2*pi*3.515769e-3;  % rad/s/nT
a1 = {}; a3 = {};
OD_list = {};
R0 = g129/g131;
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


% plot 
hf = figure( Position = [89 136 2281 1008]);
axList{1} = subplot(1,2,1); hold on; box on; grid on;
axList{2} = subplot(1,2,2); hold on; box on; grid on;

xData = Phi_rel_mat;
yData = OD_mat;
zData = {a1_mat, a3_mat};
name_str = {'$a^{(1)}$', '$a^{(3)}$'};
for ii = 1:2
    F = griddedInterpolant(xData', yData', zData{ii}', 'spline');
    [xMat, yMat] = meshgrid( linspace(min(xData(:)),max(xData(:)),2000), linspace(min(yData(:)),max(yData(:)),500)   );
    zMat = F(xMat, yMat);
    ax = axList{ii};
    sObj = pcolor(ax, xMat, yMat, abs(zMat));   hold(ax, 'on');  box(ax, 'on');  grid(ax, 'on');
    set(sObj,  'FaceColor', 'interp', 'EdgeColor', 'none');

    caxis_max = max(max(abs(zMat)));
    caxis_min = min(min(abs(zMat(1:end/2,end/2:end))));    % 右下角的最小值
    caxis_max = 10^ceil(log10(caxis_max));
    caxis_min = 10^floor(log10(caxis_min));

    contour(ax, xMat, yMat, log10(abs(zMat)),[log10(caxis_min):log10(caxis_max)], 'LineColor', 0.3*[1 1 1]);
    
    % caxis(ax, max(abs(caxis(ax)))*[-1,+1]);   % colorbar 上下对称
    caxis(ax, [caxis_min, caxis_max]);   % colorbar 上下限

    axes(ax);    % make ax become gca
    % set colormap
    % crameri('vik', 'pivot', 0);
    crameri('lapaz');
    
    set(ax, 'ColorScale', 'log');
    xlabel(ax, '$\Phi_{\rm m} / \Phi^*$', 'interpreter', 'latex');
    ylabel(ax, 'OD', 'interpreter', 'latex');
    cObj = colorbar(ax);
    title(ax, name_str{ii}, 'interpreter', 'latex');
    contour(ax, xMat, yMat, zMat, [-100 0 100], '--', 'LineColor', 'red');
end

set(hf, 'name',sprintf('w/L = %.1f', result{1}.pumpBeamProfile.w/L));
print(hf, [fileName '.png'], '-dpng', '-r300');
