% Demo Code: sweep cell temperature, compute the differential field (bA) vs. cell temperature curve
% D, lambda, nRb都随温度变化，温度依赖模型见张祥栋博士论文3.5节
%    (#Translate# D, lambda, and nRb all change with cell temperature. For the temperature-dependent models, refer to Section 3.5 of Dr. Zhang Xiangdong's doctoral dissertation.)
% ## lambda is proportional to L ##
solver = Table0.CubeSolver2();

g129 = -2*pi*11.86015e-3;  % rad/s/nT
g131 = +2*pi*3.515769e-3;  % rad/s/nT

% define cell parameters
cellPars = solver.cellPars_example;    % an example of cellPars structure is defined in the constructor of @CubeSolver2
cellPars.L = 0.8*[1, 1, 1];
cellPars.D = 0.2;
cellPars.D_alpha = 1;    % 扩散系数的温度依赖，D \propto T^alpha   (#Translate# Temperature dependence of the diffusion coefficient, D \propto T^alpha)
cellPars.OD = 6/0.8*cellPars.L(1); 
cellPars.Temp = 110;
cellPars.Rrel = 1/(14e-6);

cellPars_129 = cellPars;
cellPars_129.lambda = 1e-3*[1,1,1]/0.8*cellPars.L(1);   % lambda is proportional to L
cellPars_129.T_lambda = 0.10 * 1.602e-19/1.38e-23;      % lambda的特征温度（对应Ebar）   (#Translate# Characteristic temperature of lambda (corresponding to Ebar).)
cellPars_129.gXe = g129;
cellPars_129.G2c = 1/20;      % 1/s

cellPars_131 = cellPars;
cellPars_131.lambda = 2e-2*[1,1,1]/0.8*cellPars.L(1);   % lambda is proportional to L
cellPars_131.T_lambda = 0.15 * 1.602e-19/1.38e-23;      % lambda的特征温度（对应Ebar）   (#Translate# Characteristic temperature of lambda (corresponding to Ebar).)
cellPars_131.gXe = g131;
cellPars_131.G2c = 1/20;      % 1/s

pumpBeamProfile = struct();
pumpBeamProfile.w = 0.3*cellPars.L(1);   % cm
pumpBeamProfile.xc = 0.02;
pumpBeamProfile.yc = 0.10;
pumpBeamProfile.I_max = 0.5;
circleMask = @(x,y,r) 1./((1+((x.*x+y.*y)/r^2).^20));  % 圆形通光孔    (#Translate# Circular aperture)

cellTemp_list = unique([linspace(100, 115, 2*56)]);

fileName = sprintf('sweepCellTemperature_v2 %s', datestr(now(),'yyyymmdd_HHMMSS'));


aperture = 10.0;
pumpBeamProfile.aperture_mask = @(x,y)circleMask(x,y,aperture/2);           % 通光孔外形，通光区域为1，不通光区域为0   (#Translate# Aperture shape function. Returns 1 in the transmissive region and 0 in the blocking region.)
result = solver.get_cell_temperature_dependence_v2(cellTemp_list, pumpBeamProfile, cellPars_129, cellPars_131, false);


%save([fileName '.mat'], 'result', 'cellPars_129', 'cellPars_131', 'aperture', 'fileName');


%% plot
hf = figure('Position', [2258 84 1021 868]);
ax1 = subplot(2,2,1);
ax2 = subplot(2,2,3);
ax3 = subplot(2,2,2);
ax4 = subplot(2,2,4);

g129 = cellPars_129.gXe;
g131 = cellPars_131.gXe;
R0 = g129/g131;
re = result;
plot(ax1, re.cellTemp_list, re.T2_129, 'DisplayName', sprintf('%.1f', aperture));  hold(ax1, 'on'); grid(ax1, 'on');
plot(ax1, re.cellTemp_list, -1./real(re.ptb_129.eig0), '--'); 
plot(ax2, re.cellTemp_list, re.T2_131, 'DisplayName', sprintf('%.1f', aperture));  hold(ax2, 'on'); grid(ax2, 'on');
plot(ax2, re.cellTemp_list, -1./real(re.ptb_131.eig0), '--');  
plot(ax3, re.cellTemp_list, re.f_129, 'DisplayName', sprintf('%.1f', aperture)); hold(ax3, 'on'); grid(ax3, 'on');
plot(ax3, re.cellTemp_list, abs(imag(re.ptb_129.eig0))/2/pi, '--');  
plot(ax4, re.cellTemp_list, -cellPars_129.B0*( abs(re.R) - abs(R0))/(R0), 'b-','DisplayName', 'Numerical Diagonalization');  hold(ax4, 'on'); grid(ax4, 'on');
plot(ax4, re.cellTemp_list, -cellPars_129.B0*( abs(re.ptb_R) - abs(R0) )/(R0), 'r--', 'DisplayName', 'Perturbation Approximation');
plot(ax4, re.cellTemp_list, -1i*(R0*re.ptb_131.order1st - re.ptb_129.order1st)/g131/(R0), 'm-.', 'DisplayName', '1st order contribution'); 
plot(ax4, re.cellTemp_list, -1i*(R0*re.ptb_131.order3rd_a - re.ptb_129.order3rd_a)/g131/(R0)-1i*(R0*re.ptb_131.order3rd_b - re.ptb_129.order3rd_b)/g131/(R0), 'k-.', 'DisplayName', '3rd order contribution'); 


xlabel(ax1, 'Cell Temperature $\rm (^\circ C)$', 'interpreter', 'latex');
ylabel(ax1, '$T_{2,129}$ (s)', 'interpreter', 'latex');
title(ax1, re.type, 'interpreter', 'none');

xlabel(ax2, 'Cell Temperature $\rm (^\circ C)$', 'interpreter', 'latex');
ylabel(ax2, '$T_{2,131}$ (s)', 'interpreter', 'latex');
title(ax2, sprintf( 'Num. Basis: %d, $\\lambda_{129}=%g$, $\\lambda_{131}=%g$' ,  re.numBasis, cellPars_129.lambda(1), cellPars_131.lambda(1)), 'interpreter', 'latex');   % basis info

xlabel(ax3, 'Cell Temperature $\rm (^\circ C)$', 'interpreter', 'latex');
ylabel(ax3, '$f_{129}$ (Hz)', 'interpreter', 'latex');
% title(ax3, sprintf( 'Pump Beam: $w/L$=%g, center=$(%g,%g)$cm' ,  re.pumpBeamProfile.w/cellPars_129.L(1), re.pumpBeamProfile.xc, re.pumpBeamProfile.yc ), 'interpreter', 'latex');   % pump beam info
title(ax3, sprintf( 'Pump Beam: $w/L$=%g, center=$(%g,%g)$cm \n$I_{\\rm max}=%g~{\\rm W/cm^2}$, $P_{\\rm beam}=%g~{\\rm mW}$' ,  ...
                             re.pumpBeamProfile.w/cellPars_129.L(1), ...
                             re.pumpBeamProfile.xc, re.pumpBeamProfile.yc, ...
                             re.pumpBeamProfile.I_max, 1e3*re.pumpBeamProfile.P_beam ), 'interpreter', 'latex');   % pump beam info

xlabel(ax4, 'Cell Temperature $\rm (^\circ C)$', 'interpreter', 'latex');
ylabel(ax4, '$b_{\rm A}$ (nT)', 'interpreter', 'latex');
title(ax4, sprintf( 'OD($%g~{\\rm ^\\circ C}$) = %g, $D=%g~{\\rm cm^2/s}$' ,  cellPars_129.Temp, cellPars_129.OD, cellPars_129.D ), 'interpreter', 'latex');   % pump beam info

lg = legend(ax4, {}, 'location', 'best');

%print(hf, [fileName '.png'], '-dpng', '-r330');
