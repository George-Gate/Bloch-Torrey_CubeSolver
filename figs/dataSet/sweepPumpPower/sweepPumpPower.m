solver = Table0.CubeSolver2();

g129 = -2*pi*11.86015e-3;  % rad/s/nT
g131 = +2*pi*3.515769e-3;  % rad/s/nT

% define cell parameters
cellPars = solver.cellPars_example;    % an example of cellPars structure is defined in the constructor of @CubeSolver2
cellPars.L = [0.8, 0.8, 0.8];
cellPars.D = 0.2;
cellPars.Temp = 110;
cellPars.OD = 6;
cellPars.Rrel = 1/100e-6; 

cellPars_129 = cellPars;
cellPars_129.lambda = 1e-3*[1,1,1];
cellPars_129.gXe = g129;
cellPars_129.G2c = 1/20;      % 1/s

cellPars_131 = cellPars;
cellPars_131.lambda = 2e-2*[1,1,1];
cellPars_131.gXe = g131;
cellPars_131.G2c = 1/20;      % 1/s

pumpBeamProfile = struct();
pumpBeamProfile.w = 0.4*cellPars.L(1);   % cm
pumpBeamProfile.xc = 0.02*cellPars.L(1);
pumpBeamProfile.yc = 0.12*cellPars.L(1);
circleMask = @(x,y,r) 1./((1+((x.*x+y.*y)/r^2).^20));  % 圆形通光孔

I_max_list = unique([linspace(0, 0.15, 1*64), linspace(0, 0.05, 1*64)]);
% I_max_list = [0.3818, 0.3909];

fileName = sprintf('sweepPumpPower %s', datestr(now(),'yyyymmdd_HHMMSS'));


aperture = 5.0;
pumpBeamProfile.aperture_mask = @(x,y)circleMask(x,y,aperture/2);           % 通光孔外形，通光区域为1，不通光区域为0
result = solver.get_pump_power_dependence(I_max_list, pumpBeamProfile, cellPars_129, cellPars_131, false);


save([fileName '.mat'], 'result', 'cellPars_129', 'cellPars_131', 'aperture', 'fileName');


%%
hf = figure('Position', [100 84 1021 868]);
ax1 = subplot(2,2,1);
ax2 = subplot(2,2,3);
ax3 = subplot(2,2,2);
ax4 = subplot(2,2,4);

g129 = cellPars_129.gXe;
g131 = cellPars_131.gXe;
R0 = g129/g131;
re = result;
plot(ax1, re.I_max_list, re.T2_129, 'DisplayName', sprintf('%.1f', aperture));  hold(ax1, 'on'); grid(ax1, 'on');
plot(ax1, re.I_max_list, -1./real(re.ptb_129.eig0), '--'); 
plot(ax2, re.I_max_list, re.T2_131, 'DisplayName', sprintf('%.1f', aperture));  hold(ax2, 'on'); grid(ax2, 'on');
plot(ax2, re.I_max_list, -1./real(re.ptb_131.eig0), '--');  
plot(ax3, re.I_max_list, re.f_129, 'DisplayName', sprintf('%.1f', aperture)); hold(ax3, 'on'); grid(ax3, 'on');
plot(ax3, re.I_max_list, abs(imag(re.ptb_129.eig0))/2/pi, '--');  
plot(ax4, re.I_max_list, -cellPars_129.B0*( abs(re.R) - abs(R0))/(R0), 'b-','DisplayName', 'Numerical Diagonalization');  hold(ax4, 'on'); grid(ax4, 'on');
plot(ax4, re.I_max_list, -cellPars_129.B0*( abs(re.ptb_R) - abs(R0) )/(R0), 'r--', 'DisplayName', 'Perturbation Approximation');
plot(ax4, re.I_max_list, -1i*(R0*re.ptb_131.order1st - re.ptb_129.order1st)/g131/(R0), 'm-.', 'DisplayName', '1st order contribution'); 
plot(ax4, re.I_max_list, -1i*(R0*re.ptb_131.order3rd_a - re.ptb_129.order3rd_a)/g131/(R0)-1i*(R0*re.ptb_131.order3rd_b - re.ptb_129.order3rd_b)/g131/(R0), 'k-.', 'DisplayName', '3rd order contribution'); 


xlabel(ax1, 'Center intensity $I_0$ $\rm (W/cm^2)$', 'interpreter', 'latex');
ylabel(ax1, '$T_{2,129}$ (s)', 'interpreter', 'latex');
title(ax1, re.type, 'interpreter', 'none');

xlabel(ax2, 'Center intensity $I_0$ $\rm (W/cm^2)$', 'interpreter', 'latex');
ylabel(ax2, '$T_{2,131}$ (s)', 'interpreter', 'latex');
title(ax2, sprintf( 'Num. Basis: %d, $\\lambda_{129}=%g$, $\\lambda_{131}=%g$' ,  re.numBasis, cellPars_129.lambda(1), cellPars_131.lambda(1)), 'interpreter', 'latex');   % basis info

xlabel(ax3, 'Center intensity $I_0$ $\rm (W/cm^2)$', 'interpreter', 'latex');
ylabel(ax3, '$f_{129}$ (Hz)', 'interpreter', 'latex');
title(ax3, sprintf( 'Pump Beam: $w/L$=%g, center=$(%g,%g)$cm' ,  re.pumpBeamProfile.w/cellPars_129.L(1), re.pumpBeamProfile.xc, re.pumpBeamProfile.yc ), 'interpreter', 'latex');   % pump beam info

xlabel(ax4, 'Center intensity $I_0$ $\rm (W/cm^2)$', 'interpreter', 'latex');
ylabel(ax4, '$-B_0 \cdot \delta R / R_0$ (nT)', 'interpreter', 'latex');
title(ax4, sprintf( 'Temp.: %g${\\rm ^\\circ C}$, OD=%g, $D=%g~{\\rm cm^2/s}$' ,  cellPars_129.Temp, cellPars_129.OD, cellPars_129.D ), 'interpreter', 'latex');   % pump beam info

lg = legend(ax4, {}, 'location', 'best');

print(hf, [fileName '.png'], '-dpng', '-r330');
