% lambda 正比于 L
solver = Table0.CubeSolver2();

g129 = -2*pi*11.86015e-3;  % rad/s/nT
g131 = +2*pi*3.515769e-3;  % rad/s/nT

% define cell parameters
cellPars = solver.cellPars_example;    % an example of cellPars structure is defined in the constructor of @CubeSolver2
cellPars.L = 0.2*[1, 1, 1];
cellPars.D = 0.2;
cellPars.D_alpha = 1;    % 扩散系数的温度依赖，D \propto T^alpha
cellPars.OD = 6/0.8*cellPars.L(1); 
cellPars.Temp = 110;
cellPars.Rrel = 1/(100e-6);

if abs( cellPars.L(1) - 0.2 ) < 1e-6
    P_max = 10e-3;
elseif abs( cellPars.L(1) - 0.3 ) < 1e-6
    P_max = 15e-3;
elseif abs( cellPars.L(1) - 0.5 ) < 1e-6
    P_max = 45e-3;
elseif abs( cellPars.L(1) - 0.8 ) < 1e-6
    P_max = 90e-3;
elseif abs( cellPars.L(1) - 1.3 ) < 1e-6
    P_max = 200e-3;
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

pumpBeamProfile = struct();
% pumpBeamProfile.w = 0.3*cellPars.L(1);   % cm
pumpBeamProfile.xc = 0.02*cellPars.L(1);
pumpBeamProfile.yc = 0.03*cellPars.L(1);
% pumpBeamProfile.xc = 0;
% pumpBeamProfile.yc = 0;
pumpBeamProfile.aperture = 10.0;
circleMask = @(x,y,r) 1./((1+((x.*x+y.*y)/r^2).^20));  % 圆形通光孔
pumpBeamProfile.aperture_mask = @(x,y)circleMask(x,y,pumpBeamProfile.aperture/2);           % 通光孔外形，通光区域为1，不通光区域为0

% I_max_list = unique([linspace(0, 7, 2*56), linspace(0, 1, 3*56)])/5;
% I_max_list = [30, 50, 300, 1, 2, 3];

fileName = sprintf('sweepPumpPower_varyBeamWidth_v2 L=%.2fcm %s', cellPars.L(1), datestr(now(),'yyyymmdd_HHMMSS'));


w_L_list = [0.25, 0.3, 0.35, 0.4, 0.5, 0.7, 0.8];
% aperture_list = 3.0;
result = {};
for ii = 1:length(w_L_list)
    pumpBeamProfile.w = w_L_list(ii)*cellPars.L(1);
    I_max_list = P_max/(2*pi*pumpBeamProfile.w^2)*unique([linspace(0, 0.3, 128), linspace(0, 1, 160)]);
    result{ii} = solver.get_pump_power_dependence(I_max_list, pumpBeamProfile, cellPars_129, cellPars_131, false);
end

save([fileName '.mat'], 'result', 'cellPars_129', 'cellPars_131', 'w_L_list', 'fileName');


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
    re = result{ii};
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