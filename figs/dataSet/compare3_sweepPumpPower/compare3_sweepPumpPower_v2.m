% lambda 正比于 L
solver = Table0.CubeSolver2();

g129 = -2*pi*11.86015e-3;  % rad/s/nT
g131 = +2*pi*3.515769e-3;  % rad/s/nT

% define cell parameters
cellPars = solver.cellPars_example;    % an example of cellPars structure is defined in the constructor of @CubeSolver2
cellPars.L = 0.3*[1, 1, 1];
cellPars.D = 0.2;
cellPars.D_alpha = 1;    % 扩散系数的温度依赖，D \propto T^alpha
cellPars.OD = 6/0.8*cellPars.L(1); 
cellPars.Temp = 110;
cellPars.Rrel = 1/(100e-6);

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
pumpBeamProfile.w = 0.4*cellPars.L(1);   % cm
pumpBeamProfile.xc = 0.00;
pumpBeamProfile.yc = 0.00;
circleMask = @(x,y,r) 1./((1+((x.*x+y.*y)/r^2).^20));  % 圆形通光孔

I_star = 3.7414e-3;  % W/cm^2
I_max_list = unique([logspace(-3, -1, 3*64), linspace(0.1, 10, 3*64), logspace(1, 3, 3*64)])*I_star;
% I_max_list = [0.3818, 0.3909];

fileName = sprintf('sweepPumpPower_v2 %s', datestr(now(),'yyyymmdd_HHMMSS'));


aperture = 10.0;
pumpBeamProfile.aperture_mask = @(x,y)ones(size(x));
% pumpBeamProfile.aperture_mask = @(x,y)circleMask(x,y,aperture/2);           % 通光孔外形，通光区域为1，不通光区域为0
result = solver.get_pump_power_dependence(I_max_list, pumpBeamProfile, cellPars_129, cellPars_131, false);


save([fileName '.mat'], 'result', 'cellPars_129', 'cellPars_131', 'aperture', 'fileName');


%% calc. 近似解析式
gS = 2;
mu0 = 4*pi*1e-7;       % N/A^2     % 1 T = 1 N/(m*A)
muB = 927.400968e-26;  % J/T
h = 6.62606957e-34;    % Plank's constant, J*s
g129 = cellPars_129.gXe;
g131 = cellPars_131.gXe;
R0 = g129/g131;

re = result;
L = cellPars_129.L(1);
D = cellPars_129.D;
OD = cellPars_129.OD;
Bm = -cellPars_129.pumpPolarization/3*cellPars_129.kappa_RbXe*(1e6*re.characteristic_info.nRb)*gS*mu0*muB * 1e9;  % in unit of nT
I_star = re.characteristic_info.I_star;
base_ampl_weak = Bm*re.I_max_list/I_star;   % Bm*Phi_m/Phi_star,  in unit of nT
base_ampl_strong = Bm*I_star./re.I_max_list;   % Bm*Phi_star/Phi_m,  in unit of nT

cfs_weak = get_characteristic_functions_weakPump(re.pumpBeamProfile.w/L, OD);
cfs_strong_1 = get_characteristic_functions_strongPump(re.pumpBeamProfile.w/L);
cfs_strong_2 = get_characteristic_functions_strongPump(re.pumpBeamProfile.w/L/sqrt(2));
cfs_strong_3 = get_characteristic_functions_strongPump(re.pumpBeamProfile.w/L/sqrt(3));


bA1st_approx_weak = base_ampl_weak*(cellPars_129.lambda(1) - cellPars_131.lambda(1))*(cfs_weak.F1*cfs_weak.G1+cfs_weak.F2*cfs_weak.G2);
bA1st_approx_strong = base_ampl_strong*(cellPars_129.lambda(1) - cellPars_131.lambda(1))*(cfs_strong_1.F5);    % 强光极限近似
bA1st_approx_strong_2 = base_ampl_strong*(cellPars_129.lambda(1) - cellPars_131.lambda(1)).*(cfs_strong_1.F5 + I_star./re.I_max_list .* (OD/2 - 1).*cfs_strong_2.F5);    % 强光极限近似（2阶）
bA1st_approx_strong_3 = base_ampl_strong*(cellPars_129.lambda(1) - cellPars_131.lambda(1)).*(cfs_strong_1.F5 + I_star./re.I_max_list .* (OD/2 - 1).*cfs_strong_2.F5 ...
                                                                            + (I_star./re.I_max_list).^2.*(  (1+OD/6*(2*OD-9)).*cfs_strong_3.F5 - OD.^2/90.*cfs_strong_3.F7 ));    % 强光极限近似（3阶）


hf = figure('Position', [2258 342 1274 610]);
ax1 = subplot(1,2,1);
ax2 = subplot(1,2,2);

% plot(ax1, re.I_max_list, -cellPars_129.B0*( abs(re.R) - abs(R0))/(R0), 'b-','DisplayName', 'Numerical Diagonalization');  hold(ax4, 'on'); grid(ax4, 'on');
% plot(ax1, re.I_max_list, -cellPars_129.B0*( abs(re.ptb_R) - abs(R0) )/(R0), 'r--', 'DisplayName', 'Perturbation Approximation');

bA1st = -1i*(R0*re.ptb_131.order1st - re.ptb_129.order1st)/g131/(R0);
bA3rd = -1i*(R0*re.ptb_131.order3rd_a - re.ptb_129.order3rd_a)/g131/(R0)-1i*(R0*re.ptb_131.order3rd_b - re.ptb_129.order3rd_b)/g131/(R0);

% x_data = 1./(re.I_max_list/I_star);
x_data = (re.I_max_list/I_star);

plot(ax1, x_data, bA1st, 'k-', 'DisplayName', '$b_{\rm A}^{(1)}$'); hold(ax1, 'on'); grid(ax1, 'on');
% plot(ax1,x_data, bA3rd, 'b-', 'DisplayName', '$b_{\rm A}^{(3)}$'); 
plot(ax1, x_data, bA1st_approx_weak, 'm--', 'DisplayName', '$b_{\rm A}^{(1)}$ weak approx.');
plot(ax1, x_data, bA1st_approx_strong, 'b--', 'DisplayName', '$b_{\rm A}^{(1)}$ strong approx.');
plot(ax1, x_data, bA1st_approx_strong_2, 'r--', 'DisplayName', '$b_{\rm A}^{(1)}$ strong approx. 2');
plot(ax1, x_data, bA1st_approx_strong_3, 'g--', 'DisplayName', '$b_{\rm A}^{(1)}$ strong approx. 3');
% plot(ax1, [1,1], minmax([bA1st(:);bA3rd(:)]'), 'k-.', 'DisplayName', '$I^* \equiv h \nu \Phi^*$');

%复制 ax1, 然后放大
tmp_pos = ax2.Position;
delete(ax2);
ax2 = copyobj(ax1, hf);
ax2.Position = tmp_pos;
% 将y轴数据取绝对值
lines = ax2.findobj('type','Line');
for ii = 1:length(lines)
    lines(ii).YData = abs(lines(ii).YData);
end


xlabel(ax1, 'Center intensity $I_0/I^*$', 'interpreter', 'latex');
ylabel(ax1, '$b_{\rm A}$ (nT)', 'interpreter', 'latex');
title(ax1, sprintf( 'Num. Basis: %d, $\\lambda_{129}=%g$, $\\lambda_{131}=%g$, $L=%g~{\\rm cm}$' ,  re.numBasis, cellPars_129.lambda(1), cellPars_131.lambda(1), cellPars_129.L(1)), 'interpreter', 'latex');   % basis info

lg = legend(ax1, {}, 'location', 'best', 'interpreter', 'latex');


xlabel(ax2, 'Center intensity $I_0/I^*$', 'interpreter', 'latex');
ylabel(ax2, '$|b_{\rm A}|$ (nT)', 'interpreter', 'latex');
titleStr = sprintf( 'Temp.: %g${\\rm ^\\circ C}$, OD=%g, $D=%g~{\\rm cm^2/s}$, $I^*=%.1f~{\\rm mW/cm^2}$' ,  cellPars_129.Temp, cellPars_129.OD, cellPars_129.D, re.characteristic_info.I_star*1e3 );
titleStr = sprintf( '%s\nPump Beam: $w/L$=%g, center=$(%g,%g)$cm' ,  titleStr, re.pumpBeamProfile.w/cellPars_129.L(1), re.pumpBeamProfile.xc, re.pumpBeamProfile.yc );
title(ax2, titleStr, 'interpreter', 'latex');   % pump beam info

ylim(ax1, minmax([bA1st(:);bA3rd(:)]'));
xlim(ax1, [0, 10]);
xlim(ax2, [0, 1000]);
set(ax2, 'xscale', 'log', 'yscale', 'log');


print(hf, [fileName '.png'], '-dpng', '-r330');
