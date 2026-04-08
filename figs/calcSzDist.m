clear;

OD = 6;
Temp = 110;   % degC
w_L = 0.4;   % width of laser beam

kappa = 493;
gS = 2;
mu0 = 4*pi*1e-7;       % N/A^2     % 1 T = 1 N/(m*A)
muB = 927.400968e-26;  % J/T
nRb = getNumberDensity2(Temp);   % cm^-3

Bm = -1/3*kappa*nRb*1e6*gS*mu0*muB* 1e9;   % in unit of nT

fprintf('OD: %g, Temperature: %gdegC, n_Rb: %.3E cm^-3\n', OD, Temp, nRb);
fprintf('w/L: %g\n', w_L);

%% 1D plot (alpha_z)
z_L_list = linspace(-0.5 ,+0.5, 2000);
alpha0_list = [0.1, 1, 3, 10, 50, 100]';

alpha_z = lambertWpropagation(alpha0_list, OD, z_L_list+0.5);


sol_alpha1D.z_L_list = z_L_list;
sol_alpha1D.alpha0_list = alpha0_list;
sol_alpha1D.alpha_z = alpha_z;

figure();
plot(z_L_list, alpha_z./alpha0_list);
title(['OD = ', num2str(OD)]);
lg = legend(num2str([sol_alpha1D.alpha0_list]), 'Position',[0.371118573698493 0.293798521101819 0.0941176464890732 0.253529405225726]);
title(lg, '$\Phi_0 / \Phi^*$', 'interpreter', 'latex');
hold on;
% 强光极限 线性近似
plot(z_L_list', (1-OD./(1+sol_alpha1D.alpha0_list).*(z_L_list+0.5))', '--');
ylim([0,1]);

%% 2D plot (B_A, y=0 plane)
max_alpha0 = [3, 5, 10, 50, 100];
[z_L, x_L] = meshgrid(linspace(-0.5, 0.5, 400), linspace(0, 0.5*sqrt(2), 200));
for ii = 1:length(max_alpha0)
    alpha0 = max_alpha0(ii) * exp( - x_L(:,1).^2 ./ (2*w_L^2)  );   % 假设光斑中心在 (0,0), 入射中心的 alpha0 为 max_alpha0
    alpha_z = lambertWpropagation(alpha0, OD, z_L+0.5);
    Sz = 0.5*alpha_z./(alpha_z+1);
    B_A = -2/3*kappa*gS*nRb*1e6*mu0*muB*Sz * 1e9;   % in unit of nT
    
    sol_BA_2D(ii).max_alpha0 = max_alpha0(ii);
    sol_BA_2D(ii).z_L = z_L;
    sol_BA_2D(ii).x_L = x_L;
    sol_BA_2D(ii).Sz = Sz;
    sol_BA_2D(ii).B_A = B_A;
    
    figure('Position', [710 510 580 248]);
    sObj = pcolor(z_L, x_L, abs(B_A));   hold on;
    set(sObj,  'FaceColor', 'interp', 'EdgeColor', 'none');
    contour(z_L, x_L, abs(B_A), 'LineColor', 0.2*[1 1 1]);
    xlabel('z/L');      ylabel('sqrt(x^2+y^2)/L');
    box on; grid on;
    title(['\alpha_0(0) = ' num2str(max_alpha0(ii))]);
    set(gca,'Layer','top');
    caxis([0,45]);
end
%%
getNumberDensity(110)

%%
% calc number density
% use Killian formula, see https://csrc.yuque.com/sensing/basics/plmf5u
function n = getNumberDensity(temp_degc)
    % kB_val = kB();    % function in OpticalPumping\const
    kB_val = 1.3806488e-23;%2010 CODATA: http://physics.nist.gov/cuu/Constants/ unit: J K-1
    T = temp_degc + 273.15;
    p = 10.^(  9.55 - 4132./T  );   % in unit of Pa
    n = p./(kB_val*T)/1e6;     % in unit of cm^-3
end

% calc number density
% use Daniel Adam Steck's Rubidium data, see https://csrc.yuque.com/sensing/basics/plmf5u
function n = getNumberDensity2(temp_degc)
    % kB_val = kB();    % function in OpticalPumping\const
    kB_val = 1.3806488e-23;%2010 CODATA: http://physics.nist.gov/cuu/Constants/ unit: J K-1
    T = temp_degc + 273.15;
    p = 10.^(  9.318 - 4040./T  );   % in unit of Pa
    n = p./(kB_val*T)/1e6;     % in unit of cm^-3
end

function alpha_z = lambertWpropagation(alpha0, od, z_L)
    alpha_z = lambertw(alpha0.*exp(alpha0-od.*z_L));
end