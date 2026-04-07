function [B_A_chebfun3, auxInfo] = calc_polarization_field(obj, I0_fun, cellPars, aperture_mask, disp_check)
%CALC_POLARIZATIONFIELD 计算 极化场 分布
% I0_fun: 二维函数@(x,y), 给定入射面的光强分布，单位 W/cm^2
% cellPars: 指定气室参数，需要的项目见下方代码
% aperture_mask: 通光孔形状，二维函数@(x,y), 在可以通光的区域返回1，在不能通光的区域返回0
% disp_check: 是否进行精度检查，并输出检查结果

% (#Translate#
%     CALC_POLARIZATIONFIELD Calculate the Rb polarization field distribution.
%     I0_fun: 2D function handle @(x,y), specifies the incident laser intensity distribution on the input plane, in units of W/cm^2.
%     cellPars: A structure specifying the vapor cell parameters; required fields are listed in the code below.
%     aperture_mask: 2D function handle @(x,y) defining the aperture shape. Returns 1 in transmissive regions and 0 in blocking regions.
%     disp_check: Logical flag to control whether accuracy checks are performed and their results displayed.
% )

    if ~exist('disp_check', 'var')
        disp_check = true;
    end
    if ~exist('aperture_mask', 'var') || ~isa(aperture_mask, 'function_handle')
        aperture_mask = @(x,y)ones(size(x));
    end


    L = cellPars.L;        % cell domain [Lx, Ly, Lz], in unit of cm. 
                           % The domain of B_A_chebfun3 is [-Lx/2, +Lx/2] x [-Ly/2, +Ly/2] x [-Lz/2, +Lz/2]
                           % Light propogates from z=-Lz/2 to z=+Lz/2
    OD = cellPars.OD;
    T = cellPars.Temp;      % degC
    Rrel = cellPars.Rrel;   % in unit of 1/s
    nu = cellPars.laserFreq; % Hz
    kappa_RbXe = cellPars.kappa_RbXe;

    h = 6.62606957e-34;    % Plank's constant, J*s
    gS = 2;
    mu0 = 4*pi*1e-7;       % N/A^2     % 1 T = 1 N/(m*A)
    muB = 927.400968e-26;  % J/T
    nRb = obj.getRbNumberDensity(T);  % cm^(-3)
    sigma_abs = OD/L(3)/nRb;    % cm^2
    s = cellPars.pumpPolarization;     % pump polarization

    alpha0 = @(x,y)I0_fun(x,y)/(h*nu*Rrel/sigma_abs);  % dimensionless
    alpha_z = @(x,y,z)lambertWpropagation(alpha0(x,y), OD, z/L(3)+0.5);
    Sz = @(x,y,z) s/2 * alpha_z(x,y,z) ./ (alpha_z(x,y,z) + 1);
    B_A_chebfun3 = chebfun3( @(x,y,z) -2/3*kappa_RbXe*gS*(nRb*1e6)*mu0*muB*Sz(x,y,z) * 1e9 .* aperture_mask(x,y), ...  % in unit of nT
                            [-L(1)/2,+L(1)/2, -L(2)/2,+L(2)/2, -L(3)/2,+L(3)/2] ,'eps', 1e-8);     % the 'eps' option here means: absolute error is smaller than <eps>*max(|fun|)
                                                                                                   % 1e-8 might be not enough for very narrow beam
                                                                                                   % higher accuracy will significantly increase calculation time

    P_inc = sum2(chebfun2(@(x,y)I0_fun(x,y) , 1/2*[-L(1),+L(1),-L(2),+L(2)] ));  % in unit of W,  total laser power incident on the vapor cell.
    
    Phi_star = Rrel/sigma_abs;  % in unit of 1/(s*cm^2)
    I_star = Phi_star*h*nu;     % in unit of W/cm^2
    Bm = -s/3*kappa_RbXe*(1e6*nRb)*gS*mu0*muB * 1e9;  % in unit of nT, maximal polarization field

    % pack results
    input_cellPars = cellPars; %#ok<NASGU>
    fields = {'L', 'OD', 'T', 'Rrel', 'nu', 'kappa_RbXe', 'nRb', 'sigma_abs', 's', 'P_inc', 'input_cellPars', 'Phi_star', 'I_star', 'Bm'};
    for ii = 1:length(fields)
        auxInfo.(fields{ii}) = eval(fields{ii});
    end
    
    
    if disp_check
        fprintf('OD: %g, Temperature: %gdegC, n_Rb: %.3E cm^-3\n', OD, T, nRb);
        fprintf('Incident Pump Power: %fmW\n', 1e3*P_inc );


        % 1D plot (alpha_z)
        figure();
        plot(chebfun(@(z)alpha_z(0,0,z)./alpha0(0,0)), [-L(3)/2, +L(3)/2]);
        xlabel('z / cm');
        ylabel('\alpha(0,0,z) / \alpha_0(0,0)');

        % 3D plot
        figure();
        slice(B_A_chebfun3);
        
        % 2D plot (B_A, x=y plane)
        tmp_fun = @(z,x)B_A_chebfun3(x, 0*x, z);
        [x, z] = meshgrid(L(1)*linspace(-1/2, +1/2, 150),  L(3)*linspace(-1/2, 1/2, 100));
        
        funVal = abs(tmp_fun(z,x));
        
        hf = figure('Position', [710 387 582 371]);
        sObj = pcolor(z, x, funVal);   hold on;
        set(sObj,  'FaceColor', 'interp', 'EdgeColor', 'none');
        contour(z, x, funVal, 'LineColor', 0.2*[1 1 1]);
        xlabel('$z$ (cm)', 'interpreter', 'latex');      ylabel('$x$ (cm)', 'interpreter', 'latex');
        box on; grid on;
        title(sprintf('Polarization Field Distribution. Incident Pump Power: %.3fmW', 1e3*P_inc ));
        set(gca,'Layer','top');
        caxis([0,    ceil(abs(2/3*kappa_RbXe*gS*(nRb*1e6)*mu0*muB*0.5 * 1e9)  /  5)*5 ]);
        colorbar(gca,'Position', [0.94 0.10 0.02 0.80]);
        annotation(hf,'textbox',  [0.87 0.96 0.16 0.03], 'String',{'$B_\mathrm{A}$ (nT)'}, 'HorizontalAlignment','center', 'EdgeColor','none','interpreter','latex');
    end

end

%%
% % calc number density
% % use Killian formula, see https://csrc.yuque.com/sensing/basics/plmf5u
% function n = getNumberDensity(temp_degc)
%     kB_val = kB();    % function in OpticalPumping\const
%     T = temp_degc + 273.15;
%     p = 10.^(  9.55 - 4132./T  );   % in unit of Pa
%     n = p./(kB_val*T)/1e6;     % in unit of cm^-3
% end
% 
% % calc number density
% % use Daniel Adam Steck's Rubidium data, see https://csrc.yuque.com/sensing/basics/plmf5u
% function n = getNumberDensity2(temp_degc)
%     kB_val = kB();    % function in OpticalPumping\const
%     T = temp_degc + 273.15;
%     p = 10.^(  9.312 - 4040./T  );   % in unit of Pa
%     n = p./(kB_val*T)/1e6;     % in unit of cm^-3
% end

% calc. the light propagation using alpha_z
% z_L in [0, 1]
function alpha_z = lambertWpropagation(alpha0, od, z_L)
    % assume size(alpha0) == size(tmp)
    tmp = alpha0-od.*z_L;
    alpha_z = nan(size(tmp));
    idx = (tmp > 700);     % for points that tmp>>1, alpha_z will be +inf due to the calculation of exp(tmp).
    alpha_z(~idx) = lambertw(alpha0(~idx).*exp(tmp(~idx)));
    alpha_z(idx) = wrightOmega(tmp(idx) + log(alpha0(idx)));  % use  wrightOmega() instead of lambertw()
end
% use wrightOmega() instead of lambertw() to avoid Inf when alpha0 is very large
% this method runs very slow ......
function alpha_z = lambertWpropagation2(alpha0, od, z_L)
    alpha_z = wrightOmega(alpha0-od.*z_L + log(alpha0));
end
