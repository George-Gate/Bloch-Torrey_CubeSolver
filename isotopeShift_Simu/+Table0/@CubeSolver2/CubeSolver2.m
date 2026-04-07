% Author: GeorgeGate
% Last revision: GeorgeGate@20260208
classdef CubeSolver2 < handle
    %CUBESOLVER2 利用 chebfun 包求解立方体区域中的 Torrey 方程
    % (#Translate# Using the Chebfun package to solve the Torrey equation in a cubic domain)
    
    properties
        cellPars_example;    % an example of cellPars structure
        version = 'v1.1';
    end
    
    methods
        function obj = CubeSolver2()
            % give an example of the cellPars structure
            cellPars = struct('ID', 'CH7');
            cellPars.lambda = [1e-3, 1e-3, 1e-3]; % boundary condition at x, y, z direction
            cellPars.T_lambda = 0.10 * 1.602e-19/kB();   % lambda的特征温度（对应Ebar）   (#Translate# Characteristic temperature of lambda (corresponding to Ebar).)
            cellPars.L = [0.8, 0.6, 0.7];         % cell domain [Lx, Ly, Lz], in unit of cm. 
                                                  % The domain of cell is [-Lx/2, +Lx/2] x [-Ly/2, +Ly/2] x [-Lz/2, +Lz/2]
                                                  % Light propogates from z=-Lz/2 to z=+Lz/2
            cellPars.D = 0.2;               % diffusion constant, cm^2/s
            cellPars.D_alpha = 1;           % 扩散系数的温度依赖，D \propto T^alpha   (#Translate# Temperature dependence of the diffusion coefficient, D \propto T^alpha)
            cellPars.gXe = -2*pi*11.84e-3;  % gyromagnetic ratio of Xe, rad/(s*nT)
            cellPars.G2c = 1/20;            % colisional relaxation rate, 1/s
            cellPars.B0 = 2e4;              % main field, nT
            cellPars.OD = 6;                % optical depth along z direction
            cellPars.Temp = 110;            % cell temperature, degC
            cellPars.kappa_RbXe = 493;      % colisional enhancement factor of Rb-Xe
            cellPars.Rrel = 1/(14e-6);      % relaxation rate of Rb, in unit of 1/s
            cellPars.laserFreq = 377.107e12;  % pump laser frequency, in unit of Hz
            cellPars.pumpPolarization = +1;   % +/-1
            obj.cellPars_example = cellPars;
        end
        
        % Return the Rb polarization field distribution 
        [B_A_chebfun3, auxInfo] = calc_polarization_field(obj, I0_fun, cellPars, aperture_mask, disp_check);
        % Return the coefficient matrix
        [coeff_mat,  eigs_list, b_matrix] = generate_coeff_matrix(obj, B1_chebfun3, basis_3D, cellPars, neigs);
        [coeff_mats, eigs_list, bmats]    = generate_coeff_matrix_with_btot(obj, B1_chebfun3, Bx_chebfun3, By_chebfun3, basis_3D, cellPars, neigs, speed_up_flag)
        
        result = get_pump_power_dependence(obj, I_max_list, pumpBeamProfile, cellPars_129, cellPars_131, print_out, varargin);    % sweep pump beam power
        result = get_cell_temperature_dependence(obj, cellTemp_list, pumpBeamProfile, cellPars_129, cellPars_131, print_out, varargin);  % sweep cell temperature, only nRb and tau change vs. temperature
        result = get_cell_temperature_dependence_v2(obj, cellTemp_list, pumpBeamProfile, cellPars_129, cellPars_131, print_out, varargin);  % sweep cell temperature, both nRb ,tau, D and lambda change vs. temperature

        % benchmark methods
        result = validate_basisInnerProduct(obj);
        result = validate_sweepLinearGradient(obj, saveToFile);
        result = validate_convergencyBasisNum(obj, saveToFile);
        
        function [OD_list, sigma_abs] = get_OD_list(obj, cellTemp_list, cellPars)
            % get absorption cross-section "sigma_abs" from reference OD and cell temperature
            sigma_abs = cellPars.OD/cellPars.L(3)/obj.getRbNumberDensity(cellPars.Temp);    % cm^2
            OD_list = sigma_abs * obj.getRbNumberDensity(cellTemp_list) * cellPars.L(3);
        end
    end
    

    methods(Static)
        
        x = besselzero(n,k,kind);
        % 获取 1D, 3D 基函数 (扩散方程的本征函数)
        % (#Translate# Return 1D and 3D basis functions (eigenfunctions of the diffusion equation).)
        basis_1D = generate_basis_1D(lambda, L, num_basis, disp_check);
        basis_3D = generate_basis_3D(basis_x, basis_y, basis_z, max_basis_num, disp_check);
        % 计算微扰论修正(基模本征值)
        % (#Translate# Calculate perturbative corrections (for the fundamental mode eigenvalue).)
        ptb = calc_perturbation_corrections(b_matrix, cellPars, basis_3D);
%         ptb = calc_perturbation_corrections_old(b_matrix, cellPars, basis_3D);
        
        % 枚举入射光强分布                                  (#Translate# Enumerate the incident light intensity distribution)
        % 光束沿+z方向传播                                  (#Translate# The laser beam propagates along the +z direction)
        % I0_fun(x,y) 为入射光强分布，单位与输入的I_max相同 (#Translate# I0_fun(x,y) is the incident light intensity distribution, with units identical to the input I_max)
        % P_beam为整束光的总功率，单位为[I_max]*cm^2        (#Translate# P_beam is the total power of the entire laser beam, in units of [I_max]*cm^2)
        function [I0_fun, P_beam] = get_input_intensity_distribution(I_max, pumpBeamProfile, beamShape)
            % 枚举光斑类型
            w = pumpBeamProfile.w;                % in unit of cm
            xc = pumpBeamProfile.xc;  
            yc = pumpBeamProfile.yc;
            if strcmpi(beamShape, 'gaussian')
                I0_fun = @(x,y) I_max * exp( - ((x-xc).^2+(y-yc).^2) ./ (2*w*w)  );
                P_beam = I_max *2*pi*w*w;                               %  The total power of the Gaussian beam
            elseif strcmpi(beamShape, 'gaussian_ring')
                r0 = pumpBeamProfile.r0;
                I0_fun = @(x,y) I_max * exp( - ( sqrt((x-xc).^2+(y-yc).^2) - r0    ).^2 ./ (2*w*w)  );   % Input intensity distribution
                P_beam = I_max *2*pi*(w*w*exp(-r0^2/(2*w^2)) + sqrt(pi/2)*r0*w*(  1+erf(r0/(sqrt(2)*w))  )  );  %  The total power of the Gaussian beam
            else
                error('Unknow beamShape "%s".', beamShape);
            end
        end
        
        % 快速计算 fun_chebfun3(x,y,z) 与 ex(x)*ey(y)*ez(z) 的三维重叠积分     (#Translate# Fast algorithm for the 3D overlap integral between fun_chebfun3(x,y,z) and ex(x)*ey(y)*ez(z))
        % ex, ey, ez: chebfun对象，dims: Inf x 1                               (#Translate# ex, ey, ez: chebfun objects, dims: Inf x 1)
        function I = basisInnerProduct(fun_chebfun3, ex, ey, ez)
            cols = ex.' * fun_chebfun3.cols;
            rows = ey.' * fun_chebfun3.rows;
            tubes = ez.' * fun_chebfun3.tubes;
            I = fun_chebfun3.core .* reshape(cols, [], 1, 1).* reshape(rows, 1, [], 1).* reshape(tubes, 1, 1, []);
            I = sum(I(:));
        end
        
        % calc number density
        % use Daniel Adam Steck's Rubidium data, see https://csrc.yuque.com/sensing/basics/plmf5u
        % Assumption: temperature high enough, Rb is in the liquid phase
        function n = getRbNumberDensity(temp_degc)
            kB_val = kB();    
            T = temp_degc + 273.15;
            p = 10.^(  9.318 - 4040./T  );   % in unit of Pa
            n = p./(kB_val*T)/1e6;           % in unit of cm^-3
        end
    end
end


function val = kB()
%KB Summary of this function goes here
%   Detailed explanation goes here
    val=1.3806488e-23;%2010 CODATA: http://physics.nist.gov/cuu/Constants/ unit: J K-1
end
