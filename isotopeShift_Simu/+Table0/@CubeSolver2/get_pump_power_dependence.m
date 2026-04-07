function result = get_pump_power_dependence(obj, I_max_list, pumpBeamProfile, cellPars_129, cellPars_131, print_out, varargin)
% 扫描高斯光束的功率，看本征值的变化
% I_max_list: 扫描的格点，I_max为高斯光束中心的入射光强，单位 W/cm^2
% pumpBeamProfile: 指定高斯光束的中心、宽度、通光孔径等信息
% cellPars_129, cellPars_131: Xe129, Xe131的气室参数
%
% (#Translate#
%     Scan the power of the Gaussian beam and observe the variation in eigenvalues.
%     I_max_list: Grid points for scanning. I_max is the incident intensity at the center of the Gaussian beam, in units of W/cm^2.
%     pumpBeamProfile: Specifies information about the Gaussian beam, such as its center, width, and aperture.
%     cellPars_129, cellPars_131: Vapor cell parameters for Xe129 and Xe131.
% )

    if ~exist('print_out', 'var')
        print_out = true;
    end
    
    flag_use_btot = false;   % consider the effect of transverse magnetic field, use {b_tot} instead of {b} matrix
    
    % Check the validity of cellPars
    % The value of ID, L, B0, OD, Temp, kappa_RbXe, Rrel, laserFreq and pumpPolarization should be the same in both cellPars_129 and cellPars_131
    if ~strcmp(cellPars_129.ID, cellPars_131.ID) || ~all(cellPars_129.L == cellPars_131.L)  ...
          || cellPars_129.B0 ~= cellPars_131.B0  || cellPars_129.OD ~= cellPars_131.OD   ...
          || cellPars_129.Temp ~= cellPars_131.Temp || cellPars_129.kappa_RbXe ~= cellPars_131.kappa_RbXe  ...
          || cellPars_129.Rrel ~= cellPars_131.Rrel || cellPars_129.laserFreq ~= cellPars_131.laserFreq  ...
          || cellPars_129.pumpPolarization ~= cellPars_131.pumpPolarization
        error(sprintf('The following fields in cellPars_129 and cellPars_131 should be the same: \n    ID, L, B0, OD, Temp, kappa_RbXe, Rrel, laserFreq, pumpPolarization'));
    end
    
    p = inputParser;
    p.addParameter('beamShape', 'gaussian', @(x)ismember(lower(x), {'gaussian', 'gaussian_ring'}));    %  入射光的空间分布形式  (#Translate# Spatial distribution form of the incident laser beam.)
    p.parse(varargin{:});
    
    solver = obj;

    numBasis = 100;   % 本征模式数量   (#Translate# The number of 3D basis used in the calculation)
    N1 = ceil(5*( 8*3*numBasis/(4*pi) )^(1/3));  % 预估需要的一维基函数数量   (#Translate# Estimate the required number of 1D basis functions)
    tS = tic;
    tic;
    basis_x = solver.generate_basis_1D(cellPars_129.lambda(1), cellPars_129.L(1), N1, false);     % L(1), L(2), L(3) 不应该相差太远，否则这里的N1可能需要取得更大
    basis_y = solver.generate_basis_1D(cellPars_129.lambda(2), cellPars_129.L(2), N1, false);     % (#Translate# L(1), L(2), L(3) should not differ too much, otherwise N1 here may need to be larger)
    basis_z = solver.generate_basis_1D(cellPars_129.lambda(3), cellPars_129.L(3), N1, false); 
    basis_3D_129 = solver.generate_basis_3D(basis_x, basis_y, basis_z, numBasis, false); toc;

    basis_x = solver.generate_basis_1D(cellPars_131.lambda(1), cellPars_131.L(1), N1, false);    
    basis_y = solver.generate_basis_1D(cellPars_131.lambda(2), cellPars_131.L(2), N1, false); 
    basis_z = solver.generate_basis_1D(cellPars_131.lambda(3), cellPars_131.L(3), N1, false);
    basis_3D_131 = solver.generate_basis_3D(basis_x, basis_y, basis_z, numBasis, false); toc;

    
    I_max_list = reshape(I_max_list, [], 1);
    
    aperture_mask = pumpBeamProfile.aperture_mask;   % 通光孔外形函数，通光区域为1，不通光区域为0      (#Translate# Aperture shape function. Returns 1 in the transmissive region and 0 in the blocking region.)
%  example:     aperture_mask = @(x,y,r) 1./((1+((x.*x+y.*y)/r^2).^20));  % 圆形通光孔    (#Translate# Circular aperture)
    
    N = length(I_max_list);
    Neigs = 6;
    coeff_mats_129 = {}; coeff_mats_131 = {};
    b_matrix_129 = {};   b_matrix_131 = {};
    eigs_list_129 = nan(N,Neigs)*(1+1i);  eigs_list_131 = nan(N,Neigs)*(1+1i);
    P_inc_list =  0*I_max_list;       %  入射到气室区域的总pump功率     (#Translate# Total pump power incident on the vapor cell region.)
    P_beam_list = 0*I_max_list;       %  整束高斯光束的总功率           (#Translate# Total power of the entire Gaussian laser beam.)
    
    beamShape = p.Results.beamShape;
    
    parfor ii = 1:length(I_max_list)
        tic;
        % generate laser intensity distribution in the incident plane
        [I0_fun, P_beam_list(ii)] = solver.get_input_intensity_distribution(I_max_list(ii), pumpBeamProfile, beamShape);

        warning off;
        try
            [B_A_chebfun3, auxInfo{ii}] = solver.calc_polarization_field(I0_fun, cellPars_129, @(x,y)aperture_mask(x,y), false); %#ok<PFBNS>
            P_inc_list(ii) = auxInfo{ii}.P_inc;
            if flag_use_btot
                % used to verify if the effect of transverse field is negligible
                [coeff_mat_tot_129, eigs_list_tot_129, b_matrix_tot_129] = solver.generate_coeff_matrix_with_btot(B_A_chebfun3, B_A_chebfun3/100, B_A_chebfun3/500, basis_3D_129, cellPars_129, Neigs, true);
                [coeff_mat_tot_131, eigs_list_tot_131, b_matrix_tot_131] = solver.generate_coeff_matrix_with_btot(B_A_chebfun3, B_A_chebfun3/100, B_A_chebfun3/500, basis_3D_131, cellPars_131, Neigs, true);
                coeff_mats_129{ii,1} = coeff_mat_tot_129.coeff_mat_btot;  eigs_list_129(ii,:) = eigs_list_tot_129.eigs_list_btot;    b_matrix_129{ii,1} = b_matrix_tot_129.btot_matrix;
                coeff_mats_131{ii,1} = coeff_mat_tot_131.coeff_mat_btot;  eigs_list_131(ii,:) = eigs_list_tot_131.eigs_list_btot;    b_matrix_131{ii,1} = b_matrix_tot_131.btot_matrix;
            else
                [coeff_mats_129{ii,1}, eigs_list_129(ii,:) , b_matrix_129{ii,1}] = solver.generate_coeff_matrix(B_A_chebfun3, basis_3D_129, cellPars_129, Neigs);
                [coeff_mats_131{ii,1}, eigs_list_131(ii,:) , b_matrix_131{ii,1}] = solver.generate_coeff_matrix(B_A_chebfun3, basis_3D_131, cellPars_131, Neigs);
            end
        catch ME
            b_matrix_129{ii,1} = nan(numBasis);
            b_matrix_131{ii,1} = nan(numBasis);
            fprintf('%d: [Error]%s\n', ii, ME.message);
        end
        warning on;
        fprintf('%d: %.1fs\n', ii, toc());
    end

    % characteristic values
    cinfo = struct;
    cinfo.Phi_star = auxInfo{1}.Phi_star;
    cinfo.I_star = auxInfo{1}.I_star;
    cinfo.nRb = auxInfo{1}.nRb;
    cinfo.sigma_abs = auxInfo{1}.sigma_abs;
    cinfo.Bm = auxInfo{1}.Bm;
    
    % perturbation solution
    ptb_129 = cellfun(@(bMat)solver.calc_perturbation_corrections(bMat, cellPars_129, basis_3D_129), b_matrix_129);
    ptb_131 = cellfun(@(bMat)solver.calc_perturbation_corrections(bMat, cellPars_131, basis_3D_131), b_matrix_131);
    ptb_129 = structArray2fieldArray(ptb_129);
    ptb_131 = structArray2fieldArray(ptb_131);
    
    toc(tS);
    %% pack results
    result = struct();
    result.solverVersion = solver.version;
    result.type = 'pump_power_dependence';
    result.numBasis = numBasis;
    result.eig0_129 = eigs_list_129(:,1);
    result.eig0_131 = eigs_list_131(:,1);
    result.eigs_129 = eigs_list_129;
    result.eigs_131 = eigs_list_131;
    result.I_max_list = I_max_list;
    result.pumpBeamProfile = pumpBeamProfile;
    result.P_inc_list = P_inc_list;
    result.P_beam_list = P_beam_list;
    result.T2_129 = -1./real(result.eig0_129);
    result.T2_131 = -1./real(result.eig0_131);
    result.f_129 = abs(imag(result.eig0_129)/2/pi);
    result.f_131 = abs(imag(result.eig0_131)/2/pi);
    result.R = result.f_129./result.f_131;
    result.ptb_129 = ptb_129;
    result.ptb_131 = ptb_131;
    result.ptb_R = abs(imag(result.ptb_129.eig0))./abs(imag(result.ptb_131.eig0));
    result.characteristic_info = cinfo;
    
%     result.coeff_mats_129 = coeff_mats_129;  
%     result.b_matrix_129 = b_matrix_129;  
%     result.coeff_mats_131 = coeff_mats_131;
%     result.b_matrix_131 = b_matrix_131;

    %%
    if print_out
        hf = figure('Position', [667 377 1048 522]);
        subplot(2,2,1);
        plot(result.I_max_list, result.T2_129);  hold on;
        plot(result.I_max_list, -1./real(result.ptb_129.eig0), '--'); grid on;
        xlabel('Center intensity $I_0$ $\rm (W/cm^2)$', 'interpreter', 'latex');
        ylabel('$T_{2,129}$ (s)', 'interpreter', 'latex');
        title(result.type, 'interpreter', 'none');

        subplot(2,2,3);
        plot(result.I_max_list, result.T2_131);  hold on;
        plot(result.I_max_list, -1./real(result.ptb_131.eig0), '--');  grid on;
        xlabel('Center intensity $I_0$ $\rm (W/cm^2)$', 'interpreter', 'latex');
        ylabel('$T_{2,131}$ (s)', 'interpreter', 'latex');

        subplot(2,2,2);
        plot(result.I_max_list, result.f_129); hold on;
        plot(result.I_max_list, abs(imag(result.ptb_129.eig0))/2/pi, '--');  grid on;
        xlabel('Center intensity $I_0$ $\rm (W/cm^2)$', 'interpreter', 'latex');
        ylabel('$f_{129}$ (Hz)', 'interpreter', 'latex');
        title(sprintf( 'Pump Beam: $w/L$=%g, center=$(%g,%g)$cm' ,  result.pumpBeamProfile.w/cellPars_129.L(1), result.pumpBeamProfile.xc, result.pumpBeamProfile.yc ), 'interpreter', 'latex');   % pump beam info

        subplot(2,2,4);
        plot(result.I_max_list, cellPars_129.B0*( result.R - result.R(1)));  hold on;
        plot(result.I_max_list, cellPars_129.B0*( result.ptb_R - result.ptb_R(1) ), '--'); grid on;
        xlabel('Center intensity $I_0$ $\rm (W/cm^2)$', 'interpreter', 'latex');
        ylabel('$B_0 \cdot \delta R$ (nT)', 'interpreter', 'latex');
    end
        

end

% 将结构体数组转化成单一结构体s，s的域为相应数组
% (#Translate# Convert the struct array into a single struct, whose fields correspond to the respective arrays.)
function s = structArray2fieldArray(s_array)
    fnames = fields(s_array(1));
    sz = size(s_array);
    s = struct();
    for ii = 1:length(fnames)
        fn = fnames{ii};
        s.(fn) = reshape([s_array.(fn)], sz);
    end
end