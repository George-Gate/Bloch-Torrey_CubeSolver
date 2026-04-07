function result = get_cell_temperature_dependence(obj, cellTemp_list, pumpBeamProfile, cellPars_129, cellPars_131, print_out, varargin)
% 扫描气室温度，看本征值的变化
% cellTemp_list: 扫描的格点，气室温度，单位为℃
% pumpBeamProfile: 指定高斯光束的中心、宽度、通光孔径等信息
% cellPars_129, cellPars_131: parameters for Xe129 and Xe131
%
% 假设：边界条件 lambda 不随温度变化，温度变化只影响极化场的分布
% 假设：Rrel不随温度变化（缓冲气体远多于Rb蒸汽）
%
% (#Translate#
%    Scan the vapor cell temperature and observe the variation in eigenvalues.
%    cellTemp_list: Grid points for scanning, representing the vapor cell temperature in degrees Celsius (℃).
%    pumpBeamProfile: Specifies information about the Gaussian beam, such as its center, width, and aperture.
%    cellPars_129, cellPars_131: parameters for Xe129 and Xe131
%
%    Assumption: The boundary condition lambda does not vary with temperature. Temperature changes only affect the distribution of the polarization field.
%    Assumption: Rrel does not vary with temperature (the buffer gas greatly outnumbers the Rb vapor).
% )


    if ~exist('print_out', 'var')
        print_out = true;
    end
    assert(all(cellTemp_list > 0));

    p = inputParser;
    p.addParameter('beamShape', 'gaussian', @(x)ismember(lower(x), {'gaussian', 'gaussian_ring'}));    %  入射光的空间分布形式  (#Translate# Spatial distribution form of the incident laser beam.)
    p.parse(varargin{:});
    
    % validates the format of cellPars
    % The value of ID, L, B0, OD, Temp, kappa_RbXe, Rrel, laserFreq and pumpPolarization should be the same in both cellPars_129 and cellPars_131
    if ~strcmp(cellPars_129.ID, cellPars_131.ID) || ~all(cellPars_129.L == cellPars_131.L)  ...
          || cellPars_129.B0 ~= cellPars_131.B0  || cellPars_129.OD ~= cellPars_131.OD   ...
          || cellPars_129.Temp ~= cellPars_131.Temp || cellPars_129.kappa_RbXe ~= cellPars_131.kappa_RbXe  ...
          || cellPars_129.Rrel ~= cellPars_131.Rrel || cellPars_129.laserFreq ~= cellPars_131.laserFreq  ...
          || cellPars_129.pumpPolarization ~= cellPars_131.pumpPolarization
        error(sprintf('The following fields in cellPars_129 and cellPars_131 should be the same: \n    ID, L, B0, OD, Temp, kappa_RbXe, Rrel, laserFreq, pumpPolarization'));
    end
    
    solver = obj;

    numBasis = 100;   % 本征模式数量   (#Translate# The number of 3D basis used in the calculation)
    N1 = ceil(5*( 8*3*numBasis/(4*pi) )^(1/3));    % 预估需要的一维基函数数量   (#Translate# Estimate the required number of 1D basis functions)
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

    
    % get sigma_abs from OD and cellTemp    
    cellTemp_list = reshape(cellTemp_list, [], 1);
%     sigma_abs = cellPars_129.OD/cellPars_129.L(3)/obj.getRbNumberDensity(cellPars_129.Temp);    % cm^2
%     OD_list = sigma_abs * obj.getRbNumberDensity(cellTemp_list) * cellPars_129.L(3);
    [OD_list, sigma_abs] = obj.get_OD_list(cellTemp_list, cellPars_129);
    
    aperture_mask = pumpBeamProfile.aperture_mask;   % 通光孔外形函数，通光区域为1，不通光区域为0      (#Translate# Aperture shape function. Returns 1 in the transmissive region and 0 in the blocking region.)
    I_max = pumpBeamProfile.I_max;
%  example:     aperture_mask = @(x,y,r) 1./((1+((x.*x+y.*y)/r^2).^20));  % 圆形通光孔    (#Translate# Circular aperture)
    
    N = length(cellTemp_list);
    Neigs = 6;
    coeff_mats_129 = {}; coeff_mats_131 = {};
    b_matrix_129 = {};   b_matrix_131 = {};
    eigs_list_129 = nan(N,Neigs)*(1+1i);  eigs_list_131 = nan(N,Neigs)*(1+1i);
    polarizationFieldInfo = {};
    
    % generate laser intensity distribution in the incident plane
    beamShape = p.Results.beamShape;
    [I0_fun, P_beam] = solver.get_input_intensity_distribution(I_max, pumpBeamProfile, beamShape);

    P_inc =  0*cellTemp_list;       %  入射到气室区域的总pump功率     (#Translate# Total pump power incident on the vapor cell region.)

    parfor ii = 1:length(cellTemp_list)
        tS = tic;
        try
            cP = cellPars_129;
            cP.Temp = cellTemp_list(ii);   % 根据当前温度更新各种参数   (#Translate# Update various parameters based on the current temperature.)
            cP.OD = OD_list(ii);

            [B_A_chebfun3, auxInfo] = solver.calc_polarization_field(I0_fun, cP, @(x,y)aperture_mask(x,y), false); %#ok<PFBNS>
            P_inc(ii) = auxInfo.P_inc;
            polarizationFieldInfo{ii,1} = auxInfo;
            % 温度只影响极化场，所以下面两句中的 cellPars 和 basis_3D 不需要更新
            % (#Translate# Temperature only affects the polarization field, so cellPars and basis_3D in the following two lines do not need to be updated.)
            [coeff_mats_129{ii,1}, eigs_list_129(ii,:) , b_matrix_129{ii,1}] = solver.generate_coeff_matrix(B_A_chebfun3, basis_3D_129, cellPars_129, Neigs);
            [coeff_mats_131{ii,1}, eigs_list_131(ii,:) , b_matrix_131{ii,1}] = solver.generate_coeff_matrix(B_A_chebfun3, basis_3D_131, cellPars_131, Neigs);
        catch ME
            b_matrix_129{ii,1} = nan(numBasis);
            b_matrix_131{ii,1} = nan(numBasis);
            fprintf('%d: [Error]%s\n', ii, ME.message);
        end
        fprintf('%d: %.1fs\n', ii, toc(tS));
    end
    P_inc = mean(P_inc);
    
    % 微扰解
    ptb_129 = cellfun(@(bMat)solver.calc_perturbation_corrections(bMat, cellPars_129, basis_3D_129), b_matrix_129);
    ptb_131 = cellfun(@(bMat)solver.calc_perturbation_corrections(bMat, cellPars_131, basis_3D_131), b_matrix_131);
    polarizationFieldInfo = cellfun(@(x)x, polarizationFieldInfo);   % convert cell array to struct array
    ptb_129 = structArray2fieldArray(ptb_129);
    ptb_131 = structArray2fieldArray(ptb_131);
    polarizationFieldInfo = structArray2fieldArray(polarizationFieldInfo);
    
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
    result.cellTemp_list = cellTemp_list;
    result.sigma_abs = sigma_abs;
    result.OD_list = OD_list;
    result.pumpBeamProfile = pumpBeamProfile;
    result.pumpBeamProfile.P_inc = P_inc;
    result.pumpBeamProfile.P_beam = P_beam;
    result.T2_129 = -1./real(result.eig0_129);
    result.T2_131 = -1./real(result.eig0_131);
    result.f_129 = abs(imag(result.eig0_129)/2/pi);
    result.f_131 = abs(imag(result.eig0_131)/2/pi);
    result.R = result.f_129./result.f_131;
    result.ptb_129 = ptb_129;
    result.ptb_131 = ptb_131;
    result.ptb_R = abs(imag(result.ptb_129.eig0))./abs(imag(result.ptb_131.eig0));
    result.polarizationFieldInfo = polarizationFieldInfo;

    %%
    if print_out
        hf = figure('Position', [667 377 1048 522]);
        subplot(2,2,1);
        plot(result.cellTemp_list, result.T2_129);  hold on;
        plot(result.cellTemp_list, -1./real(result.ptb_129.eig0), '--'); grid on;
        xlabel('Cell Temperature $\rm (^\circ C)$', 'interpreter', 'latex');
        ylabel('$T_{2,129}$ (s)', 'interpreter', 'latex');
        title(result.type, 'interpreter', 'none');

        subplot(2,2,3);
        plot(result.cellTemp_list, result.T2_131);  hold on;
        plot(result.cellTemp_list, -1./real(result.ptb_131.eig0), '--');  grid on;
        xlabel('Cell Temperature $\rm (^\circ C)$', 'interpreter', 'latex');
        ylabel('$T_{2,131}$ (s)', 'interpreter', 'latex');
        title(sprintf( 'Num. Basis: %d' ,  result.numBasis), 'interpreter', 'latex');   % basis info
        
        subplot(2,2,2);
        plot(result.cellTemp_list, result.f_129); hold on;
        plot(result.cellTemp_list, abs(imag(result.ptb_129.eig0))/2/pi, '--');  grid on;
        xlabel('Cell Temperature $\rm (^\circ C)$', 'interpreter', 'latex');
        ylabel('$f_{129}$ (Hz)', 'interpreter', 'latex');
        title(sprintf( 'Pump Beam: $w/L$=%g, center=$(%g,%g)$cm \n$I_{\\rm max}=%g~{\\rm W/cm^2}$, $P_{\\rm beam}=%g~{\\rm mW}$' ,  result.pumpBeamProfile.w/cellPars_129.L(1), ...
                                                                                                      result.pumpBeamProfile.xc, result.pumpBeamProfile.yc, ...
                                                                                                      result.pumpBeamProfile.I_max, 1e3*result.pumpBeamProfile.P_beam ), 'interpreter', 'latex');   % pump beam info

        subplot(2,2,4);
        plot(result.cellTemp_list, cellPars_129.B0*( result.R - result.R(1)));  hold on;
        plot(result.cellTemp_list, cellPars_129.B0*( result.ptb_R - result.ptb_R(1) ), '--'); grid on;
        xlabel('Cell Temperature $\rm (^\circ C)$', 'interpreter', 'latex');
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
        try
            s.(fn) = reshape([s_array.(fn)], sz);
        end
    end
end