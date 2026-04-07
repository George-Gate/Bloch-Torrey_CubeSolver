function result = get_cell_temperature_dependence_v2(obj, cellTemp_list, pumpBeamProfile, cellPars_129, cellPars_131, print_out, varargin)
% 扫描气室温度，看本征值的变化
% nRb, tau, D, lambda都会随温度变化
% cellTemp_list: 扫描的格点，气室温度，单位为℃
% pumpBeamProfile: 指定高斯光束的中心、宽度、通光孔径等信息
% cellPars_129, cellPars_131: parameters for Xe129 and Xe131
%
% 假设：Rrel不随温度变化（缓冲气体远多于Rb蒸汽）

% (#Translate#
%     Scan the vapor cell temperature and observe the variation in eigenvalues.
%     Both nRb, tau, D, and lambda will change vs. temperature.
%     cellTemp_list: Grid points for scanning, representing the vapor cell temperature in degrees Celsius (℃).
%     pumpBeamProfile: Specifies information about the Gaussian beam, such as its center, width, and aperture.
%     cellPars_129, cellPars_131: parameters for Xe129 and Xe131
%     
%     Assumption: Rrel does not vary with temperature (the buffer gas greatly outnumbers the Rb vapor).
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

    % 获取温度列表，并计算不同温度下的各种参数  (#Translate# Obtain the temperature list and compute various parameters at different temperatures.)
    cellTemp_list = reshape(cellTemp_list, [], 1);
    Twork = cellPars_129.Temp;   % degC
    Tlist = cellTemp_list + 273.15;   % K
    T0 = Twork + 273.15;   % K
    % get sigma_abs from OD and cellTemp    
    [OD_list, sigma_abs] = obj.get_OD_list(cellTemp_list, cellPars_129);
    % Diffusion constant
    D_list = cellPars_129.D ./ (T0^cellPars_129.D_alpha).*(Tlist.^cellPars_129.D_alpha);
    % Boundary conditions
    lambda129_list = cellPars_129.lambda(1)* exp(  cellPars_129.T_lambda./Tlist -  cellPars_129.T_lambda/T0   );
    lambda131_list = cellPars_131.lambda(1)* exp(  cellPars_131.T_lambda./Tlist -  cellPars_131.T_lambda/T0   );
    % 不同温度下的nRb会在 Table0.CubeSolver2.calc_polarization_field() 中自动根据温度计算
    %   (#Translate# nRb at different temperatures will be automatically computed based on temperature within Table0.CubeSolver2.calc_polarization_field().)
    

    aperture_mask = pumpBeamProfile.aperture_mask;   % 通光孔外形函数，通光区域为1，不通光区域为0      (#Translate# Aperture shape function. Returns 1 in the transmissive region and 0 in the blocking region.)
    I_max = pumpBeamProfile.I_max;
%  example:     aperture_mask = @(x,y,r) 1./((1+((x.*x+y.*y)/r^2).^20));  % 圆形通光孔    (#Translate# Circular aperture)
    
    N = length(cellTemp_list);
    Neigs = 6;
    coeff_mats_129 = {}; coeff_mats_131 = {};
    b_matrix_129 = {};   b_matrix_131 = {};
    eigs_list_129 = nan(N,Neigs)*(1+1i);  eigs_list_131 = nan(N,Neigs)*(1+1i);
    
    % generate laser intensity distribution in the incident plane
    beamShape = p.Results.beamShape;
    [I0_fun, P_beam] = solver.get_input_intensity_distribution(I_max, pumpBeamProfile, beamShape);

    P_inc =  0*cellTemp_list;       %  入射到气室区域的总pump功率  (#Translate# Total pump power incident on the vapor cell region.)
    
    
    numBasis = 100;   % 本征模式数量   (#Translate# The number of 3D basis used in the calculation)
    N1 = ceil(5*( 8*3*numBasis/(4*pi) )^(1/3));     % 预估需要的一维基函数数量  (#Translate# Estimate the required number of 1D basis functions)
    
    tS = tic;
    parfor ii = 1:length(cellTemp_list)
        tic;
        % 更新各种依赖于温度的参数，假设sigma_abs是固定不变的，根据nRb来估算OD的变化
        %   (#Translate# Update various temperature-dependent parameters. It is assumed that sigma_abs remains constant; the change of OD is estimated based on nRb.)
        cP_129 = cellPars_129;            cP_131 = cellPars_131;
        cP_129.Temp = cellTemp_list(ii);  cP_131.Temp = cellTemp_list(ii);  
        cP_129.OD = OD_list(ii);          cP_131.OD = OD_list(ii);
        cP_129.D = D_list(ii);            cP_131.D = D_list(ii);
        cP_129.lambda = lambda129_list(ii)*[1,1,1];
        cP_131.lambda = lambda131_list(ii)*[1,1,1];

        % 因为lambda在随温度变化，所以每个温度下的本征模式需要重新计算
        %   (#Translate# Since lambda varies with temperature, the eigenmodes need to be recalculated for each temperature.)
        basis_x = solver.generate_basis_1D(cP_129.lambda(1), cP_129.L(1), N1, false);     % L(1), L(2), L(3) 不应该相差太远，否则这里的N1可能需要取得更大
        basis_y = solver.generate_basis_1D(cP_129.lambda(2), cP_129.L(2), N1, false);     % (#Translate# L(1), L(2), L(3) should not differ too much, otherwise N1 here may need to be larger)
        basis_z = solver.generate_basis_1D(cP_129.lambda(3), cP_129.L(3), N1, false);   
        basis_3D_129 = solver.generate_basis_3D(basis_x, basis_y, basis_z, numBasis, false);

        basis_x = solver.generate_basis_1D(cP_131.lambda(1), cP_131.L(1), N1, false);      
        basis_y = solver.generate_basis_1D(cP_131.lambda(2), cP_131.L(2), N1, false); 
        basis_z = solver.generate_basis_1D(cP_131.lambda(3), cP_131.L(3), N1, false);
        basis_3D_131 = solver.generate_basis_3D(basis_x, basis_y, basis_z, numBasis, false);
        % 温度太高，极化场太大时，计算B_A_chebfun3可能会报错，所以这里需要一个try结构
        %   (#Translate#  A try-catch block is required here because computing B_A_chebfun3 may fail at very high temperatures or when the polarization field is too strong.)
        try
            % update polarization field distribution
            [B_A_chebfun3, auxInfo] = solver.calc_polarization_field(I0_fun, cP_129, @(x,y)aperture_mask(x,y), false); %#ok<PFBNS>   极化场只依赖于 OD 和 nRb，不需要区分cP_129和cP_131 
            P_inc(ii) = auxInfo.P_inc;                                                                                         %  (#Translate# The polarization field depends only on OD and nRb, so there is no need to distinguish between cP_129 and cP_131.)
            
            % update coefficient matrix
            [coeff_mats_129{ii,1}, eigs_list_129(ii,:) , b_matrix_129{ii,1}] = solver.generate_coeff_matrix(B_A_chebfun3, basis_3D_129, cP_129, Neigs);
            [coeff_mats_131{ii,1}, eigs_list_131(ii,:) , b_matrix_131{ii,1}] = solver.generate_coeff_matrix(B_A_chebfun3, basis_3D_131, cP_131, Neigs);
        catch ME
            b_matrix_129{ii,1} = nan(numBasis);
            b_matrix_131{ii,1} = nan(numBasis);
            fprintf('%d: [Error]%s\n', ii, ME.message);
        end
        % perturbation solution
        ptb_129(ii,1) = solver.calc_perturbation_corrections(b_matrix_129{ii,1}, cP_129, basis_3D_129);
        ptb_131(ii,1) = solver.calc_perturbation_corrections(b_matrix_131{ii,1}, cP_131, basis_3D_131);
        fprintf('%d: %.1fs\n', ii, toc());
    end
    P_inc = mean(P_inc);
    
    % perturbation solution
%     ptb_129 = cellfun(@(bMat)solver.calc_perturbation_corrections(bMat, cellPars_129, basis_3D_129), b_matrix_129);
%     ptb_131 = cellfun(@(bMat)solver.calc_perturbation_corrections(bMat, cellPars_131, basis_3D_131), b_matrix_131);
    ptb_129 = structArray2fieldArray(ptb_129);
    ptb_131 = structArray2fieldArray(ptb_131);
    
    toc(tS);
    %%
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
    result.D_list = D_list;
    result.lambda129_list = lambda129_list;
    result.lambda131_list = lambda131_list;
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
        s.(fn) = reshape([s_array.(fn)], sz);
    end
end