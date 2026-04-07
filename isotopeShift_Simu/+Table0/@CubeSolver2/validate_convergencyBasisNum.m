function result = validate_convergencyBasisNum(obj, saveToFile)
% 变化计算中用到的三维基函数的数量，看结果的收敛性
% 用于验证算法的收敛性，并演示不同基函数数量可以得到的最终计算精度
%
% (#Translate#
%    Vary the number of 3D basis used in the calculation to observe the convergence of the results.
%    This method validates the convergence of the algorithm and demonstrates the final computation accuracy achievable with different numbers of basis functions.
% )
    
    if ~exist('saveToFile', 'var')
        saveToFile = true;
    end

    solver = obj;
    cellPars = solver.cellPars_example;
    
    % define "number of basis" list
    numBasis = unique([round(linspace(30, 500, 20)), 1000]);
    numBasis = numBasis(end:-1:1);
    
    tS = tic;
    % define the distribution of B1
    I_max = 0.5;            % W/cm^2
    xc = 0.05;  yc = 0.03;  % cm
    w = 0.25*cellPars.L(1);                % cm
    circleMask = @(x,y,r) 1./((1+((x.*x+y.*y)/r^2).^20));                % 圆形通光孔    (#Translate# Circular aperture)
    I0_fun = @(x,y) I_max * exp( - ((x-xc).^2+(y-yc).^2) ./ (2*w*w)  );  % Gaussian beam
    [B_A_chebfun3, auxInfo] = solver.calc_polarization_field(I0_fun, cellPars, @(x,y)circleMask(x,y, cellPars.L(1)/3), false);

    coeff_mats = {};
    eigs_list = [];
    parfor ii = 1:length(numBasis)
        N = numBasis(ii);
        N1 = ceil(5*( 8*3*N/(4*pi) )^(1/3));
        tic;
        basis_x = solver.generate_basis_1D(cellPars.lambda(1), cellPars.L(1), N1, false); 
        basis_y = solver.generate_basis_1D(cellPars.lambda(2), cellPars.L(2), N1, false); 
        basis_z = solver.generate_basis_1D(cellPars.lambda(3), cellPars.L(3), N1, false); 
        basis_3D = solver.generate_basis_3D(basis_x, basis_y, basis_z, N, false); 

        [coeff_mats{ii,1}, eigs_list(ii,:), b_matrix ] = solver.generate_coeff_matrix(B_A_chebfun3, basis_3D, cellPars, 15);
        fprintf('%d: %.1fs\n', ii, toc());
    end
    toc(tS);
    standard_eigs = eigs_list(1,:);
    rel_err = (eigs_list(2:end,:) - standard_eigs)./abs(standard_eigs);
    
    %% pack result
    
    result = struct('type', 'validate_convergencyBasisNum');
    result.cellPars = cellPars;
    result.B1.auxInfo = auxInfo;
    result.B1.chebfun3 = B_A_chebfun3;
    result.B1.w = w;
    result.B1.xc = xc;
    result.B1.yc = yc;
    result.B1.I_max = I_max;
    result.numBasis_list = numBasis;
    result.eigs_list = eigs_list;
    result.relative_error = rel_err;
    result.standard_eigs = standard_eigs;
    
    if saveToFile
        testName = sprintf('validate_convergencyBasisNum_%s', datestr(now, 'yyyymmdd_HHMMSS'));
        save([testName '.mat'], 'result', 'coeff_mats');
    end

%%

    hf = figure('Position', [724 141 678 849]);
    semilogy(numBasis(2:end), abs(rel_err));
    grid on;
    xlabel('Number of 3D basis');
    ylabel('|relative error of eigenvalues|');
    title({['validate\_convergencyBasisNum(), standard: ', num2str(numBasis(1)), ' 3D bases'], 'B_1: polarization field of Gaussian beam'});
    
    if saveToFile
        print(hf, [testName '.png'], '-dpng', '-r330');
        savefig(hf, [testName '.fig']);
    end
    
end
