function result = validate_sweepLinearGradient(obj, saveToFile)
% 令B1为线性梯度场 B1 = G*z，对比解析解和数值解的差别，验证CubeSolver2的正确性
% (#Translate# Set B1 as a linear gradient field, B1 = G*z, and compare the difference between the analytical solution and the numerical solution to verify the correctness of @CubeSolver2.)
    
    if ~exist('saveToFile', 'var')
        saveToFile = true;
    end

    % init cellPars
    cellPars = obj.cellPars_example;
    cellPars.lambda = 1e-8*[1,1,1];


    %% 确定线性梯度范围（EP点附近加密网格）
    % (#Translate# Determine the range of the linear gradient (refine the mesh near the EP).)
    
    % 精确的第一EP点位置 PRA 44, 7459 (1991) Eq. (6.20)
    % (#Translate# The exact value of the first exceptional point  [PRA 44, 7459 (1991) Eq. (6.20)])
    g0 = 27*sqrt(3)/32*real(obj.besselzero(-2/3, 1, 1)).^2;
    s0 = g0/sqrt(3);
    G_EP = 8*cellPars.D*g0/abs(cellPars.gXe)/cellPars.L(3)^3;
    Gamma_EP = 4*cellPars.D*s0/cellPars.L(3)^2;

    gList2 = unique([0:0.001:0.02 0:0.2:5, g0(1), g0(1) + linspace(-0.1, + 0.05, 50)]).';
    gList = unique([0:0.2:5, g0(1), g0(1) + linspace(-0.1, + 0.05, 200), gList2']).';
    GList = 8*cellPars.D*gList/abs(cellPars.gXe)/cellPars.L(3)^3;
    GList2 = 8*cellPars.D*gList2/abs(cellPars.gXe)/cellPars.L(3)^3;

    %% calc. the exact eigenvalues using chebop 
    tic;
    eigenvalues_exact = getEigenvalues(gList, 3);
    toc;
    eigenvalues_exact = -4*cellPars.D*eigenvalues_exact/cellPars.L(3)^2;


    %% calc. the eigenvalues using @CubeSolver2 
    
    tic;
    basis_x = obj.generate_basis_1D(cellPars.lambda(1), cellPars.L(1), 50, false); toc;
    basis_y = obj.generate_basis_1D(cellPars.lambda(2), cellPars.L(2), 50, false); toc;
    basis_z = obj.generate_basis_1D(cellPars.lambda(3), cellPars.L(3), 50, false); toc;
    basis_3D = obj.generate_basis_3D(basis_x, basis_y, basis_z, 100, false); toc;

    % sweep linear gradient
    BG = @(G)(@(x,y,z)G*z);

    coeff_mats = {};
    eigs_list = [];
    parfor ii = 1:length(GList2)
        tic;
        [coeff_mats{ii,1}, eigs_list(ii,:) ] = obj.generate_coeff_matrix(BG(GList2(ii)), basis_3D, cellPars, 3);
        fprintf('%d: %.1fs\n', ii, toc());
    end
    
    eigenvalues_solver = eigs_list + cellPars.D* (( basis_x.kappaL(1)/basis_3D.L(1) )^2 + ( basis_y.kappaL(1)/basis_3D.L(2) )^2) + cellPars.G2c + 1i*cellPars.gXe*cellPars.B0;
    
    
    %% calc error
    imag_part = griddedInterpolant(GList, imag(eigenvalues_exact(:,1)));
    real_part = griddedInterpolant(GList, real(eigenvalues_exact(:,1)));
    err_imag = (imag(eigenvalues_solver(:,1)) - imag_part(GList2));
    err_real = (real(eigenvalues_solver(:,1)) - real_part(GList2));
    
    %% pack result
    result = struct('type', 'validate_sweepLinearGradient');
    result.cellPars = cellPars;
    result.basis_3D = basis_3D;
    result.exact_sol.GList = GList;
    result.exact_sol.eigenvalue = eigenvalues_exact;
    result.exact_sol.G_EP = G_EP;
    result.exact_sol.Gamma_EP = Gamma_EP;
    result.solver_sol.GList = GList2;
    result.solver_sol.eigenvalue = eigenvalues_solver;
    result.error.imag_part = err_imag;
    result.error.real_part = err_real;
    
    if saveToFile
        testName = sprintf('validate_sweepLinearGradient_%s', datestr(now, 'yyyymmdd_HHMMSS'));
        save([testName '.mat'], 'result', 'coeff_mats');
    end
    %% plot

    hf = figure('Position', [165 178 1293 673]);
    subplot(2,2,1);
    plot(result.exact_sol.GList,       imag(result.exact_sol.eigenvalue(:,1:2)), '.-'); hold on; grid on;
    plot(result.solver_sol.GList,      imag(result.solver_sol.eigenvalue(:,1:3)), '*');
    plot(G_EP, 0, 'ok');
    xlabel('G (nT/cm)');
    ylabel('\omega correction (rad/s)');
    title(sprintf('number of basis: %d', basis_3D.num_basis));
    subplot(2,2,2);
    plot(result.exact_sol.GList,  -real(result.exact_sol.eigenvalue(:,1:2)), '.-'); hold on; grid on;
    plot(result.solver_sol.GList, -real(result.solver_sol.eigenvalue(:,1:3)), '*');
    plot(G_EP, Gamma_EP, 'ok');
    xlabel('G (nT/cm)');
    ylabel('\Gamma correction (1/s)');
    title(sprintf('\\lambda: [%g, %g, %g]', cellPars.lambda));

    subplot(2,2,3);
    semilogy(result.solver_sol.GList, abs(result.error.imag_part), '*'); hold on; grid on;
    xlabel('G (nT/cm)');
    ylabel('err. of \omega_0 (rad/s)');
    subplot(2,2,4);
    semilogy(result.solver_sol.GList, abs(result.error.real_part), '*'); hold on; grid on;
    xlabel('G (nT/cm)');
    ylabel('err. of \Gamma_0 (1/s)');
    
    if saveToFile
        print(hf, [testName '.png'], '-dpng', '-r330');
        savefig(hf, [testName '.fig']);
    end
end



%% solvers
% get the eigenvalue using chebop() and eigs()
function eigenvalues_exact = getEigenvalues(gList, neigs)
    sol = [];
    parfor ii = 1:length(gList)
        sol(ii) = getSol(gList(ii), neigs);
    end
    
    eigenvalues_exact = arrayfun(@(sol)sort_eigs(diag(sol.D).', @isAscending), sol, 'uniformOutput', false);   % in unit of (q_k)
    eigenvalues_exact = vertcat(eigenvalues_exact{:});
end

% g: Gradient
% neigs: number of eigenvalues to calculate
function sol = getSol(g, neigs)
    L = chebop(@(x,u) -diff(u,2) +1j*g*x*u, [-1,1], 'neumann');
    [sol.EV, sol.D] = eigs(L,neigs);
end

function v = sort_eigs(v, cmpFun)
    sz = size(v);
    v = v(:);
    N = numel(v);
    for ii = N:-1:2
        for kk = 1:ii-1
            if ~cmpFun( v(kk), v(ii) )
                tmp = v(kk);
                v(kk) = v(ii);
                v(ii) = tmp;
            end
        end
    end
    v = reshape(v, sz);
end

function flag = isAscending(a, b)
    th = 1e-6;
    if abs( real(a) - real(b) ) > th && real(a)<real(b)
        flag = true;
    elseif abs( real(a) - real(b) ) < th
        flag = (imag(a) <= imag(b));
    else
        flag = false;
    end
end
    
    


