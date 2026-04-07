function basis = generate_basis_1D(lambda, L, num_basis, disp_check)
%GET_BASIS_1D 获取一维的非微扰基函数 (拉普拉斯算子的一维本征函数，带边界条件)
% (#Translate# GET_BASIS_1D Obtains the one-dimensional non-perturbative basis functions (eigenfunctions of the Laplace operator with boundary conditions).)
%  lambda: boundary condition, lambda/L*fun + dfun/dn = 0 on boundary.
%  L: cubic side length, in unit of cm
%  num_basis: number of basis to generate
%  disp_check: wether to perform and display accuracy check
%
    if ~exist('disp_check', 'var')
        disp_check = true;
    end

    kappaL = get_kappaL(lambda, num_basis);
    delta = atan(kappaL./lambda) + kappaL/2;
    
    % 1D unperturbed basis
    % matlab function handle
    basis_handle = arrayfun(@(kL, d)(@(x)1/sqrt(  L/2 + lambda*L/(kL*kL+lambda*lambda)  ) * sin(kL*x/L+d) ), kappaL, delta, 'UniformOutput', false);
    % cell array of chebfun
    basis_1D = arrayfun(@(hdl)chebfun(hdl , [-L/2, +L/2] ), basis_handle, 'UniformOutput', false);
    % basis{i} x basis{j}
    basis_product = {};
    N = length(basis_1D);
    for ii = 1:N
        for jj = 1:ii
            basis_product{ii, jj} = basis_1D{ii} .* basis_1D{jj};
            basis_product{jj, ii} = basis_product{ii, jj};
        end
    end
    
    if disp_check
        % chebfun matrix
        basis_1D_matrix = horzcat(basis_1D{:});
        % check orthonormalization
        basis_int = basis_1D_matrix' * basis_1D_matrix;
        max_normalization_err = max(abs(diag(basis_int) - 1))
        max_orth_err = max(max(abs(basis_int - diag(diag(basis_int)))))

        % check interpolation error
        bID = num_basis;
        xList = linspace(-L/2, +L/2, 200);
        interpolation_err = max(abs(  basis_handle{bID}(xList) - basis_1D{bID}(xList)  ))
    end
    
    % pack results
    basis.fun_hdl = basis_handle;
    basis.chebfun = basis_1D;
    basis.basis_product = basis_product;   % basis_product{i,j} = basis_1D{i} * basis_1D{j}
    basis.domain = [-L/2, +L/2];
    basis.num_basis = num_basis;
    basis.lambda = lambda;
    basis.kappaL = kappaL;
    basis.delta = delta;
    basis.type = 'Laplacian_1D';
end

%%
% use roots() function of chebfun package to find kappaL
% this method is more accurate than the eigs() method
function kappaL = get_kappaL(lambda, neigs)
    domain = [sqrt(lambda), pi*(1/2:neigs) ,neigs*pi+10*eps];    % 手动给出定义域中的发散点 (#Translate# Manually specify divergence points in the domain.)
    targetFun = @(kL)(tan(kL) - 2*lambda.*kL./(kL.*kL-lambda.*lambda));
    fun1 = chebfun( targetFun, domain ,'exps',[0, -ones(1, neigs), 0]);   % 指定发散点的发散指数 (#Translate# Specify the divergence exponent for the divergence points.)
    kappaL = roots(fun1, 'nojump');  % 'nojump' 参数可以将 -inf --> +inf 的跳跃点排除在根的范围外 (#Translate# The 'nojump' parameter excludes jump points from -inf to +inf from the root-finding range.)
%     figure();
%     plot(fun1);  hold on;
%     plot(kappaL,0*kappaL, 'ro');  grid on;
%     ylim([-5 5]);
    kappaL = reshape(kappaL(1:neigs),1,[]);
    if max(abs(targetFun(kappaL))) > 1e-9
        warning('CubeSolver2::generate_basis_1D.get_kappaL(): The solution of kappaL may not accurate.');
    end
end