function [basis_3D] = generate_basis_3D(basis_x, basis_y, basis_z, max_basis_num, disp_check)
%GENERATE_BASIS_3D 通过 basis_1D 直积构造出 basis_3D
%  (#Translate# GENERATE_BASIS_3D Constructs the 3D basis via the tensor product of 1D basis functions.)
    chebfun3eps eps;
    if ~exist('max_basis_num', 'var')
        max_basis_num = Inf;
    end
    if ~exist('disp_check', 'var')
        disp_check = true;
    end
    Nx = basis_x.num_basis;  Lx = range(basis_x.domain);
    Ny = basis_y.num_basis;  Ly = range(basis_y.domain);
    Nz = basis_z.num_basis;  Lz = range(basis_z.domain);
    if max_basis_num > Nx*Ny*Nz
        max_basis_num = Nx*Ny*Nz;
    end
    % calc. eigenvalue list
    kappa2_x = (basis_x.kappaL/Lx).^2;    % kappa^2
    kappa2_y = (basis_y.kappaL/Ly).^2;
    kappa2_z = (basis_z.kappaL/Lz).^2;
    cnt = 0;
    eigenvalue_3D = [];  ix = [];  iy = []; iz = [];
    for ii = 1:Nx
        for jj = 1:Ny
            for kk = 1:Nz
                cnt = cnt+1;
                eigenvalue_3D(cnt) = kappa2_x(ii) + kappa2_y(jj) + kappa2_z(kk);    % the eigenvalue of 3D basis
                ix(cnt) = ii;   iy(cnt) = jj; iz(cnt) = kk;
            end
        end
    end
    % sort eigenvalue
    [eigenvalue_3D, sID] = sort(eigenvalue_3D);   
    ix = ix(sID);  iy = iy(sID);  iz = iz(sID);
    
    % construct 3D basis
    bx = basis_x.chebfun;  by = basis_y.chebfun;  bz = basis_z.chebfun;
    new_basis = {};
    for ii = 1:max_basis_num
        new_basis{ii} = chebfun3(@(x,y,z)bx{ix(ii)}(x).*by{iy(ii)}(y).*bz{iz(ii)}(z), [bx{ix(ii)}.domain, by{iy(ii)}.domain, bz{iz(ii)}.domain]);        % chebfun3 object
    end
    
    basis_3D.chebfun = new_basis.';
    basis_3D.sum_kappa2 = eigenvalue_3D(1:max_basis_num).';   % the eigenvalue of basis
    basis_3D.num_basis = max_basis_num;
    basis_3D.lambda = [basis_x.lambda, basis_y.lambda, basis_z.lambda];
    basis_3D.type = 'Laplacian_3D';
    basis_3D.domain = [basis_x.domain, basis_y.domain, basis_z.domain];
    basis_3D.L = [Lx, Ly, Lz];
    basis_3D.mnp = [ix(1:max_basis_num)', iy(1:max_basis_num)', iz(1:max_basis_num)'];
    basis_3D.basis_x = basis_x;
    basis_3D.basis_y = basis_y;
    basis_3D.basis_z = basis_z;
    
    if disp_check
        % eigenvalue check
        calc_kappa2 = cellfun(@(basis) -sum3(basis.*lap(basis)),  new_basis);
        max_rel_kappa2_err = max(abs(calc_kappa2 - eigenvalue_3D(1:max_basis_num))./calc_kappa2)

        % normalization check
        max_normalization_err = max(abs(cellfun(@(basis) sum3(basis.*basis),  new_basis) - 1))

        % function value check
        bID = max_basis_num;  
        npts = 1000;
        tx = Lx*(rand(npts, 1) - 0.5);  ty = Ly*(rand(npts, 1) - 0.5);  tz = Lz*(rand(npts, 1) - 0.5);  
        max_funVal_err = max( abs(  new_basis{bID}(tx, ty, tz)  -  basis_x.fun_hdl{ix(bID)}(tx).*basis_y.fun_hdl{iy(bID)}(ty).*basis_z.fun_hdl{iz(bID)}(tz)     ) )
    end
end

