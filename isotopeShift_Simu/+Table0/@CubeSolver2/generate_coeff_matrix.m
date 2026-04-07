function [coeff_mat, eigs_list, b_matrix] = generate_coeff_matrix(obj, B1_chebfun3, basis_3D, cellPars, neigs)
%GENERATE_COEFF_MATRIX 计算横向Torrey方程的系数矩阵
%   (#Translate# GENERATE_COEFF_MATRIX Computes the coefficient matrix for the transverse Torrey equation.)
% Transverse Torrey eq: dM/dt = D*nabla^2(M) - (G2c + 1i*gXe*B0 + 1i*gXe*B1)*M
    N = basis_3D.num_basis;
    if ~exist('neigs', 'var')
        neigs = min(6, N);
    end

    L = cellPars.L;
    D = cellPars.D;
    gXe = cellPars.gXe;
    G2c = cellPars.G2c;
    B0 = cellPars.B0;
    
    % convert function_handle to @chebfun3
    if ~isa(B1_chebfun3, 'chebfun3') && isa(B1_chebfun3, 'function_handle')
        B1_chebfun3 = chebfun3(@(x,y,z)B1_chebfun3(x,y,z), 1/2*[-L(1),+L(1), -L(2),+L(2), -L(3),+L(3)]);
    end
    
    % calc. overlapping intergral (b matrix)
    coeff_mat = diag(-D*basis_3D.sum_kappa2 - G2c - 1i*gXe*B0);    % Bz = B0 + B1_chebfun3
    b_matrix = nan(size(coeff_mat));
    for ii = 1:N
        for jj = 1:ii
            e1 = basis_3D.mnp(ii,:);
            e2 = basis_3D.mnp(jj,:);
%             b_ii_jj = ( basis_1D.basis_product{e1(1), e2(1)} ).' * B1_chebfun3;    % integrate over x dir  !!! this syntex only valid for B1_chebfun3 with 1x1x1 core !!
%             b_ii_jj = ( basis_1D.basis_product{e1(2), e2(2)} ).' * b_ii_jj;        % integrate over y dir
%             b_ii_jj = ( basis_1D.basis_product{e1(3), e2(3)} ).' * (b_ii_jj.');    % integrate over z dir
            b_ii_jj = obj.basisInnerProduct(B1_chebfun3, ...
                                        basis_3D.basis_x.basis_product{e1(1), e2(1)}, ...    % ex
                                        basis_3D.basis_y.basis_product{e1(2), e2(2)}, ...    % ey
                                        basis_3D.basis_z.basis_product{e1(3), e2(3)});       % ez
            % update matrix
            b_matrix(ii,jj) = b_ii_jj;
            b_matrix(jj,ii) = b_ii_jj;
        end
    end
    coeff_mat = coeff_mat - 1i*gXe*b_matrix;
    
    eigs_list = eigs(coeff_mat, neigs, 'largestreal');
    eigs_list = -sort_eigs(-eigs_list, @isAscending);
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
