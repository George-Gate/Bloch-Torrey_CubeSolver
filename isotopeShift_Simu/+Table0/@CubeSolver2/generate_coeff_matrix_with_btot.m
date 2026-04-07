function [coeff_mats, eigs_list, bmats] = generate_coeff_matrix_with_btot(obj, B1_chebfun3, Bx_chebfun3, By_chebfun3, basis_3D, cellPars, neigs, speed_up_flag)
%GENERATE_COEFF_MATRIX_WITH_BTOT 计算横向Torrey方程的系数矩阵 b^(tot)
%    (#Translate# GENERATE_COEFF_MATRIX_WITH_BTOT Computes the coefficient matrix b^(tot) for the transverse Torrey equation.)
% 考虑横向磁场的效应
%    (#Translate# Accounts for the effect of the transverse magnetic field.)
% Bz(x,y,z) = B0 + B1_chebfun3(x,y,z)
% Torrey eq: dM/dt = D*nabla^2(M) - (G2c + 1i*gXe*B0 + 1i*gXe*B1)*M
    N = basis_3D.num_basis;
    if ~exist('neigs', 'var')
        neigs = min(6, N);
    end
    if ~exist('speed_up_flag', 'var')
        speed_up_flag = true;    % used to skip the calculation of bx2_matrix and by2_matrix, which is not necessary. approx. 5x faster if speed_up_flag = true
    end

    L = cellPars.L;
    D = cellPars.D;
    gXe = cellPars.gXe;
    G2c = cellPars.G2c;
    B0 = cellPars.B0;
    
    % convert B1_chebfun3 to chebfun3
    if ~isa(B1_chebfun3, 'chebfun3') && isa(B1_chebfun3, 'function_handle')
        B1_chebfun3 = chebfun3(@(x,y,z)B1_chebfun3(x,y,z), 1/2*[-L(1),+L(1), -L(2),+L(2), -L(3),+L(3)]);
    end
    if ~isa(Bx_chebfun3, 'chebfun3') && isa(Bx_chebfun3, 'function_handle')
        Bx_chebfun3 = chebfun3(@(x,y,z)Bx_chebfun3(x,y,z), 1/2*[-L(1),+L(1), -L(2),+L(2), -L(3),+L(3)]);
    end
    if ~isa(By_chebfun3, 'chebfun3') && isa(By_chebfun3, 'function_handle')
        By_chebfun3 = chebfun3(@(x,y,z)By_chebfun3(x,y,z), 1/2*[-L(1),+L(1), -L(2),+L(2), -L(3),+L(3)]);
    end
    
    if ~speed_up_flag
%         Bx2_chebfun3 = Bx_chebfun3 .* Bx_chebfun3;
%         By2_chebfun3 = By_chebfun3 .* By_chebfun3;
        Bx2_chebfun3 = chebfun3( @(x,y,z) Bx_chebfun3(x,y,z).^2, Bx_chebfun3.domain);  % this is 20% faster than Bx_chebfun3 .* Bx_chebfun3
        By2_chebfun3 = chebfun3( @(x,y,z) By_chebfun3(x,y,z).^2, By_chebfun3.domain);
    end

    coeff_mat = diag(-D*basis_3D.sum_kappa2 - G2c - 1i*gXe*B0);   % Bz = B0 + B1_chebfun3
    b_matrix = nan(size(coeff_mat));    % b_{alpha beta}
    bp_matrix = nan(size(coeff_mat));   % b^+_{alpha beta}
    bm_matrix = nan(size(coeff_mat));   % b^-_{alpha beta}
    if ~speed_up_flag
        bx2_matrix = nan(size(coeff_mat));   % bx^2_{alpha beta}
        by2_matrix = nan(size(coeff_mat));   % by^2_{alpha beta}
    end
%     btot_matrix = nan(size(coeff_mat));   % b^{tot}_{alpha beta}
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
            bx_ii_jj = obj.basisInnerProduct(Bx_chebfun3, ...
                                        basis_3D.basis_x.basis_product{e1(1), e2(1)}, ...    % ex
                                        basis_3D.basis_y.basis_product{e1(2), e2(2)}, ...    % ey
                                        basis_3D.basis_z.basis_product{e1(3), e2(3)});       % ez             
            by_ii_jj = obj.basisInnerProduct(By_chebfun3, ...
                                        basis_3D.basis_x.basis_product{e1(1), e2(1)}, ...    % ex
                                        basis_3D.basis_y.basis_product{e1(2), e2(2)}, ...    % ey
                                        basis_3D.basis_z.basis_product{e1(3), e2(3)});       % ez          
            
            if ~speed_up_flag                        
                bx2_ii_jj = obj.basisInnerProduct(Bx2_chebfun3, ...
                                            basis_3D.basis_x.basis_product{e1(1), e2(1)}, ...    % ex
                                            basis_3D.basis_y.basis_product{e1(2), e2(2)}, ...    % ey
                                            basis_3D.basis_z.basis_product{e1(3), e2(3)});       % ez        
                by2_ii_jj = obj.basisInnerProduct(By2_chebfun3, ...
                                            basis_3D.basis_x.basis_product{e1(1), e2(1)}, ...    % ex
                                            basis_3D.basis_y.basis_product{e1(2), e2(2)}, ...    % ey
                                            basis_3D.basis_z.basis_product{e1(3), e2(3)});       % ez        
            end
                                    
            % update matrix
            b_matrix(ii,jj) = b_ii_jj;          bp_matrix(ii,jj) = bx_ii_jj + 1i*by_ii_jj;  bm_matrix(ii,jj) = bx_ii_jj - 1i*by_ii_jj;   
            b_matrix(jj,ii) = b_matrix(ii,jj);  bp_matrix(jj,ii) = bp_matrix(ii,jj);        bm_matrix(jj,ii) = bm_matrix(ii,jj);         
            if ~speed_up_flag
                bx2_matrix(ii,jj) = bx2_ii_jj ;          by2_matrix(ii,jj) = by2_ii_jj;
                bx2_matrix(jj,ii) = bx2_matrix(ii,jj);   by2_matrix(jj,ii) = by2_matrix(ii,jj);
            end
%             coeff_mat(ii, jj) = coeff_mat(ii, jj) - 1i*gXe*b_ii_jj;
%             if ii ~= jj
%                 coeff_mat(jj, ii) = coeff_mat(jj, ii) - 1i*gXe*b_ii_jj;
%             end
        end
    end
    
    k2 = reshape(basis_3D.sum_kappa2, 1, []);
    btot_matrix = b_matrix + bp_matrix*bm_matrix*gXe./(  2*(gXe*B0 -1i*(D*k2 + G2c)) );
    coeff_mat_b = coeff_mat - 1i*gXe*b_matrix;         % The coefficient matrix calculated using b_{alpha beta}, transverse field effect not counted. 
    coeff_mat_btot = coeff_mat - 1i*gXe*btot_matrix;   % The coefficient matrix calculated using b^{tot}_{alpha beta}, transverse field effect counted.
    
    eigs_list_b = eigs(coeff_mat_b, neigs, 'largestreal');   % eigenvalue list of coeff_mat_b
    eigs_list_b = -sort_eigs(-eigs_list_b, @isAscending);
    eigs_list_btot = eigs(coeff_mat_btot, neigs, 'largestreal');      % eigenvalue list of coeff_mat_btot
    eigs_list_btot = -sort_eigs(-eigs_list_btot, @isAscending);
    
    % pack results
    bmats.b_matrix = b_matrix;                % b_{alpha beta} matrix
    bmats.bp_matrix = bp_matrix;              % b_{alpha beta}^+ matrix
    bmats.bm_matrix = bm_matrix;              % b_{alpha beta}^- matrix
    if ~speed_up_flag
        bmats.bx2_matrix = bx2_matrix;        % int(phi_alpha * Bx^2 * phi_beta) matrix
        bmats.by2_matrix = by2_matrix;
    else
        bmats.bx2_matrix = '(Disabled by speed_up_flag)';
        bmats.by2_matrix = '(Disabled by speed_up_flag)';
    end
    bmats.btot_matrix = btot_matrix;
    coeff_mats.coeff_mat_b = coeff_mat_b;         % coefficient matrix for the transverse Torrey eq.
    coeff_mats.coeff_mat_btot = coeff_mat_btot;  
    eigs_list.eigs_list_b = eigs_list_b;          
    eigs_list.eigs_list_btot = eigs_list_btot;
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
