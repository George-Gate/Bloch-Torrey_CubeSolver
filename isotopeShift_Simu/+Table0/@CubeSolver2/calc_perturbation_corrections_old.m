function ptb = calc_perturbation_corrections_old(b_matrix, cellPars, basis_3D)
%calc_perturbation_corrections 计算微扰论对基模本征值的修正
% 注意，这里算出来的修正与文章中的 s^(k)_0 相差一个负号，例如，order1st = -s^(1)_0
    k2  = basis_3D.sum_kappa2;
    D   = cellPars.D;
    G2c = cellPars.G2c;
    gXe = cellPars.gXe;
    B0  = cellPars.B0;
    ptb = struct();
    ptb.order0th = -D*k2(1) - G2c - 1i*gXe*B0;
    ptb.order1st = -1i*gXe*b_matrix(1,1);
    ptb.order2nd = -gXe*gXe/D * sum( b_matrix(1,2:end).^2 ./ reshape(k2(2:end) - k2(1),1,[]) );
    ptb.order3rd_a = +1i*gXe*gXe*gXe/D/D * ( b_matrix(1,2:end)./ reshape(k2(2:end) - k2(1),1,[])   )*b_matrix(2:end,2:end)*(   b_matrix(2:end,1)./ reshape(k2(2:end) - k2(1),[],1)  );
    ptb.order3rd_b = -1i*gXe*gXe*gXe/D/D * b_matrix(1,1) * sum(  b_matrix(1,2:end).^2 ./ reshape(k2(2:end) - k2(1),1,[]).^2  );
    
    ptb.eig0 = ptb.order0th + ptb.order1st + ptb.order2nd + ptb.order3rd_a + ptb.order3rd_b;

end

