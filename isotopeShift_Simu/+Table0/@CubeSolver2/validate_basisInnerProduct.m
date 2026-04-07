function result = validate_basisInnerProduct(obj)
%  验证 CubeSolver2.basisInnerProduct() 方法计算重叠积分的正确性
%  basisInnerProduct()方法使用了 chebfun3 的内部数据结构来加速计算，未来 Chebfun 包更新后方法可能失效
%  建议每次更新 Chebfun 包后，再验证一次basisInnerProduct()方法。
%
% (#Translate#
%    Validate the correctness of the CubeSolver2.basisInnerProduct() method for computing overlap integrals.
%    The basisInnerProduct() method utilizes the internal data structures of chebfun3 object to accelerate computation. This approach may become invalid after future updates to the Chebfun package.
%    It is recommended to re-validate the basisInnerProduct() method after each Chebfun package update.
% )

    Dx = [-2, 1];
    Dy = [-3, 1/2];
    Dz = [0, 5];

    fx = chebfun(@(x)x.^2, Dx);
    fy = chebfun(@(y)exp(y), Dy);
    fz = chebfun(@(z)z, Dz);

    fun = @(a)chebfun3(@(x,y,z)x.*y.*y./(a+z.*z)+cos(z),  [Dx, Dy, Dz]);

    anaSol = @(a)15*( 5*exp(7/2)-68 )*( log(a)-log(25+a) )/( 32*exp(3) ) + 3*(exp(7/2)-1)*(cos(5)+5*sin(5)-1)/exp(3);

    aList = [linspace(1, 100, 100), 2+1i*linspace(-90, 90, 100)];
    numSol = nan(size(aList));

    tic;
    for ii = 1:length(aList)
        numSol(ii) = obj.basisInnerProduct(fun(aList(ii)), fx, fy, fz);
    end
    toc;
    anaSol = anaSol(aList);
    err = numSol - anaSol;
    
    result = struct('type', 'validate_basisInnerProduct');
    result.aList = aList;
    result.sol_by_basisInnerProduct = numSol;
    result.sol_analytical = anaSol;
    result.error = err;
    
    
    %%
    figure();
    subplot(2,1,1);
    plot(abs(aList), anaSol, '-');  hold on;
    plot(abs(aList), numSol, '.');
    xlabel('|a|');
    ylabel('I');
    title('validation of basisInnerProduct()');

    subplot(2,1,2);
    semilogy(abs(aList), abs(err));
    xlabel('|a|');
    ylabel('abs error');

end