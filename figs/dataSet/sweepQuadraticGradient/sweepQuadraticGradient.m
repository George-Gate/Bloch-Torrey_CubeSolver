solver = Table0.CubeSolver2();

g129 = -2*pi*11.86015e-3;  % rad/s/nT
g131 = +2*pi*3.515769e-3;  % rad/s/nT

% define cell parameters
cellPars = solver.cellPars_example;    % an example of cellPars structure is defined in the constructor of @CubeSolver2
cellPars.L = [0.8, 0.8, 0.8];
cellPars.D = 0.2;

cellPars_129 = cellPars;
cellPars_129.lambda = 1e-3*[1,1,1];
cellPars_129.gXe = g129;
cellPars_129.G2c = 1/20;      % 1/s

cellPars_131 = cellPars;
cellPars_131.lambda = 2e-2*[1,1,1];
cellPars_131.gXe = g131;
cellPars_131.G2c = 1/20;      % 1/s


fileName = sprintf('sweepQuadraticGradient %s', datestr(now(),'yyyymmdd_HHMMSS'));


% solver = obj;

numBasis = 200;
N1 = ceil(5*( 8*3*numBasis/(4*pi) )^(1/3));
tS = tic;
tic;
basis_x = solver.generate_basis_1D(cellPars_129.lambda(1), cellPars_129.L(1), N1, false);     % L(1), L(2), L(3) 不应该相差太远，否则这里的N1就不能这么取
basis_y = solver.generate_basis_1D(cellPars_129.lambda(2), cellPars_129.L(2), N1, false);
basis_z = solver.generate_basis_1D(cellPars_129.lambda(3), cellPars_129.L(3), N1, false); 
basis_3D_129 = solver.generate_basis_3D(basis_x, basis_y, basis_z, numBasis, false); toc;

basis_x = solver.generate_basis_1D(cellPars_131.lambda(1), cellPars_131.L(1), N1, false);      % L(1), L(2), L(3) 不应该相差太远，否则这里的N1就不能这么取
basis_y = solver.generate_basis_1D(cellPars_131.lambda(2), cellPars_131.L(2), N1, false); 
basis_z = solver.generate_basis_1D(cellPars_131.lambda(3), cellPars_131.L(3), N1, false);
basis_3D_131 = solver.generate_basis_3D(basis_x, basis_y, basis_z, numBasis, false); toc;


B1_G2 = @(x,y,z,G2)0*x  +  0*y  +  G2.*z.*z;
G2 = linspace(-200, +200, 50);    % nT/cm^2

N = length(G2);
Neigs = 6;
coeff_mats_129 = {}; coeff_mats_131 = {};
b_matrix_129 = {};   b_matrix_131 = {};
eigs_list_129 = nan(N,Neigs)*(1+1i);  eigs_list_131 = nan(N,Neigs)*(1+1i);

parfor ii = 1:N
    tic;

    try
        B_A_chebfun3 = chebfun3(@(x,y,z)B1_G2(x,y,z,G2(ii)), basis_3D_129.domain);

        [coeff_mats_129{ii,1}, eigs_list_129(ii,:) , b_matrix_129{ii,1}] = solver.generate_coeff_matrix(B_A_chebfun3, basis_3D_129, cellPars_129, Neigs);
        [coeff_mats_131{ii,1}, eigs_list_131(ii,:) , b_matrix_131{ii,1}] = solver.generate_coeff_matrix(B_A_chebfun3, basis_3D_131, cellPars_131, Neigs);
    catch ME
        b_matrix_129{ii,1} = nan(numBasis);
        b_matrix_131{ii,1} = nan(numBasis);
        fprintf('%d: [Error]%s\n', ii, ME.message);
    end
    fprintf('%d: %.1fs\n', ii, toc());
end

% 微扰解
ptb_129 = cellfun(@(bMat)solver.calc_perturbation_corrections(bMat, cellPars_129, basis_3D_129), b_matrix_129);
ptb_131 = cellfun(@(bMat)solver.calc_perturbation_corrections(bMat, cellPars_131, basis_3D_131), b_matrix_131);
ptb_129 = structArray2fieldArray(ptb_129);
ptb_131 = structArray2fieldArray(ptb_131);

%%
result = struct();
result.solverVersion = solver.version;
result.type = 'sweep_quadratic_gradient';
result.numBasis = numBasis;
result.eig0_129 = eigs_list_129(:,1);
result.eig0_131 = eigs_list_131(:,1);
result.G2_list = G2;
result.T2_129 = -1./real(result.eig0_129);
result.T2_131 = -1./real(result.eig0_131);
result.f_129 = abs(imag(result.eig0_129)/2/pi);
result.f_131 = abs(imag(result.eig0_131)/2/pi);
result.R = result.f_129./result.f_131;
result.ptb_129 = ptb_129;
result.ptb_131 = ptb_131;
result.ptb_R = abs(imag(result.ptb_129.eig0))./abs(imag(result.ptb_131.eig0));


save([fileName '.mat'], 'result', 'cellPars_129', 'cellPars_131', 'fileName');


%%
hf = figure('Position', [2258 84 1021 868]);
ax1 = subplot(2,2,1);
ax2 = subplot(2,2,3);
ax3 = subplot(2,2,2);
ax4 = subplot(2,2,4);

g129 = cellPars_129.gXe;
g131 = cellPars_131.gXe;
R0 = g129/g131;
re = result;
plot(ax1, re.G2_list, re.T2_129);  hold(ax1, 'on'); grid(ax1, 'on');
plot(ax1, re.G2_list, -1./real(re.ptb_129.eig0), '--'); 
plot(ax2, re.G2_list, re.T2_131);  hold(ax2, 'on'); grid(ax2, 'on');
plot(ax2, re.G2_list, -1./real(re.ptb_131.eig0), '--');  
plot(ax3, re.G2_list, re.f_129); hold(ax3, 'on'); grid(ax3, 'on');
plot(ax3, re.G2_list, abs(imag(re.ptb_129.eig0))/2/pi, '--');  
plot(ax4, re.G2_list, cellPars_129.B0*( abs(re.R) - abs(R0)), 'b-','DisplayName', 'Numerical Diagonalization');  hold(ax4, 'on'); grid(ax4, 'on');
plot(ax4, re.G2_list, cellPars_129.B0*( abs(re.ptb_R) - abs(R0) ), 'r--', 'DisplayName', 'Perturbation Approximation');
plot(ax4, re.G2_list, 1i*(R0*re.ptb_131.order1st - re.ptb_129.order1st)/g131, 'm-.', 'DisplayName', '1st order contribution'); 
plot(ax4, re.G2_list, 1i*(R0*re.ptb_131.order3rd_a - re.ptb_129.order3rd_a)/g131+1i*(R0*re.ptb_131.order3rd_b - re.ptb_129.order3rd_b)/g131, 'k-.', 'DisplayName', '3rd order contribution'); 
% plot(ax4, re.G2_list, -1i*(R0*re.ptb_131.order3rd_b - re.ptb_129.order3rd_b)/g131, '--'); 


xlabel(ax1, 'Quadratic Gradient $G_2$ $\rm (nT/cm^2)$', 'interpreter', 'latex');
ylabel(ax1, '$T_{2,129}$ (s)', 'interpreter', 'latex');
title(ax1, re.type, 'interpreter', 'none');

xlabel(ax2, 'Quadratic Gradient $G_2$ $\rm (nT/cm^2)$', 'interpreter', 'latex');
ylabel(ax2, '$T_{2,131}$ (s)', 'interpreter', 'latex');

xlabel(ax3, 'Quadratic Gradient $G_2$ $\rm (nT/cm^2)$', 'interpreter', 'latex');
ylabel(ax3, '$f_{129}$ (Hz)', 'interpreter', 'latex');
% title(ax3, sprintf( 'Pump Beam: $w/L$=%g, center=$(%g,%g)$cm' ,  re.pumpBeamProfile.w/cellPars_129.L(1), re.pumpBeamProfile.xc, re.pumpBeamProfile.yc ), 'interpreter', 'latex');   % pump beam info

xlabel(ax4, 'Quadratic Gradient $G_2$ $\rm (nT/cm^2)$', 'interpreter', 'latex');
ylabel(ax4, '$B_0 \cdot \delta R$ (nT)', 'interpreter', 'latex');

lg = legend(ax4, {}, 'location', 'best');

print(hf, [fileName '.png'], '-dpng', '-r330');

%%
% 将结构体数组转化成单一结构体s，s的域为相应数组
function s = structArray2fieldArray(s_array)
    fnames = fields(s_array(1));
    sz = size(s_array);
    s = struct();
    for ii = 1:length(fnames)
        fn = fnames{ii};
        s.(fn) = reshape([s_array.(fn)], sz);
    end
end