%% standard use
solver = Table0.CubeSolver2();

% define cell parameters
cellPars = solver.cellPars_example;    % an example of cellPars structure is defined in the constructor of @CubeSolver2
cellPars.L = [0.8, 0.8, 0.6];
cellPars.lambda = [1e-3, 1e-3, 1e-4];


% Although there is a "cellPars_example" property in solver object, the methods of @CubeSolver2 will never directly use this property.
% All methods use the cellPars input via the function call.

% generate basis functions
tic;
basis_x = Table0.CubeSolver2.generate_basis_1D(cellPars.lambda(1), cellPars.L(1), 20, false); toc;
basis_y = Table0.CubeSolver2.generate_basis_1D(cellPars.lambda(2), cellPars.L(2), 50, false); toc;
basis_z = Table0.CubeSolver2.generate_basis_1D(cellPars.lambda(3), cellPars.L(3), 40, false); toc;
basis_3D = Table0.CubeSolver2.generate_basis_3D(basis_x, basis_y, basis_z, 100, false); toc;

% generate polarization field distribution
w = 0.25*cellPars.L(1);
I_max = 0.2;  
I0_fun = @(x,y) I_max * exp( - (x.*x+y.*y) ./ (2*w*w)  ); % in unit of W/cm^2

circleMask = @(x,y,r) 1./((1+((x.*x+y.*y)/r^2).^20));    % 圆形通光孔    (#Translate# Circular aperture)
[B_A_chebfun3, auxInfo] = solver.calc_polarization_field(I0_fun, cellPars, @(x,y)circleMask(x,y,0.3), true); toc;
% [B_A_chebfun3, auxInfo] = solver.calc_polarization_field(I0_fun, cellPars, @(x,y)ones(size(x)), true); toc;

% generate coefficient matrix, and find its eigenvalues
neigs = 6;
[coeff_mat, eigs_list, b_matrix] = solver.generate_coeff_matrix(B_A_chebfun3, basis_3D, cellPars, neigs); toc;
[coeff_mat_tot, eigs_list_tot, b_matrix_tot] = solver.generate_coeff_matrix_with_btot(B_A_chebfun3, B_A_chebfun3/100, B_A_chebfun3/500, basis_3D, cellPars, neigs); toc;

% calc. perturbation corrections
ptb = solver.calc_perturbation_corrections(b_matrix, cellPars, basis_3D);
ptb_with_btot = solver.calc_perturbation_corrections(b_matrix_tot.btot_matrix, cellPars, basis_3D);

%
figure();
bar3(b_matrix);
title('b matrix');
zlabel('matrix element (nT)');


%% sweep pump power
solver = Table0.CubeSolver2();

g129 = -2*pi*11.86015e-3;  % rad/s/nT
g131 = +2*pi*3.515769e-3;  % rad/s/nT

% define cell parameters
cellPars = solver.cellPars_example;    % an example of cellPars structure is defined in the constructor of @CubeSolver2
cellPars.L = [0.8, 0.8, 0.8];

cellPars_129 = cellPars;
cellPars_129.lambda = 1e-3*[1,1,1];
cellPars_129.gXe = g129;
cellPars_129.G2c = 1/20;      % 1/s


cellPars_131 = cellPars;
cellPars_131.lambda = 2e-2*[1,1,1];
cellPars_131.gXe = g131;
cellPars_131.G2c = 1/20;      % 1/s

pumpBeamProfile = struct();
pumpBeamProfile.w = 0.5*cellPars.L(1);   % cm
pumpBeamProfile.xc = 0;
pumpBeamProfile.yc = 0;
circleMask = @(x,y,r) 1./((1+((x.*x+y.*y)/r^2).^20));  % 圆形通光孔    (#Translate# Circular aperture)
pumpBeamProfile.aperture_mask = @(x,y)circleMask(x,y,0.5);     % 通光孔外形函数，通光区域为1，不通光区域为0      (#Translate# Aperture shape function. Returns 1 in the transmissive region and 0 in the blocking region.)

I_max_list = unique([linspace(0, 2, 16), linspace(0, 0.3, 16)]);

result = solver.get_pump_power_dependence(I_max_list, pumpBeamProfile, cellPars_129, cellPars_131, true, 'beamShape', 'Gaussian');
