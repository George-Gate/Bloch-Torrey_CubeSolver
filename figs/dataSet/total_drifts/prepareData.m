clear; clc;

warning off

% bA, dbA/dT, dbA/dP
dataC1 = load(['sweepPumpPower2' filesep 'sweepPumpPower_varyBeamWidth_v3 L=0.30cm 20250626_210140.mat']);   % L=0.3cm, (xc,yc)/L=(0.02,0.03)
dataC2 = load(['sweepPumpPower2' filesep 'sweepPumpPower_varyBeamWidth_v3 L=0.50cm 20250627_023431.mat']);   % L=0.5cm, (xc,yc)/L=(0.02,0.03)
dataC3 = load(['sweepPumpPower2' filesep 'sweepPumpPower_varyBeamWidth_v3 L=0.80cm 20250627_074617.mat']);   % L=0.8cm, (xc,yc)/L=(0.02,0.03)

% dbA/dxc, dbA/dsigma_abs
dataA1 = load(['sweepPumpPower4' filesep 'sweepPumpPower_varyBeamWidth_v4 L=0.30cm 20250628_151847.mat']);   % L=0.3cm, (xc,yc)/L=(0.02,0.03)
dataA2 = load(['sweepPumpPower4' filesep 'sweepPumpPower_varyBeamWidth_v4 L=0.50cm 20250628_215718.mat']);   % L=0.5cm, (xc,yc)/L=(0.02,0.03)
dataA3 = load(['sweepPumpPower4' filesep 'sweepPumpPower_varyBeamWidth_v4 L=0.80cm 20250629_041309.mat']);   % L=0.8cm, (xc,yc)/L=(0.02,0.03)
% dbA/dRrel
dataB1 = load(['sweepPumpPower5' filesep 'sweepPumpPower_varyBeamWidth_v5 L=0.30cm 20250701_023748.mat']);   % L=0.3cm, (xc,yc)/L=(0.02,0.03)
dataB2 = load(['sweepPumpPower5' filesep 'sweepPumpPower_varyBeamWidth_v5 L=0.50cm 20250701_070342.mat']);   % L=0.5cm, (xc,yc)/L=(0.02,0.03)
dataB3 = load(['sweepPumpPower5' filesep 'sweepPumpPower_varyBeamWidth_v5 L=0.80cm 20250701_110630.mat']);   % L=0.8cm, (xc,yc)/L=(0.02,0.03)

warning on;


%% construct data
data1A = constructData1(dataC1);
data2A = constructData1(dataC2);
data3A = constructData1(dataC3);

data1B = constructData2(dataA1, dataB1);
data2B = constructData2(dataA2, dataB2);
data3B = constructData2(dataA3, dataB3);

dataList = {{data3A, data3B};
            {data2A, data2B};
            {data1A, data1B}};

save('preparedData.mat', 'dataList');

%% aux funs
function data = constructData2(dataA, dataB)
    P_beam_mat = {};
    w_L_mat = {};
    R_mat = {};
    bA_mat_xc = {};
    bA_mat_sigma = {};
    bA_mat_Rrel = {};

    g129 = dataA.cellPars_129.gXe;
    g131 = dataA.cellPars_131.gXe;
    R0 = g129/g131;

    for ii = 1:2
        tmp = cellfun(@(re)re.R, dataA.result_xc(:,ii), 'UniformOutput', false);
        R_mat{ii} = cell2mat(tmp');
        bA_mat_xc{ii} = -dataA.cellPars_129.B0*( abs(R_mat{ii}) - abs(R0))/(R0);
    
        tmp = cellfun(@(re)re.R, dataA.result_sigma_abs(:,ii), 'UniformOutput', false);
        R_mat{ii} = cell2mat(tmp');
        bA_mat_sigma{ii} = -dataA.cellPars_129.B0*( abs(R_mat{ii}) - abs(R0))/(R0);
    
        tmp = cellfun(@(re)re.R, dataB.result_Rrel(:,ii), 'UniformOutput', false);
        R_mat{ii} = cell2mat(tmp');
        bA_mat_Rrel{ii} = -dataB.cellPars_129.B0*( abs(R_mat{ii}) - abs(R0))/(R0);

        tmp = cellfun(@(re)re.P_beam_list, dataA.result_xc(:,ii), 'UniformOutput', false);
        P_beam_mat{ii} = cell2mat(tmp');
        w_L_mat{ii} = repmat(dataA.w_L_list,  size(bA_mat_xc{ii},1), 1);
    end

    % xc导数
    dbA_dxc_mat = ( bA_mat_xc{2} - bA_mat_xc{1} ) ./ dataA.dxc_list;   % nT/cm
    % sigma_abs导数
    rel_dbA_dsigma_mat = ( bA_mat_sigma{2} - bA_mat_sigma{1} ) / dataA.rel_dsigma;   % sigma_abs * dbA/dsigma_abs, in unit of nT
    % Rrel导数
    rel_dbA_dRrel_mat = ( bA_mat_Rrel{2} - bA_mat_Rrel{1} ) / dataB.rel_dGamma;   % Rrel * dbA/dRrel, in unit of nT
        
    data = dataA;
    data.P_beam_mat = P_beam_mat{1};
    data.w_L_mat = w_L_mat{1};
    data.dbA_dxc_mat = dbA_dxc_mat;
    data.rel_dbA_dsigma_mat = rel_dbA_dsigma_mat;
    data.rel_dbA_dRrel_mat = rel_dbA_dRrel_mat;
    data.I_star = dataA.result_xc{1}.characteristic_info.I_star;

end





%% aux funs
function data = constructData1(data)
    P_beam_mat = {};
    w_L_mat = {};
    R_mat = {};
    bA_mat = {};

    g129 = data.cellPars_129.gXe;
    g131 = data.cellPars_131.gXe;
    R0 = g129/g131;

    for iT = 1:length(data.cellTemp_list)
        tmp = cellfun(@(re)re.P_beam_list, data.result(:,iT), 'UniformOutput', false);
        P_beam_mat{iT} = cell2mat(tmp');
        tmp = cellfun(@(re)re.R, data.result(:,iT), 'UniformOutput', false);
        R_mat{iT} = cell2mat(tmp');
        bA_mat{iT} = -data.cellPars_129.B0*( abs(R_mat{iT}) - abs(R0))/(R0);
        w_L_mat{iT} = repmat(data.w_L_list,  size(bA_mat{iT},1), 1);
    end

    % 温度导数
    dbA_dT_mat = ( bA_mat{3} - bA_mat{1} ) / 2/data.dT;   % nT/K
    % 功率导数
    iT=2;
    dbA_dPeam_mat = nan(size(bA_mat{iT}));
    for ii = 1:size(bA_mat{iT}, 2)
        bA_fun = griddedInterpolant(P_beam_mat{iT}(:,ii), bA_mat{iT}(:,ii), 'spline');
        bA_chebfun = chebfun(@(x)bA_fun(x), minmax(P_beam_mat{iT}(:,ii)'));
        dbA_dPeam_chebfun = diff(bA_chebfun,1);
        dbA_dPeam_mat(:,ii) = dbA_dPeam_chebfun(P_beam_mat{iT}(:,ii));    % nT/W
    end
    
    data.P_beam_mat = P_beam_mat;
    data.w_L_mat = w_L_mat;
    data.R_mat = R_mat;
    data.bA_mat = bA_mat;
    data.dbA_dPeam_mat = dbA_dPeam_mat;
    data.dbA_dT_mat = dbA_dT_mat;
    data.I_star = data.result{1}.characteristic_info.I_star;

end


