%% 按所有变量提取点
clc, clear all
close all
load .\var\Landcover_2020.mat
load .\var\X_predict(35)_inter.mat
load .\var\BRU_frac_100.mat
% 变量名
Variable_Name = {'MAT','MAT_season','MAP','MAP_season','AI','AET','VPD','srad','tmax','tmin',...
    'CEC','BD','pH','Sand','Silt','Clay','SWC','SOC','TN','TP',...
    'MBC','MBN','FB_ratio','NPP','BNPP','NDVI','Ndep',...
    'LDMC','SLA_MODIS','LNC_MODIS','LPC_MODIS','Vcmax','fLNR',...
    'EM_tree','AM_tree'};
Variable_Name = string(Variable_Name)';
variablesNum = height(Variable_Name);

BNF_frac = reshape(BNF_frac,[259200,1]);
NRE_frac = reshape(NRE_frac,[259200,1]);
Nup_frac = reshape(Nup_frac,[259200,1]);

landcover = reshape(Landcover_2020,[259200,1]);
landcover = double(landcover);

ResultTable = [MAT,MAP,BNF_frac,NRE_frac,Nup_frac];
ResultTable = [MAT,MAT_season,MAP,MAP_season,AI, ...
    AET,VPD,srad,tmax,tmin,CEC,BD,	pH,	Sand,Silt,Clay, ...
    SWC,SOC,TN,TP,	MBC,MBN,FB_ratio,NPP,BNPP,NDVI,Ndep, ...
    LDMC,SLA_MODIS,LNC_MODIS,LPC_MODIS,Vcmax,fLNR,EM_tree,AM_tree,BNF_frac,NRE_frac,Nup_frac,landcover];
% 检查第38列是否存在NaN值（Nup_frac)
nan_indices = isnan(ResultTable(:, 38));

% 删除包含NaN值的行
ResultTable_cleaned = ResultTable(~nan_indices, :);

Variable = {'MAT','MAT_season','MAP','MAP_season','AI','AET','VPD','srad','tmax','tmin',...
    'CEC','BD','pH','Sand','Silt','Clay','SWC','SOC','TN','TP',...
    'MBC','MBN','FB_ratio','NPP','BNPP','NDVI','Ndep',...
    'LDMC','SLA_MODIS','LNC_MODIS','LPC_MODIS','Vcmax','fLNR',...
    'EM_tree','AM_tree','BNF','NRE','Nup','LandID'};
Table = array2table(ResultTable_cleaned, 'VariableNames', Variable);
% writetable(Table,'PCAtable2.xlsx');

% save PCAdata_100 ResultTable_cleaned

%% 主成分分析
clear
load PCAdata_100.mat
pca_X = ResultTable_cleaned(:,[1:35]);
X_filled = fillmissing(pca_X, 'movmean', 49469);
X_normalized = zscore(X_filled);
[coeff, score, latent, ~, explained] = pca(X_normalized,'NumComponents',3,'VariableWeights','variance');
Variable = {'MAT','MAT_season','MAP','MAP_season','AI','AET','VPD','srad','tmax','tmin',...
    'CEC','BD','pH','Sand','Silt','Clay','SWC','SOC','TN','TP',...
    'MBC','MBN','FB_ratio','NPP','BNPP','NDVI','Ndep',...
    'LDMC','SLA_MODIS','LNC_MODIS','LPC_MODIS','Vcmax','fLNR',...
    'EM_tree','AM_tree','BNF','NRE','Nup','LandID',...
    'PCA1','PCA2','PCA3',};
PCA_shapelydata = [ResultTable_cleaned,score];
Table = array2table(PCA_shapelydata, 'VariableNames', Variable);
writetable(Table,'PCA_shapelydata.xlsx');

figure;
scatter(score(:, 1), score(:, 2));
xlabel('主成分1');
ylabel('主成分2');
title('前两个主成分的PCA图');

% 绘制主成分碎石图
figure;
plot(1:length(explained), explained, 'bo-', 'LineWidth', 2);
xlabel('主成分');
ylabel('方差解释比例（%）');
title('主成分碎石图');
grid on;
