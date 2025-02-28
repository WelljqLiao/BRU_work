% -*- coding: GBK -*-
% Created on Feb 20 2025 by Jiaqiang Liao
clear all, clc

% 初始化并行计算环境
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool();
end

%导入预测数据：259200*1
load .\var\X_predict(35)_inter.mat
load .\var\Landcover_2020.mat
load .\var\mycolor.mat

% 变量名
Variable_Name = {'MAT','MAT_season','MAP','MAP_season','AI','AET','VPD','srad','tmax','tmin',...
    'CEC','BD','pH','Sand','Silt','Clay','SWC','SOC','TN','TP',...
    'MBC','MBN','FB_ratio','NPP','BNPP','NDVI','Ndep',...
    'LDMC','SLA_MODIS','LNC_MODIS','LPC_MODIS','Vcmax','fLNR',...
    'EM_tree','AM_tree'};
Variable_Name = string(Variable_Name)';
variablesNum = height(Variable_Name);

% % 批量reshape
% newSize = [259200,1];
%
% %
% for i = 1:numel(Variable_Name)
%     varName = Variable_Name{i};
%     if exist(varName, 'var') % 检查变量是否存在
%         % 使用 evalin 从工作区获取变量，然后使用 reshape 重塑
%         evalin('base', [varName, ' = reshape(', varName, ', newSize);']);
%     end
% end

%% -------  SNF ------- %%
load SNFmodel\Feature.mat

varNames = Variable_Name(Features)
predict_X = zeros(length(eval(varNames{1})), length(varNames));
for i = 1:length(varNames)
    predict_X(:, i) = eval(varNames{i});
end

%% 循环读取100次模型进行预测
Y_predict_all = zeros(size(predict_X, 1), 100);
Landcover_2020 = reshape(Landcover_2020,[259200 1]);
tic;
parfor j = 1:100
    Modelname = ['\1-code\SNFmodel\BestModel_',num2str(j),'.mat'];
    matObj = matfile(Modelname);
    modelsave = matObj.modelsave;
    Y_predict = predict(modelsave, predict_X);
    Y_predict(Landcover_2020 <1 | Landcover_2020 >14) = nan;
    Y_predict_all(:, j) = Y_predict;
end
toc;

Y_predict_sum = sum(Y_predict_all,2);
Y_predict_avg = Y_predict_sum/100;
Y_predict_std = std(Y_predict_all, 0, 2);
Y_predict_std = reshape(Y_predict_std,[360,720]);
anss = prctile(Y_predict_avg ,[5,95],'all')  % rate分位数
meanY = mean(Y_predict_avg,"all","omitnan"); % rate 均值
Y_predict_avg = reshape(Y_predict_avg,[360,720]);
Y_predict_avg(Landcover_2020 <1 | Landcover_2020 >14) = nan;
histogram(Y_predict_avg)
BNF_cv = Y_predict_std./Y_predict_avg;
Y_predict_avg = exp(Y_predict_avg);
histogram(Y_predict_avg)
BNF_predict = Y_predict_avg;

% global size
load .\var\Area_WGS_1984_720_360.mat
Area = Area_WGS_1984/10000;
Area = reshape(Area,[259200 1]);
SNF_rate = exp(Y_predict_all);
area_SNF = SNF_rate.*Area;
total_BNF = sum(area_SNF,1,'omitnan');
total_BNF = total_BNF*1000*1e-12;
anss = prctile(total_BNF,[5,95],'all')
meanY = mean(total_BNF,"all","omitnan")

%% -------  NRE ------- %%
load NREmodel\Feature.mat
varNames = Variable_Name(Features)
predict_X = zeros(length(eval(varNames{1})), length(varNames));
for i = 1:length(varNames)
    predict_X(:, i) = eval(varNames{i});
end

%% 循环读取100次模型
Y_predict_all = zeros(size(predict_X, 1), 100);
tic;
parfor j = 1:100
    Modelname = ['\1-code\NREmodel\BestModel_',num2str(j),'.mat'];
    matObj = matfile(Modelname);
    modelsave = matObj.modelsave;
    Y_predict = predict(modelsave, predict_X);
    Y_predict(Landcover_2020 <1 | Landcover_2020 >14) = nan;
    Y_predict_all(:, j) = Y_predict;
end
toc;

%
Y_predict_sum = sum(Y_predict_all,2);
Y_predict_avg = Y_predict_sum/100;
Y_predict_std = std(Y_predict_all, 0, 2);
Y_predict_std = reshape(Y_predict_std,[360,720]);
anss = prctile(Y_predict_avg ,[5,95],'all')
meanY = mean(Y_predict_avg,"all","omitnan")
Y_predict_avg = reshape(Y_predict_avg,[360,720]);
Y_predict_avg(Landcover_2020 <1 | Landcover_2020 >14) = nan;
histogram(Y_predict_avg)
NRE_predict = Y_predict_avg;
NRE_cv = Y_predict_std./Y_predict_avg;
load leaf_litterN.mat
leaf_litterN = reshape(leaf_litterN,[259200,1]);
Y_predict_all = Y_predict_all*0.01;
resorp_N = (leaf_litterN.*Y_predict_all)./(1-Y_predict_all);
load Area_WGS_1984_720_360.mat
Area = Area_WGS_1984/10000;
Area = reshape(Area,[259200 1]);
area_NRE = resorp_N.*Area;
total_NRE = sum(area_NRE,1,'omitnan');
total_NRE = total_NRE*1000*1e-12;
anss = prctile(total_NRE,[5,95],'all')
meanY = mean(total_NRE,"all","omitnan")


%% -------  Nup ------- %%
load Nupmodel\Feature.mat
varNames = Variable_Name(Features)
predict_X = zeros(length(eval(varNames{1})), length(varNames));
for i = 1:length(varNames)
    predict_X(:, i) = eval(varNames{i});
end

%% 循环读取100次模型
Y_predict_all = zeros(size(predict_X, 1), 100);

tic;
parfor j = 1:100
    Modelname = ['\1-code\Nupmodel\BestModel_',num2str(j),'.mat'];
    matObj = matfile(Modelname);
    modelsave = matObj.modelsave;
    Y_predict = predict(modelsave, predict_X);
    Y_predict_all(:, j) = Y_predict;
end
toc;

%
Y_predict_sum = sum(Y_predict_all,2);
Y_predict_avg = Y_predict_sum/100;
Y_predict_std = std(Y_predict_all, 0, 2);
Y_predict_std = reshape(Y_predict_std,[360,720]);
anss = prctile(Y_predict_avg ,[5,95],'all')
meanY = mean(Y_predict_avg,"all","omitnan");
Y_predict_avg(Y_predict_avg <= anss(1)) = meanY;
Y_predict_avg(Y_predict_avg >= anss(2)) = meanY;
Y_predict_avg = reshape(Y_predict_avg,[360,720]);
Y_predict_avg(Landcover_2020 <1 | Landcover_2020 >14) = nan;
Y_predict_avg = exp(Y_predict_avg);
Y_predict_std = exp(Y_predict_std);
Nup_cv = Y_predict_std./Y_predict_avg;
load BNPP_200.mat
BNPP_200 = reshape(BNPP_200,[360,720]);%0~20cm
%ug*g-1*h-1 to kg*ha-1*yr-1
Nup_predict = BNPP_200.*Y_predict_avg*1e-3*24*365;
histogram(Nup_predict)
load Area_WGS_1984_720_360.mat
Area = Area_WGS_1984/10000;
area_Nup = Nup_predict.*Area;
total_Nup = sum(area_Nup,'all','omitnan');
total_Nup = total_Nup*1000*1e-12
Y_predict_all(Y_predict_all <= anss(1)) = meanY;
Y_predict_all(Y_predict_all >= anss(2)) = meanY;
BNPP_200 = reshape(BNPP_200,[259200,1]);
Nup_100rate = BNPP_200.*exp(Y_predict_all)*1e-3*24*365;
Area = reshape(Area,[259200 1]);
area_Nup = Nup_100rate.*Area;
total_Nup = sum(area_Nup,1,'omitnan');
total_Nup = total_Nup*1000*1e-12;
anss = prctile(total_Nup,[5,95],'all')
meanY = mean(total_Nup,"all","omitnan")

%% save model global map result
% save Model_100cycle_output BNF_predict BNF_cv NRE_predict NRE_cv Nup_predict Nup_cv

