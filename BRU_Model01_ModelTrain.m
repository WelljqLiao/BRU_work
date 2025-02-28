%% ------- RF model training ------- %%
% -*- coding: GBK -*-
% Created on Feb 20 2025 by Jiaqiang Liao
% To create the final community wood density maps, we used an ensemble approach,
% whereby we averaged the global predictions from the 100 best random-forest models
% based on our bootstrap procedure.

clc,clear all

Variable_Name = {'MAT','MAT_season','MAP','MAP_season','AI','AET','VPD','srad','tmax','tmin',...
    'CEC','BD','pH','Sand','Silt','Clay','SWC','SOC','TN','TP',...
    'MBC','MBN','FB_ratio','NPP','BNPP','NDVI','Ndep',...
    'LDMC','SLA_MODIS','LNC_MODIS','LPC_MODIS','Vcmax','fLNR',...
    'EM_tree','AM_tree'};
Variable_Name = string(Variable_Name)';
variablesNum = height(Variable_Name);


pool = gcp('nocreate');
if isempty(pool)
    pool = parpool();
end

%% ------- 1 choose dataset ------ %%
% SNF
BNF_ob = readtable('\data\SNFdata_extract.csv');
BNF = table2array(BNF_ob);
BNF = fillmissing(BNF,"movmean",60);
X = BNF(:,5:39);
Y = log(BNF(:,4));

%% NRE
NRE_ob = readtable('\data\NRE_extract.csv');
NRE = table2array(NRE_ob);
NRE = fillmissing(NRE,"movmean",1857);
X = NRE(:,5:39);
Y = NRE(:,4);

%% Nup
Nup_ob = readtable('\data\Nup_extract.csv');
Nup = table2array(Nup_ob);
Nup = fillmissing(Nup,"movmean",922);
X= Nup(:,5:39);
Y = log(Nup(:,4));

%% 1.1 特征选择

for i = 1:100
    % 训练集测试集划分
    cv = cvpartition(size(X, 1), 'HoldOut', 0.2);
    idxTrain = training(cv);
    Xtrain = X(idxTrain,:);
    ytrain = Y(idxTrain,:);
    Xtest = X(~idxTrain,:);
    ytest = Y(~idxTrain,:);

    % 建立随机森林模型
    numTrees = 50;
    RFModel = TreeBagger(numTrees, Xtrain, ytrain, 'Method', 'regression','OOBPredictorImportance','on');

    % 模型验证
    ypred1 = predict(RFModel, Xtest);
    mse = mean((ytest - ypred1).^2);
    rsquared = 1 - mse/var(ytest);

    r2(i,2) = rsquared;
    r2(i,1) = i;

    % 显示结果
    disp(['run: ', num2str(i)]);
    disp(['R-squared: ', num2str(rsquared)]);

    % 特征选择（这里使用特征重要性作为选择依据）
    importance = RFModel.OOBPermutedPredictorDeltaError;
    [~,idx] = sort(importance,'descend');
    numFeatures = 10; % 选择重要性最高的前10个特征
    selectedFeatures(i,:) = idx(1:numFeatures);

    % 保存特征的重要性信息
    VarIm(i,:) = RFModel.OOBPermutedPredictorDeltaError;

    % 保存全因子模型
    % Modelname = ['D:\研究生学习\氮循环\开搞BRU\1-code\SNFmodel\全因子RFModel_',num2str(i),'.mat'];
    % Modelname = ['D:\研究生学习\氮循环\开搞BRU\1-code\NREmodel\全因子RFModel_',num2str(i),'.mat'];
    % Modelname = ['D:\研究生学习\氮循环\开搞BRU\1-code\Nupmodel\全因子RFModel_',num2str(i),'.mat'];
    % save (Modelname,'RFModel')
end
[maxValue, index] = max(r2(:,2))

% 筛选出最优的十个特征
[counts,binEdges] = histcounts(selectedFeatures);
[sortedCounts, idx] = sort(counts, 'descend');
Features = idx(1:10)

% save Nupmodel\最优模型\Feature selectedFeatures Features r2

%% 1.2 基于10个特征的最优模型寻找
% load SNFmodel\最优模型\Feature.mat
% load NREmodel\最优模型\Feature.mat
% load Nupmodel\最优模型\Feature.mat

bestModels = cell(100, 1); % 用于保存最优模型的单元格数组
bestR2s = zeros(100, 1); % 用于保存最优模型的R-squared值

tic;    % 开始计时
parfor j = 1:100
    bestR2 = -inf; % 初始最优R-squared值为负无穷
    bestModel = []; % 初始最优模型为空

    for i = 1:100
        % 训练集测试集划分
        cv = cvpartition(size(X, 1), 'HoldOut', 0.2);
        idxTrain = training(cv);
        Xtrain = X(idxTrain,:);
        ytrain = Y(idxTrain,:);
        Xtest = X(~idxTrain,:);
        ytest = Y(~idxTrain,:);

        numTrees = 50;
        RFModelSelected = TreeBagger(numTrees, Xtrain(:, Features), ytrain, ...
            'Method', 'regression','OOBPredictorImportance','on');

        % 模型验证
        ypred2 = predict(RFModelSelected, Xtest(:, Features));
        mse = mean((ytest - ypred2).^2);
        rsquared = 1 - mse/var(ytest);

        % 保存最优模型和R-squared值
        if rsquared > bestR2
            bestR2 = rsquared;
            bestModel = RFModelSelected;
        end
    end

    % 保存最优模型和R-squared值到结果数组中
    bestModels{j} = bestModel;
    bestR2s(j) = bestR2;

    % 显示每次最优模型的R-squared值
    disp(['Best R-squared for iteration ', num2str(j), ': ', num2str(bestR2)]);
end
toc;

%%
% 保存所有最优模型
for j = 1:100
    Modelname = ['\1-code\Nupmodel\BestModel_',num2str(j),'.mat'];
    modelsave = bestModels{j};
    save(Modelname, 'modelsave');
end

% 找到最优的R-squared值和对应的模型索引
[maxR2, maxIndex] = max(bestR2s);
bestModel = bestModels{maxIndex};
% 显示最优的R-squared值和对应的模型索引
disp(['Best overall R-squared: ', num2str(maxR2)]);
disp(['Index of best model: ', num2str(maxIndex)]);

xlswrite("Nupmodel\bestR2Index.csv",bestR2s)

% % 关闭并行池
% delete(pool);
