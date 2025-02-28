%% 1. 计算氮获取占比后的三元相映射图
clc, clear all
close all
load .\var\Landcover_2020.mat
load .\var\BRU_frac_100.mat

BNF_frac(Landcover_2020 <1 |Landcover_2020 > 14) = 0;
NRE_frac(Landcover_2020 <1 |Landcover_2020 > 14) = 0;
Nup_frac(Landcover_2020 <1 |Landcover_2020 > 14) = 0;
filled_BNF = fillmissing(BNF_frac, 'movmean', 5); 
filled_NRE = fillmissing(NRE_frac, 'movmean', 5); 
filled_Nup = fillmissing(Nup_frac, 'movmean', 5); 
filled_BNF(Landcover_2020 <1 |Landcover_2020 > 14) = nan;
filled_NRE(Landcover_2020 <1 |Landcover_2020 > 14) = nan;
filled_Nup(Landcover_2020 <1 |Landcover_2020 > 14) = nan;
BNF_frac = flipud(filled_BNF);
NRE_frac = flipud(filled_NRE);
Nup_frac = flipud(filled_Nup);
%% 绘图
close all
figure('Position', [50, 100, 1200, 600])
% 创建地理坐标轴
% ax1 = worldmap('World'); 
gca = axesm('MapProjection','robinson','Frame','off')
% 定义纬度和经度的范围
latlim = [-90, 90]; % 纬度范围，从南极到北极
lonlim = [-180, 180]; % 经度范围，从西经到东经
% 创建地理参考对象
R = georefcells(latlim, lonlim, size(Landcover_2020));
% 使用 geoshow 函数显示栅格数据
BRUmap = geoshow(Landcover_2020, R); % grid 是栅格数据矩阵，R 是地理参考对象

% % 调用工具函数生成图例和映射表
A=BNF_frac;
B=NRE_frac;
C=Nup_frac;

[CMapData,CMapHdl]=multiVarMapTri(A,B,C,'colorList',1,'pieceNum',40);
CMapData = reshape(CMapData,[259200 3]);
A = reshape(A,[259200,1]);
nanMask = isnan(A);
CMapData(nanMask,:) = 1;  % 设置NAN值为浅灰色
% B= flipud(Landcover_2020);
% B = reshape(B,[259200,1]);
% landMask = find(B <1);
% CMapData(landMask,:) = 1;
CMapData = reshape(CMapData,[360 720 3]);
BRUmap.CData = CMapData;
gca = axesm('MapProjection','robinson','MapLatLimit',[-90 90],'Frame','off','Grid','on', ...
    'FontName','Times','FontSize',24,'FEdgeColor','white');
% setm(gca, 'FLineWidth', 1);
setm(gca, 'GLineWidth', 0.2);
setm(gca,'GLineStyle','-.');
setm(gca, 'GColor', '#545454');  
load coastlines
plotm(coastlat,coastlon,'Color','k')
tightmap;
axis off; 

%% 绘图方案二
close all
figure('Position', [50, 100, 1200, 600])

gca = axesm('MapProjection','eqdcylin','Frame','off')
latlim = [-90, 90];
lonlim = [-180, 180]; 
R = georefcells(latlim, lonlim, size(Landcover_2020));
BRUmap = geoshow(Landcover_2020, R);
A=BNF_frac;
B=NRE_frac;
C=Nup_frac;
[CMapData,CMapHdl]=multiVarMapTri(A,B,C,'colorList',1,'pieceNum',40);
CMapData = reshape(CMapData,[259200 3]);
A = reshape(A,[259200,1]);
nanMask = isnan(A);
CMapData(nanMask,:) = 1;
CMapData = reshape(CMapData,[360 720 3]);
BRUmap.CData = CMapData;

gca = axesm('MapProjection','eqdcylin','MapLatLimit',[-60 90],'Frame','on','Grid','off', ...
    'FontName','Times','FontSize',24,'FEdgeColor','k');
setm(gca, 'FLineWidth', 2);
mlabel('MLineLocation', -180:60:180, 'MLabelParallel','south','MlabelLocation',-180:60:180); 
plabel('PLineLocation', -60:30:90, 'PLabelMeridian','west','PlabelLocation',-60:30:90); 
child1 = gca.Children;
set(child1([7,8,9,10,11,12,13]), 'VerticalAlignment', 'baseline')
load coastlines
plotm(coastlat,coastlon,'Color','k')
tightmap;
axis off; 



