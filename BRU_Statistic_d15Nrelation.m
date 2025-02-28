clc, clear all
close all
load .\var\mycolor.mat
load .\var\Landcover_2020.mat
load .\var\BRU_frac_100.mat
load .\var\Vegetation_data.mat LDMC
x = reshape(LDMC,[259200,1]);
y = reshape(total_N,[259200,1]);
df = [x,y];
nan_indices = isnan(df(:,1));
df_nonan = df(~nan_indices, :);
nan_indices = isnan(df_nonan(:,2));
df_nonan = df_nonan(~nan_indices, :);

BNF_frac = reshape(BNF_frac,[259200,1]);
NRE_frac = reshape(NRE_frac,[259200,1]);
Nup_frac = reshape(Nup_frac,[259200,1]);

load d15N_data.mat
close all
EF = N15_leaf_in - N15_soil;
figure(),imagesc(EF)
figure(),histogram(EF)

d15N = reshape(N15_soil,[259200,1]);
% d15N = reshape(N15_leaf_in,[259200,1]);
df = [BNF_frac,d15N,Nup_frac];
nan_indices = isnan(df(:,1));
df_nonan = df(~nan_indices, :);
nan_indices = isnan(df_nonan(:,2));
df_nonan = df_nonan(~nan_indices, :);

%% soil 15N
mapdata = N15_leaf_in;
mapdata = flipud(mapdata);

figure()
lat = [-89.5:0.5:90];
lon = [-179.5:0.5:180];
[lon, lat] = meshgrid(lon,lat);
ax1 = axesm('MapProjection','eqdcylin','MapLatLimit',[-60 90],'MapLonLimit',[-180 180],'Frame','on','Grid','on', ...
    'FontName','Times','FontSize',12,'FEdgeColor','none', ...
    'MLineLocation',60,'MLabelRound', 0, 'MeridianLabel','on',...
    'PLineLocation',30,'PLabelRound', 0,'ParallelLabel','on', ...
    'MLabelParallel','south');
i1 = surfm(lat, lon, mapdata);
load coastlines
plotm(coastlat,coastlon,'Color','k')
title('Leaf δ15N','FontName','Times', 'FontSize',18,...
    'Units','normalized','Position', [0.5, 1.05])
h1 = colorbar('FontName', 'Times', 'FontSize', 10);
ylabel(h1,'‰', 'FontName', 'Times New Roman', 'FontSize', 12, ...
    'Rotation', 0,'Units','normalized','Position', [0.5, 1.08, 0.5]);
colormap(mycolor)
axis off;
tightmap;

%% bin
close all
r = 1;
dailyT = d15N;
dailyET = Nup_frac;
minT = floor(min(dailyT));
maxT = ceil(max(dailyT));
bins = (maxT-minT-2);
for i = 1:bins
    T = minT + (i-1);
    datapoints = find(dailyT <= (T+5) & dailyT >= T);
    movingT(i) = nanmean(dailyT(datapoints));
    movingET(i) = nanmean(dailyET(datapoints));
    movingBNF(i) = nanmean(BNF_frac(datapoints));
    movingNRE(i) = nanmean(NRE_frac(datapoints));
    movingT_Se(i) = nanstd(dailyT(datapoints))/sqrt(length(datapoints));
    movingET_Se(i) = nanstd(dailyET(datapoints))/sqrt(length(datapoints));
    movingBNF_Se(i) = nanstd(BNF_frac(datapoints))/sqrt(length(datapoints));
    movingNRE_Se(i) = nanstd(NRE_frac(datapoints))/sqrt(length(datapoints));
end
moving(1:length(movingT),r) = movingT;
moving(1:length(movingET),(r+1)) = movingET;
moving(1:length(movingNRE),(r+2)) = movingNRE;
moving(1:length(movingBNF),(r+3)) = movingBNF;

%% 绘图
x = moving(:,1)';
y = moving(:,2:4)';
err = [movingET_Se;movingNRE_Se;movingBNF_Se];
err = 0.5*err;
linecolor = [[192,0,0];
    [33,115,70];
    [43,87,154]]/255;

figure()

for diu = 1:3
    X = double(x);
    Y = double(y(diu,:));
    hold on
    E(diu) = errorbar(X,Y,err(diu,:));


    n = length(X);
    p = 0.5;
    bootstrap_num = 1000;
    pp_set = cell(bootstrap_num, 1);
    for i = 1:bootstrap_num
        ind = randi(n, n, 1);
        pp_set{i} = csaps(X(ind), Y(ind), p);
    end

    Xfit = linspace(min(X), max(X), 1000);
    Yfit_set = zeros(bootstrap_num, length(Xfit));
    for i = 1:bootstrap_num
        Yfit_set(i, :) = ppval(pp_set{i}, Xfit);
    end
    Yfit_median = median(Yfit_set);
    Yfit_ci = quantile(Yfit_set, [0.01, 0.99]);

    % 画出平滑曲线和置信区间
    % plot(Xfit, Yfit_median, '-');
    C = linecolor(4-diu,:);
    ciplot(Yfit_ci(1,:),Yfit_ci(2,:),Xfit,C,0.5);
    hold on
    set(E(diu),  'LineStyle', '-', 'Color', C,...
        'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 6, ...
        'MarkerEdgeColor', [.2 .2 .2], 'MarkerFaceColor' , C)
end
hXLabel = xlabel('Soil δ^1^5N (‰)');
hYLabel = ylabel('Relative contribution (%)');
set(gca, 'Box', 'off', ...
    'LineWidth',1,...
    'XGrid', 'off', 'YGrid', 'off', ...
    'TickDir', 'out', 'TickLength', [.015 .015], ...
    'XMinorTick', 'off', 'YMinorTick', 'off', ...
    'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1],...
    'YTick', 0:20:100,...
    'YLim', [-20 100], ...
    'XTick', 1:2:14,...
    'XLim', [1 14],...
    'XDir','reverse')

XL = get(gca,'xlim'); XR = XL(1);
YL = get(gca,'ylim'); YT = YL(2);
xc = get(gca,'XColor');
yc = get(gca,'YColor');
plot(XL,YT*ones(size(XL)),'color', xc,'linewidth',1)
plot(XR*ones(size(YL)),YL,'color', yc,'linewidth',1)
hLegend = legend(E, ...
    'Uptake', 'Resorption','Fixation', ...
    'Location', 'northeast');
P = hLegend.Position;
set(gca, 'FontName', 'Times', 'FontSize', 16)
set([hXLabel, hYLabel], 'FontSize', 18, 'FontName', 'Times')
set([hLegend], 'FontSize', 12, 'FontName', 'Times')
set([hLegend], 'Color', 'none');
set(gca, 'color', 'none');
set(gcf, 'color', 'none');

close all
