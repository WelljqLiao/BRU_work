clc, clear all
load .\var\Model_100cycle_output.mat
load .\var\Landcover_2020.mat

mycolorpoint=[[253 231 36];...
    [91 200 98];
    [32 143 140];
    [28 82 139];
    [68 1 84]];

mycolorposition=[1 32 64 96 128];
mycolormap_r=interp1(mycolorposition,mycolorpoint(:,1),1:128,'linear','extrap');
mycolormap_g=interp1(mycolorposition,mycolorpoint(:,2),1:128,'linear','extrap');
mycolormap_b=interp1(mycolorposition,mycolorpoint(:,3),1:128,'linear','extrap');
mycolor=[mycolormap_r',mycolormap_g',mycolormap_b']/255;
mycolor=round(mycolor*10^4)/10^4;

linecolor = [[192,0,0];
    [33,115,70];
    [43,87,154]]/255;

VarName = {'N Fixation','N Resorption','N uptake'};

%% CV-uncertainty
anss = prctile(BNF_cv ,[5,95],'all')
meanY = mean(BNF_cv,"all","omitnan");
BNF_cv(BNF_cv <= anss(1)) = meanY;
BNF_cv(BNF_cv >= anss(2)) = meanY;

anss = prctile(Nup_cv ,[5,95],'all')
meanY = mean(Nup_cv,"all","omitnan");
Nup_cv(Nup_cv <= anss(1)) = meanY;
Nup_cv(Nup_cv >= anss(2)) = meanY;

BNF_cv(Landcover_2020 <1 | Landcover_2020 >14) = nan;
NRE_cv(Landcover_2020 <1 | Landcover_2020 >14) = nan;
Nup_cv(Landcover_2020 <1 | Landcover_2020 >14) = nan;

BNF_cv = 100*flipud(BNF_cv)/30;
NRE_cv = 100*flipud(NRE_cv)/10;
Nup_cv = 100*flipud(Nup_cv)/200;

input = zeros(360,720,3);
input(:,:,1) = BNF_cv;
input(:,:,2) = NRE_cv;
input(:,:,3) = Nup_cv;

%% 100cycles CV
for i = 1:3
    mapdata = input(:,:,i);

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

    title(VarName{i},'FontName','Times', 'FontSize',18,...
        'Units','normalized','Position', [0.5, 1.05])
    h1 = colorbar('FontName', 'Times', 'FontSize', 10,'TickLabels','');
    ylabel(h1,'%', 'FontName', 'Times New Roman', 'FontSize', 12, ...
        'Rotation', 0,'Units','normalized','Position', [0.5, 1.08, 0.5]);

    mycolorpoint=[[0 78 203];...
        [255 255 255];
        [192 0 0]];
    mycolorposition=[1 64 128];
    mycolormap_r=interp1(mycolorposition,mycolorpoint(:,1),1:128,'linear','extrap');
    mycolormap_g=interp1(mycolorposition,mycolorpoint(:,2),1:128,'linear','extrap');
    mycolormap_b=interp1(mycolorposition,mycolorpoint(:,3),1:128,'linear','extrap');
    redblue =[mycolormap_r',mycolormap_g',mycolormap_b']/255;
    redblue =round(redblue*10^4)/10^4;

    colormap(redblue)
    colorbar
    caxis([0,1]);
    %     switch i
    %         case 1
    %             caxis([0,30]);
    %         case 2
    %             caxis([0,10]);
    %         case 3
    %             caxis([0,200]);
    %     end
    axis off;
    tightmap;
end

%% global
load .\var\leaf_litterN.mat
load .\var\Model_100cycle_output.mat
NRE_predict = NRE_predict*0.01;
resorp_N = (leaf_litterN.*NRE_predict)./(1-NRE_predict);
total_N = BNF_predict + Nup_predict + resorp_N;
load .\var\Area_WGS_1984_720_360.mat  %m2
Area = Area_WGS_1984/10000;

area_Nup = Nup_predict.*Area;
total_Nup = sum(area_Nup,'all','omitnan');
total_Nup = total_Nup*1000*1e-12;
result(1) = total_Nup

area_BNF = BNF_predict.*Area;
total_BNF = sum(area_BNF,'all','omitnan');
total_BNF = total_BNF*1000*1e-12;
result(2) = total_BNF;

area_NRE = resorp_N .*Area;
total_NRE = sum(area_NRE,'all','omitnan');
total_NRE= total_NRE*1000*1e-12;
result(3) = total_NRE;

name = {'Uptake';'Fixation';'Resorption'};
Global_Nacq = table(result','RowNames',name,'VariableNames',{'N acquisition(Tg/yr)'});
disp(Global_Nacq)
sum(result)

%% latitude pattern
lat = linspace(90, -90, 360);
plotsum = {area_Nup,area_BNF,area_NRE};

low_lat_idx = find(lat >= -23.5 & lat <= 23.5);
mid_lat_idx = find((lat > -66.5 & lat < -23.5) | (lat > 23.5 & lat < 66.5));
high_lat_idx = find(abs(lat) >= 66.5);

for i = 1:3

    lat_data = plotsum{i};
    low_lat_map = lat_data(low_lat_idx, :);
    mid_lat_map = lat_data(mid_lat_idx, :);
    high_lat_map = lat_data(high_lat_idx, :);

    Nacqsum(i,1) = sum(low_lat_map,'all','omitnan');
    Nacqsum(i,2) = sum(mid_lat_map,'all','omitnan');
    Nacqsum(i,3) = sum(high_lat_map,'all','omitnan');
end
Nacqsum = Nacqsum*1000*1e-12

colNames = {'Low','Mid','High'};
T = array2table(Nacqsum,'VariableNames', colNames,'RowNames', name);
disp(T)

labels = {'Uptake','Fixation','Resorption'};
figure('Position', [100, 100, 800, 600]);
b = bar(Nacqsum');
% ylim([-23 0])

for i = 1:3
    xtips1 = b(i).XEndPoints;
    ytips1 = b(i).YEndPoints;
    labels1 = string(b(i).YData);
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','baseline')
end

legend(labels,'Location','northeast');
set(findobj(gca,'Type','text'),'FontSize',15);
xticklabels({'Tropical° ','Temperature','Boreal'});
xlabel('Latitude')
ylabel('N acquisition amount,Tg/yr');

set(findall(gcf,'-property','FontName'),'FontName','Times');
set(findobj(gcf,'Type','text'),'FontSize',10);

%% 计算每种策略的占比
BNF_frac = 100*(BNF_predict./total_N);
NRE_frac = 100*(resorp_N./total_N);
Nup_frac = 100*(Nup_predict./total_N);
% save BRU_frac_100 BNF_predict resorp_N Nup_predict BNF_frac NRE_frac Nup_frac total_N

% case 1
BNF_frac = flipud(BNF_frac);
NRE_frac = flipud(NRE_frac);
Nup_frac = flipud(Nup_frac);

% input = zeros(360,720,3);
% input(:,:,1) = BNF_frac;
% input(:,:,2) = NRE_frac;
% input(:,:,3) = Nup_frac;
%
% globalmean = nanmean(input,1);
% globalmean = nanmean(globalmean,2);
%
% labels = {'Fixation','Resorption','Uptake'}
% figure('Position', [100, 100, 800, 600]);
% pie(globalmean);
% colormap(linecolor);
% legend(labels,'Location','southeastoutside');
% title('Globe','FontSize',18);
% set(findobj(gca,'Type','text'),'FontSize',15);
% set(findall(gcf,'-property','FontName'),'FontName','Times');

% % case 2
BNF_predict = flipud(BNF_predict);
resorp_N = flipud(resorp_N);
Nup_predict = flipud(Nup_predict);

input = zeros(360,720,3);
input(:,:,1) = BNF_predict;
input(:,:,2) = resorp_N;
input(:,:,3) = Nup_predict;

% Global pattern
close all
for i = 1:3
    mapdata = input(:,:,i);

    figure()
    lat = [-89.5:0.5:90];
    lon = [-179.5:0.5:180];
    [lon, lat] = meshgrid(lon,lat);
    gca = axesm('MapProjection','robinson','MapLatLimit',[-90 90],'Frame','off','Grid','on', ...
        'FontName','Times','FontSize',24,'FEdgeColor','white');
    setm(gca, 'FLineWidth', 0.5);
    setm(gca, 'GLineWidth', 0.2);
    setm(gca,'GLineStyle','-.');
    setm(gca, 'GColor', '#545454');
    i1 = surfm(lat, lon, mapdata);
    load coastlines
    plotm(coastlat,coastlon,'Color','k')
    h1 = colorbar('FontName', 'Times', 'FontSize', 14,'Location','southoutside',...
        'Position',[0.25 0.065 0.53 0.05]);
    colormap(mycolor)
    switch i
        case 1
            caxis([0,40]);
        case 2
            caxis([0,70]);
        case 3
            caxis([0,120]);
    end
    axis off;
    tightmap;
end

%% total N
mapdata = total_N;
mapdata = flipud(mapdata);

figure()
lat = [-89.5:0.5:90];
lon = [-179.5:0.5:180];
[lon, lat] = meshgrid(lon,lat);

gca = axesm('MapProjection','robinson','MapLatLimit',[-90 90],'Frame','off','Grid','on', ...
    'FontName','Times','FontSize',24,'FEdgeColor','white');
setm(gca, 'FLineWidth', 0.5);
setm(gca, 'GLineWidth', 0.2);
setm(gca,'GLineStyle','-.');
setm(gca, 'GColor', '#545454');

i1 = surfm(lat, lon, mapdata);
load coastlines
plotm(coastlat,coastlon,'Color','k')

h1 = colorbar('FontName', 'Times', 'FontSize', 14,'Location','southoutside',...
    'Position',[0.25 0.065 0.53 0.05]);
colormap(mycolor)
caxis([0,200])
axis off;
tightmap;

%% Latitude pattern
for i = 1:3
    mapdata = input(:,:,i);
    mapdata = flipud(mapdata);

    figure()
    data_mean = mean(mapdata,2,'omitnan');
    data_sd = nanstd(mapdata, 0, 2);
    lat = [89.75:-0.5:-89.75];

    hold on
    e = errorbar(data_mean, lat, data_sd, 'horizontal', 'Color', linecolor(i,:));
    e.CapSize = 0;
    e.LineWidth = 2;
    line(data_mean,lat,'LineWidth',2,'Color','k');

    set(gca, 'XAxisLocation', 'bottom', 'YAxisLocation', 'left');
    set(gca, 'Box', 'off', ...
        'LineWidth', 1,...
        'XMinorTick', 'off', 'YMinorTick', 'off', ...
        'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1])
    box on
    hLegend = legend(VarName{i},'Location', 'northeast');
    P = hLegend.Position;
    hLegend.AutoUpdate = 0;
    set(gca, 'YTick', -60:30:90, 'Ylim' ,[-60 90],'Yticklabel',{'60°S','30°S','0°','30°N','60°N','90°N'})
    set(gca, 'FontName', 'Times', 'FontSize', 15)
    set(gcf,'Color',[1 1 1])
    view([90, -90]);
    set(gca, 'YDir', 'reverse');

end
close all

%%
for i = 1:3
    mapdata = input(:,:,i);
    mapdata = flipud(mapdata);
    data_mean(:,i) = mean(mapdata,2,'omitnan');
end

figure;
b = bar(1:360, [data_mean(:,1), data_mean(:,2), data_mean(:,3)], 'stacked');
b(1).FaceColor = linecolor(1,:);
b(2).FaceColor = linecolor(2,:);
b(3).FaceColor = linecolor(3,:);
xlim([0,360]);
ylim([0,100]);
xticks(1:60:360);
xticklabels({'90°N','60°N','30°N','0°','30°S','60°S','90°S'});
xlabel('Latitude');
ylabel('Percentage (%)');
legend('Fixation', 'Resorption', 'Uptake', 'Location', 'southeast');
title('Stacked Bar Chart with Cumulative Percentage');
set(gca, 'FontName', 'Times', 'FontSize', 15)

close all
%%
for i = 1:3
    subplot(1,3,i)
    mapdata = input(:,:,i);
    histogram(mapdata,'FaceAlpha',0.8,FaceColor= linecolor(i,:))
    xlabel('Percentage (%)');
    title(VarName{i})
    set(gca, 'FontName', 'Times', 'FontSize', 15)
end
close all

%% ecosystem
fracdata = {'BNF_frac','NRE_frac','Nup_frac'};
% landcover
for n = 1:17
    cover_idx = find(Landcover_2020 == n);
    for i = 1:3
        Eco_data = eval(fracdata{i});
        cover_map = Eco_data(cover_idx);
        Ecomean(i,n) = nanmean(cover_map,'all');
    end
end
sum(Ecomean,1)
nanmean(Ecomean,2)
colNames = {'ENF','EBF','DNF','DBF','MF','CS','OS',...
    'WSava','Sava','Grass','Perma','Crop','Urban','Crop&Vet','Snow','Barren','Unclass'};
EcoT = array2table(Ecomean, 'VariableNames', colNames,'RowNames', fracdata);
disp(EcoT)

% pie
Ecomean(:,[3,17])= [];
colNames([3,17]) = [];
labels = {'Fixation','Resorption','Uptake'}
figure('Position', [100, 100, 800, 600]);
for i = 1:size(Ecomean, 2)
    subplot(3,5,i);
    pie(Ecomean(:, i));
    colormap(linecolor);
    title(colNames{i},'FontSize',18);
    set(findobj(gca,'Type','text'),'FontSize',15);
end
set(findall(gcf,'-property','FontName'),'FontName','Times');

close all

%% ecosystem
load landtype.mat
clear Ecomean
fracdata = {'BNF_frac','NRE_frac','Nup_frac'};
landtype = {'Forest','Shrub','Grass','Crop'};

for n = 1:4
    landtype_id = eval(landtype{n});
    cover_idx = find(landtype_id > 0);
    for i = 1:3
        Eco_data = eval(fracdata{i});
        Eco_data = reshape(Eco_data,[259200,1]);
        cover_map = Eco_data(cover_idx);
        Ecomean(i,n) = nanmean(cover_map,'all');
    end
end
sum(Ecomean,1)
nanmean(Ecomean,2)
EcoT = array2table(Ecomean, 'VariableNames', landtype,'RowNames', fracdata);
disp(EcoT)

%% climate zone
lat = linspace(90, -90, 360);
fracdata = {'BNF_frac','NRE_frac','Nup_frac'};
low_lat_idx = find(lat >= -30 & lat <= 30);
mid_lat_idx = find((lat > -60 & lat < -30) | (lat > 30 & lat < 60));
high_lat_idx = find(abs(lat) >= 60);

for i = 1:3
    lat_data = eval(fracdata{i});
    low_lat_map = lat_data(low_lat_idx, :);
    mid_lat_map = lat_data(mid_lat_idx, :);
    high_lat_map = lat_data(high_lat_idx, :);
    fracmean(i,1) = nanmean(low_lat_map, 'all');
    fracmean(i,2) = nanmean(mid_lat_map, 'all');
    fracmean(i,3) = nanmean(high_lat_map, 'all');
end
sum(fracmean,1)
colNames = {'Low','Mid','High'};
T = array2table(fracmean, 'VariableNames', colNames,'RowNames', fracdata);
disp(T)

% pie
labels = {'Fixation','Resorption','Uptake'}
figure('Position', [100, 100, 800, 600]);
for i = 1:size(fracmean, 2)
    subplot(1,3, i);
    pie(fracmean(:, i));
    colormap(linecolor); %
    title(colNames{i},'FontSize',18);
    set(findobj(gca,'Type','text'),'FontSize',15);
end
set(findall(gcf,'-property','FontName'),'FontName','Times');

close all

