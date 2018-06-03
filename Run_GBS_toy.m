% Experiments on the two-moon dataset
% Multi-view clustering method on graph-based system (GBS)
%
%% Notes:
% We construct a toy dataset named two-moon dataset. In this dataset, there
% are two are two clusters with each having one hundred samples. We randomly
% generate two sets of two-moon data and each set is seen as a view
% respectively. Our goal is to verify that the number of connected
% components in the learnt unified matrix is exactly two and the data
% points are partitioned into two clusters.
%%

clc;  close all; clear all;
currentFolder = pwd;
addpath(genpath(currentFolder));
dataname = {'ToyData'};

newdata = 0; % if you need to generate a new toy dataset, set newdata = 1
if newdata == 1
    %% Data Generalization
    m = 2;
    data = cell(1,m);
    num0 = 100;
    for i = 1:m
        oneView = twomoon_gen(num0);
        data{i} = oneView';
        truelabel{i} = [ones(num0,1);2*ones(num0,1)];
    end
else
    %% read dataset
    idata = 1;
    datadir = 'Dataset/';
    dataf = [datadir, cell2mat(dataname(idata))];
    load(dataf);
end
X = data;
y0 = truelabel{1};
num = size(X{1},2); % number of instances
c = length(unique(y0)); % number of clusters
m = length(X); % number of views
choice_graph = 2; % 1: 'Complete', and 2: 'k-nearest'
choice_metric = 4; % 1: 'Binary', 2: 'Cosine', 3: 'Gaussina-kernel', and 4: 'Our-method'
lambda = 1;
normData = 0;

%% Multi-view clustering method on graph-based system (GBS)
[F, y, U, S0, evs] = GBS(X, c, choice_graph, choice_metric, lambda, normData); % c: the # of clusters
U2 = (U+U')/2;

%% Original data
for v = 1:m
    lab = y0;
    figure; 
    plot(X{v}(1,:),X{v}(2,:),'.k', 'MarkerSize', 20); hold on;
    plot(X{v}(1,lab==1),X{v}(2,lab==1),'.r', 'MarkerSize', 20); hold on;
    plot(X{v}(1,lab==2),X{v}(2,lab==2),'.', 'MarkerSize', 20); hold on;
%     set(gca,'xlim',[-1.7,1.7],'xtick',[-1.5:0.5:1.5]); % set x-axis
%     set(gca,'ylim',[-1.2,1.2],'ytick',[-1:0.5:1]); % set y-axis
    set(gca,'FontName','Times New Roman','FontSize',24,'LineWidth',1.5);
    axis equal;
end

%% Original connected graph with probabilistic neighbors, line width denotes similarity
for v = 1:m
    figure; 
    plot(X{v}(1,:),X{v}(2,:),'.k', 'MarkerSize', 20); hold on;
    plot(X{v}(1,lab==1),X{v}(2,lab==1),'.r', 'MarkerSize', 20); hold on;
    plot(X{v}(1,lab==2),X{v}(2,lab==2),'.', 'MarkerSize', 20); hold on;
    for ii = 1 : num;
        for jj = 1 : ii
            weight = S0{v}(ii, jj);
            if weight > 0
                plot([X{v}(1, ii), X{v}(1, jj)], [X{v}(2, ii), X{v}(2, jj)], '-g', 'LineWidth', 20*weight), hold on;
            end
        end;
    end;
%     set(gca,'xlim',[-1.7,1.7],'xtick',[-1.5:0.5:1.5]); % set x-axis
%     set(gca,'ylim',[-1.2,1.2],'ytick',[-1:0.5:1]); % set y-axis
    set(gca,'FontName','Times New Roman','FontSize',24,'LineWidth',1.5);
    axis equal;
end

%% the learned unified graph by GBS, line width denotes similarity
for v = 1:m
    figure; 
    plot(X{v}(1,:),X{v}(2,:),'.k', 'MarkerSize', 20); hold on;
    plot(X{v}(1,lab==1),X{v}(2,lab==1),'.r', 'MarkerSize', 20); hold on;
    plot(X{v}(1,lab==2),X{v}(2,lab==2),'.', 'MarkerSize', 20); hold on;
    for ii = 1 : num;
        for jj = 1 : ii
            weight = U2(ii, jj);
            if weight > 0
                plot([X{v}(1, ii), X{v}(1, jj)], [X{v}(2, ii), X{v}(2, jj)], '-g', 'LineWidth', 20*weight), hold on;
            end
        end;
    end;
%     set(gca,'xlim',[-1.7,1.7],'xtick',[-1.5:0.5:1.5]); % set x-axis
%     set(gca,'ylim',[-1.2,1.2],'ytick',[-1:0.5:1]); % set y-axis
    set(gca,'FontName','Times New Roman','FontSize',24,'LineWidth',1.5);
    axis equal;
end
