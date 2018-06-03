%
% GBS: Graph-based System for Multi-view Clustering
%%
clc;  close all; clear all;

%%========================================================================
%% settings (Choice method here)
%
runtimes = 1; % run-times on each dataset, default: 1
choice_graph = 2; % 1: 'Complete', and 2: 'k-nearest'
choice_metric = 4; % 1: 'Binary', 2: 'Cosine', 3: 'Gaussina-kernel', and 4: 'Our-method'
lambda = 1; %  initial parameter, which is tuned automatically
%
%%
graph = {'Complete','k-nearest'};
metric = {'Binary', 'Cosine', 'Gaussian-kernel', 'Our-method'};
dataname = {'100leaves','3sources','BBC','BBCSport','HW','HW2sources','NGs','WebKB'};
numdata = length(dataname);



currentFolder = pwd;
addpath(genpath(currentFolder));
resultdir = 'Results/';
if(~exist('Results','file'))
    mkdir('Results');
    addpath(genpath('Results/'));
end

for cdata = 1:numdata
%% read dataset
disp(char(dataname(cdata)));
datadir = 'Dataset/';
dataf = [datadir, cell2mat(dataname(cdata))];
load(dataf);

X = data;
y0 = truelabel{1};
c = length(unique(truelabel{1}));
%% iter ...
for rtimes = 1:runtimes

% Multi-view clustering method on graph-based system (GBS)
lambda = 1;
[F, y, U, S0, evs] = GBS(X, c, choice_graph, choice_metric, lambda); % c: the # of clusters
measure = CalcMeasures(y0, y);
ACC(rtimes) = measure(1);
NMI(rtimes) = measure(2);
error_cnt(rtimes) = measure(4);
fprintf('...Runtime %d> ACC:%.4f\tNMI:%.4f\terror_cnt:%d\n',rtimes,measure(1),measure(2),measure(4));
end;
    Result(1,:) = ACC;
    Result(2,:) = NMI;
    Result(4,1) = mean(ACC);
    Result(4,2) = mean(NMI);
    Result(5,1) = std(ACC);
    Result(5,2) = std(NMI);
save([resultdir,char(dataname(cdata)),'_result.mat'],'Result','U','y0','y');
clear ACC NMI measure Result U y0 y;
end;
