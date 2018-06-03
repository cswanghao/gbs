%
% min sum_v w_v||U - Sv||^2 + lambda*trace(F'*Lu*F)
% s.t U>=0, Ui^T*1=1, F'*F=I
function [F, y, U, S0, evs] = GBS(X, c, choice_graph, choice_metric, lambda, normData)
%% Input
% X{}: multi-view dataset, each cell corresponds to a view, each column corresponds to a data point
% choice_graph: --> 1: 'Complete', and 2: 'k-nearest'
% Choice_metric: --> 1: 'Binary', 2: 'Cosine', 3: 'Gaussina-kernel', and 4: 'Our-method'
% lambda: initial parameter (default value is 1), which is tuned automatically
% normData: whether need to normalize data
%% Output
% F: the embedding matrix
% y: the final clustering result, i.e., cluster indicator vector
% U: the learned unified matrix
% S0: the constructed SIG matrix, each row corresponds to a data point
% evs: eigenvalues of learned graph Laplacian in the iterations
%%

method_info(choice_graph, choice_metric);
NITER = 200;
zr = 1e-10;
islocal = 1; % default: only update the similarities of neighbors if islocal=1
if nargin < 3
    choice_graph = 2; % suggest using k-nearest graph
end
if nargin < 4
    choice_metric = 4; % suggest using our method
end
if nargin < 5
    lambda = 1;
end
if nargin < 6
    normData = 1;
end;

num = size(X{1},2); % number of instances
m = length(X); % number of views
%% Normalization: Z-score
if normData == 1
    for i = 1:m
        for  j = 1:num
            normItem = std(X{i}(:,j));
            if (0 == normItem)
                normItem = eps;
            end;
            X{i}(:,j) = (X{i}(:,j)-mean(X{i}(:,j)))/(normItem);
        end;
    end;
end;

%% Constructing the SIG matrices
pn = 15; % pn: number of adaptive neighbours
options = [];
options.k = 5;

S0 = cell(1,m);
for i = 1:m
    if 1 == choice_graph % complete graph
        options.k = 0;
        if 1 == choice_metric
            options.WeightMode = 'Binary';
            S0{i} = constructS_KNG(X{i}', options);
        elseif 2 == choice_metric
            options.WeightMode = 'Cosine';
            S0{i} = constructS_KNG(X{i}', options);
        elseif 3 == choice_metric
            options.WeightMode = 'HeatKernel';
            S0{i} = constructS_KNG(X{i}', options);
        else
            error('Invalid input: check choice_metric');
        end
    elseif 2 == choice_graph % k-nearest graph
        if 1 == choice_metric
            options.WeightMode = 'Binary';
            S0{i} = constructS_KNG(X{i}', options);
        elseif 2 == choice_metric
            options.WeightMode = 'Cosine';
            S0{i} = constructS_KNG(X{i}', options);
        elseif 3 == choice_metric
            options.WeightMode = 'HeatKernel';
            S0{i} = constructS_KNG(X{i}', options);
        elseif 4 == choice_metric
            [S0{i}, distX_i] = constructS_PNG(X{i}, pn, 0);
        else
            error('Invalid input: check choice_metric');
        end
    else
        error('Invalid input: check choice_graph');
    end
end

%% initialize U, F and w
U0 = zeros(num);
for i = 1:m
    U0 = U0 + S0{i};
end
U0 = U0/m;
for j = 1:num
    d_sum = sum(U0(j,:));
    if d_sum == 0
        d_sum = eps;
    end
    U0(j,:) = U0(j,:)/d_sum;
end
U = (U0+U0')/2;

D = diag(sum(U));
L = D - U;
[F, ~, evs]=eig1(L, c, 0);

w = ones(1,m)/m;

%%  update ...
for iter = 1:NITER
    % calculate the objective value
    for v = 1:m
        tempF(v) = w(v)*norm(U - S0{v}, 'fro')^2;
    end
    fLf = F'*L*F;
    obj_value(iter) = sum(tempF) + lambda*trace(fLf);
    %obj_value = 0;
    % update W
    for v = 1:m
        US = U - S0{v};
        distUS = norm(US, 'fro')^2;
        if distUS == 0
            distUS = eps;
        end;
        w(v) = 0.5/sqrt(distUS);
    end
    % update U
    dist = L2_distance_1(F',F');
    U = zeros(num);
    for i=1:num
        idx = zeros();
        for v = 1:m
            s0 = S0{v}(i,:);
            idx = [idx,find(s0>0)];
        end
        idxs = unique(idx(2:end));
        if islocal == 1
            idxs0 = idxs;
        else
            idxs0 = 1:num;
        end;
        for v = 1:m
            s1 = S0{v}(i,:);
            si = s1(idxs0);
            di = dist(i,idxs0);
            mw = m*w(v);
            lmw = lambda/mw;
            q(v,:) = si-0.5*lmw*di;
        end
        U(i,idxs0) = SloutionForP20(q,m);
        clear q;
    end
    % choose the top-k neighbors
%     [~, ids] = sort(U,2,'descend');
%     ts = zeros(num);
%     for i =1:num
%         ts(i,ids(i,1:pn)) = U(i,ids(i,1:pn));
%     end
%     U = ts;
    % update F
    sU = U;
    sU = (sU+sU')/2;
    D = diag(sum(sU));
    L = D-sU;
    F_old = F;
    [F, ~, ev]=eig1(L, c, 0);
    evs(:,iter+1) = ev;
    % update lambda and the stopping criterion
    fn1 = sum(ev(1:c));
    fn2 = sum(ev(1:c+1));
    if fn1 > zr
        lambda = 2*lambda;
    elseif fn2 < zr
        lambda = lambda/2;
        F = F_old;
    else
        disp(['iter:',num2str(iter),' lambda:',num2str(lambda)]);
        break;
    end;
end;
%% generating the clustering result
%[labv, tem, y] = unique(round(0.1*round(1000*F)),'rows');
[clusternum, y]=graphconncomp(sparse(sU)); y = y';
if clusternum ~= c
    fprintf('Can not find the correct cluster number: %d\n', c)
end; 


