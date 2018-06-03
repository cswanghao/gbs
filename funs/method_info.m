function method_info(choice_graph, choice_metric)

if 1 == choice_graph % complete graph
    if 2 == choice_metric
        fprintf('Method: GBS-CC --> ');
    elseif 3 == choice_metric
        fprintf('Method: GBS-CG --> ');
    end
elseif 2 == choice_graph % k-nearest graph
    if 1 == choice_metric
        fprintf('Method: GBS-KB --> ');
    elseif 2 == choice_metric
        fprintf('Method: GBS-KC --> ');
    elseif 3 == choice_metric
        fprintf('Method: GBS-KG --> ');
    elseif 4 == choice_metric
       fprintf('Method: GBS-KO --> ');
    end
end