function [total_data, node_data, orig_subspace] = create_data_RSR(numb_of_inliers, numb_of_outliers, numb_of_nodes, min_inliers, min_outliers,  D, d, eps_var)
    % This function creates the dataset for the RSR and splits it between
    % processors
    % numb_of_inliers - total number of inliers, numb_of_outliers - total number of outliers
    % numb_of_nodes -number of nodes for the network, min_inliers - minimum number of inliers per node
    % min_outliers - minimum number of outliers per note 
    % D - ambient space dimension, d - the dimension of the linear subspace, 
    % eps_var - variance of the inliers around the linear subspace
    
    if(nargin < 8)
        eps_var = 0;
    end
    inliers_per_node = divide_number(numb_of_inliers, numb_of_nodes, min_inliers);
    outliers_per_node = divide_number(numb_of_outliers, numb_of_nodes, min_outliers);
    
    [x, ~] = svd(rand(D));
    orig_subspace = x(:, 1:d);
    total_inliers = zeros(D, numb_of_inliers);
    
    for j = 1:numb_of_inliers
        for k = 1:d
            total_inliers(:, j) = total_inliers(:, j) + randn * orig_subspace(:, k);
        end
    end
    total_inliers = total_inliers + eps_var * randn(size(total_inliers));
    
% Uniform outliers    
    total_outliers = 2 * rand(D, numb_of_outliers);

% Gaussian outliers
    node_data{1} = [total_inliers(:, 1:inliers_per_node(1)) total_outliers(:, 1:outliers_per_node(1))];
    total_data = node_data{1};
    cur_int_in = inliers_per_node(1);
    cur_int_out = outliers_per_node(1);
    for i = 2:numb_of_nodes
        node_data{i} = [total_inliers(:, (cur_int_in + 1):(cur_int_in + inliers_per_node(i))) ...
        total_outliers(:, (cur_int_out + 1):(cur_int_out + outliers_per_node(i)))];
        cur_int_in = cur_int_in + inliers_per_node(i); 
        cur_int_out = cur_int_out + outliers_per_node(i);
        total_data = [total_data node_data{i}];
    end
end

function [div_res] = divide_number(nb, count, min_nb)
    f = rand(1, count);
    f = f./ sum(f); % sum(f) = 1, approx.
    div_res = round(f .* nb);
    div_res(div_res < min_nb) = min_nb; % adjust the small numbers with min_data_points_node
    div_res(1) = nb - sum(div_res(2:end)); % deal with rounding errors etc.
    while(div_res(1) < min_nb)
        [max_val, max_ind] = max(div_res);
        div_res(max_ind) = min_nb;
        div_res(1) = div_res(1) + max_val - min_nb;
    end
end
