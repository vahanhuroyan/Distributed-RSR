function [ est_subspaces, Q_new ] = d_gms( node_data, nb_mat, d, step_size)
%   D_GMS Summary of this function goes here
    
% for high dimensional data the delta is very small, that's why sometimes
% we need large stepsize

    D = size(node_data{1}, 1);
    numb_of_nodes = size(node_data, 2);
    
    if(nargin < 4) 
        dlta = 1e-10;
        proc_numb = 1;        
        temp = sum(node_data{proc_numb}'.^2, 2);
        scaled_cov = (node_data{proc_numb}.*repmat(min(temp.^-0.5,1 / dlta), 1, D)' * node_data{proc_numb}');
        step_size = D^2 / (10 * (numb_of_nodes - 1) * trace(inv(scaled_cov))); 
    end
    
    if(nargin < 3)
        d = D;
    end
    if(nargin < 2)
        nb_mat = ones(numb_of_nodes);
    end
    
    Q_old = cell(numb_of_nodes, 1);%zeros(D, D, numb_of_nodes);
    
    for i = 1:numb_of_nodes
        Q_old{i} = consGMSF(node_data{i}, zeros(D));
    end

    lambdas = zeros(D, D, numb_of_nodes, numb_of_nodes);
%     iterNumb = 400;
    iterNumb = 50;
            traceCOEFFS = zeros(D, D, numb_of_nodes);

    for it = 1:iterNumb

        
        for i = 1:numb_of_nodes
            for j = 1:numb_of_nodes
                if(nb_mat(i, j) ~= 0)
                    traceCOEFFS(:, :, i) = traceCOEFFS(:, :, i) + step_size * (Q_old{i} - Q_old{j});
                end
            end
        end


        Q_new = cell(numb_of_nodes, 1);

        for i = 1:numb_of_nodes
            Q_new{i} = consGMSF(node_data{i}, traceCOEFFS(:, :, i));
        end
        Q_old = Q_new;
    end
    
    est_subspaces = cell(numb_of_nodes, 1);
    for i = 1:numb_of_nodes
    	[eVec, ~] = svd(Q_new{i});
        est_subspaces{i}  = eVec(:, (D - d + 1):D);
    end
end