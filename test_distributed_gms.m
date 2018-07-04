d = 3; % The dimensionality of inliers 3
D = 50; % The dimensionality of outlier 50
numb_of_inliers = 1000; % Number of inliers
numb_of_outliers = 5000; % Number of outliers
numb_of_nodes = 5; % Number of nodes in the network
min_inliers = 200;
min_outliers = 1000;
% outVariance = 20; % Variance
eps = 1e-1;
iterNumb = 5000;

%create data set;
[total_data, node_data, orig_subspace] = create_data_RSR(numb_of_inliers, numb_of_outliers, numb_of_nodes, min_inliers, min_outliers,  D, d, eps);

% create the neighborhood matrix

% nb_mat = genNeigMat(numb_of_nodes); % Network's weight matrix
% nb_mat = [1 1 1 0 0 0 1 0 0 1; 1 1 1 0 0 1 0 1 0 0; 1 1 1 1 0 0 0 0 0 1; 0 0 1 1 1 1 0 1 1 0; 0 0 0 1 1 1 0 0 0 0; 0 1 0 1 1 1 1 0 0 0; 1 0 0 0 0 1 1 1 0 0; 0 1 0 1 0 0 1 1 1 0; 0 0 0 1 0 0 0 1 1 1; 1 0 1 0 0 0 0 0 1 1];
nb_mat = ones(numb_of_nodes);

%%
% 
% tic
% est_subspaces_GMS = d_GMS( node_data, nb_mat, d);
% toc

% run the distribtued gms algorithm

tic
est_subspaces_GMS = d_gms( node_data, nb_mat, d);
toc

% compare the results at each processor with the original subspace

disp(subspace(est_subspaces_GMS{1}, orig_subspace));

disp(subspace(est_subspaces_GMS{2}, orig_subspace));

disp(subspace(est_subspaces_GMS{3}, orig_subspace));

disp(subspace(est_subspaces_GMS{4}, orig_subspace));

disp(subspace(est_subspaces_GMS{5}, orig_subspace));
