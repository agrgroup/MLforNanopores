% A script to generate the physical features in CSV format

N = 21;
FOLDER_TO_SAVE_FEATURES = 'features';

% Path to the folder of nanopores with a fixed vacancy value (6 in 
% this case), which has all the KMC-generated data
PATH_TO_PORE_DATA = '../catalog/without_edge_diffusion/pore' + string(N) + '/'; 

% Paths to all the TXT file in the above folder
path = []; 
for j=1:10000
    path = [path, strcat(PATH_TO_PORE_DATA, 'pore', string(j), '.txt')];
end

% Feature matrix for all the nanopores in the above folder
feat_matrix = [];
for i=1:length(path)
    feat_array = reshape(generateFeatures(path(i)), [22, 1]);
    feat_matrix = [feat_matrix, feat_array];
end

% Create a directory to save the features in the form of CSV files
mkdir(FOLDER_TO_SAVE_FEATURES)

feat_matrix = transpose(feat_matrix);
writematrix(feat_matrix, string(FOLDER_TO_SAVE_FEATURES) + '/Physical_features-' + string(N) + '.csv','WriteMode', 'append')