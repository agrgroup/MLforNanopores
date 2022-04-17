% This program repeatedly calls the KMC routine and stochastically
% generates a new isomer each time. The properties of the isomer are stored
% in a CSV file

pause

% Re-initialize the random number generator in MATLAB
rng('shuffle')

% Define a list of pore sizes for which to generate isomers
pore_size_list  = [21 22];

% Number of isomers to generate for each size
Niso=10000;

% Directory in which the isomers will be stored. Each pore size N will have
% a subdirectory called poreN, e.g., pore6 and pore8 in the current version
% of the code
basedir = 'C:\Users\Rahul Sheshanarayana\Desktop\catalog/without_edge_diffusion/';

% Cycle through all pore sizes
for j=pore_size_list
    
    % Create directory porej
    dirname = [basedir,'pore',num2str(j)];
    status = mkdir(dirname);
    mkdir([dirname,'/path']);
    
    % Obtain properties of nanopore isomer generated
    tf_list = zeros(Niso,1);
    num_dangling_bonds_list = tf_list;
    num_dangling_bonds_CH_list = tf_list;
    num_dangling_bonds_CH2_list = tf_list;
    tknock_list = tf_list;
    num_AC_list = tf_list;
    num_ZZ_list = tf_list;
    num_UA_list = tf_list;
    num_5R_list = tf_list;
    
    % If directory was not created, it already exists, so some isomers have
    % already been generated
    if (status == 0)
       data = csvread([dirname,'/Analysis.csv']);
       num_isomers_done  = data(end,1);
    else
        
    % If directory was created, no isomers have yet been generated
       num_isomers_done = 0;
    end
    
    % Generate the remaining number of isomers
    for i=(num_isomers_done+1):Niso
        pore_index = j
        iso_index = i
        
        % Run the KMC algorithm
        [tf,tknock_timeseries,num_dangling_bonds_timeseries,num_dangling_bonds_CH_timeseries,num_dangling_bonds_CH2_timeseries,num_AC,num_ZZ,num_UA,num_5R]=kmc_isomers(i,j,basedir);
        
        % Store the isomer properties in lists
        tf_list(i) = tf;
        tknock_list(i) = tknock_timeseries(end);
        num_dangling_bonds_list(i) = num_dangling_bonds_timeseries(end);
        num_dangling_bonds_CH_list(i) = num_dangling_bonds_CH_timeseries(end);
        num_dangling_bonds_CH2_list(i) = num_dangling_bonds_CH2_timeseries(end);
        num_AC_list(i) = num_AC;
        num_ZZ_list(i) = num_ZZ;
        num_UA_list(i) = num_UA;
        num_5R_list(i) = num_5R;
       
        % Write isomer properties in CSV file
        dlmwrite([dirname,'/Analysis.csv'],[i,tf_list(i),tknock_list(i),num_dangling_bonds_list(i),num_dangling_bonds_CH_list(i),num_dangling_bonds_CH2_list(i),num_AC,num_ZZ,num_UA,num_5R],'-append');
    end
end