%function analyze_isomers_directed_isomorphism(sizes)
% This function analyzes the generated isomers using graph isomorphism and
% counts the number of isomers of each type generated (for the case without
% edge diffusion considered)

% Input: List of pore sizes to be considered, e.g., [6,8]
% Output: None

%pore_size_list  = sizes;

pore_size_list=[21 22];

% Required to write MS-Excel files on Ubuntu, comment following lines if
% using on Windows 
% Initialisation of POI Libs
% Add Java POI Libs to matlab javapath
% javaaddpath('poi_library/poi-3.8-20120326.jar');
% javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
% javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
% javaaddpath('poi_library/xmlbeans-2.3.0.jar');
% javaaddpath('poi_library/dom4j-1.6.1.jar');
% javaaddpath('poi_library/stax-api-1.0.1.jar');

% Cycle through all pore sizes
for j=pore_size_list
    tic
    % Directory in which to search for pore files
    dirname = ['../catalog/without_edge_diffusion/pore',num2str(j)];
    adjmatListing = dir([dirname,'/adjmat_antimol_*.txt']);
    mkdir([dirname,'/isomers']);
    
    dataIsomer = [];
    isomNum = [];
    data = [];
    idx = [];
    disp(['Current pore size: ', num2str(j)]);
    
    dataPore = cell(size(adjmatListing,1),j+3);
    for i=1:size(adjmatListing,1)
        i
        antimolAdjmatFilename = [dirname,'/',adjmatListing(i).name];
        xyzFilename =  strrep(antimolAdjmatFilename,'adjmat_antimol_','pore');
        xyzFilename =  strrep(xyzFilename,'.txt','.xyz');
        txtFilename =  strrep(xyzFilename,'.xyz','.txt');

        % Get the index of the pore        
        indexCurrentPore = strrep(adjmatListing(i).name,'adjmat_antimol_','');
        indexCurrentPore = num2str(strrep(indexCurrentPore,'.txt',''));
        
        pores{i}.antimolAdjmat = sparse(dlmread(antimolAdjmatFilename));

        % This is the first isomer        
        if (i==1)
            isomers{1}.antimolAdjmat = pores{1}.antimolAdjmat;
            isomers{1}.numCopies = 1;
            numIsomers(j) = 1;
            
            isomerNumber = 1;
            copyNumber = 1;
        else
           foundIsomer = 0;
           k=1;
           while(k<=numIsomers(j) && ~foundIsomer)
               if (compare_pores_with_weights(pores{i}.antimolAdjmat, isomers{k}.antimolAdjmat) )
                   % isomer exists in database. rename file with isomer number
                   foundIsomer = 1;
                   isomers{k}.numCopies = isomers{k}.numCopies+1;
                   
                   isomerNumber = k;
                   copyNumber = isomers{k}.numCopies;
               end
               k=k+1;
           end
           if (foundIsomer == 0)
               % new isomer found
               isomers{numIsomers(j)+1}.antimolAdjmat = pores{i}.antimolAdjmat;
               
               isomers{numIsomers(j)+1}.numCopies = 1;
               numIsomers(j) = numIsomers(j) + 1;
               
               isomerNumber = numIsomers(j);
               copyNumber = 1;
           end
        end
        
        dataIsomer = [dataIsomer; str2double(indexCurrentPore), isomerNumber];
        isomNum = [isomNum; isomerNumber];
        idx = [idx; str2double(indexCurrentPore)];
        
        % determined if pore is new or already found. copy associated pore
        % files and give them appropriate name based on allotted isomer
        % number and copy number
        copyfile(antimolAdjmatFilename,[dirname,'/isomers/adjmat_antimol_',sprintf('%.3d',j),'_',sprintf('%.3d',isomerNumber),'_',sprintf('%.3d',copyNumber),'.txt']);
        copyfile(xyzFilename,          [dirname,'/isomers/',               sprintf('%.3d',j),'_',sprintf('%.3d',isomerNumber),'_',sprintf('%.3d',copyNumber),'.xyz']);
        copyfile(txtFilename,          [dirname,'/isomers/',               sprintf('%.3d',j),'_',sprintf('%.3d',isomerNumber),'_',sprintf('%.3d',copyNumber),'.txt']);
       
    % end looping through pores of a given size    
    end
    for i = 1:length(isomNum)
        data = [data; isomers{isomNum(i)}.numCopies/10000];
    end
    
    % write the mapping of each isomer to isomer type to a CSV file   
    csvwrite([dirname, '/isomerData_newAlgorithm_', sprintf('%.3d',j) ,'.csv'],dataIsomer);
    csvwrite([dirname, '/isomerCounts_', sprintf('%.3d',j) ,'.csv'],[idx, isomNum, data]);
% end looping through pore sizes
number_of_isomers = numIsomers(j)
toc
end

% % Write list of isomers with their respective count to CSV file
% csvwrite([dirname,'/pore',num2str(j),'_count.csv'],[pore_size_list',numIsomers']);
%end