function adjnew=add_weighing_nodes_in_between(adj)
% To represent directionalities within the graph isomorphism framework, 
% different edge lengths are assigned to differently-oriented bonds. 
% Accordingly, a pre-determined number of fictitious atoms are introduced into 
% the antimolecule and are inserted in between any two neighboring carbon atoms,
% as follows: (i) for a top-oriented bond, zero fictitious atoms are introduced,
% (ii) for a bottom-left-oriented bond, one fictitious atom is introduced, 
% and (iii) for a bottom-right-oriented bond, two fictitious atoms are introduced.

% Input: Adjacency matrix without weighing nodes added
% Output: Adjacency matrix with weighing nodes added in between nodes

% Obtain number of atoms in antimolecule, N
N=size(adj,1);
currN=N;
adjnew=adj;

for l=1:N
    for m=1:N
        % For a top oriented bond, do nothing
        
        % For a bottom-left-oriented bond, add 1 node in between. Just to
        % differentiate, the connections between the fictitious nodes and the 
        % original nodes are represented using "50", instead of "5"
        if (adjnew(l,m)==5) 
            adjnew(currN+1,:) = zeros(1,currN);
            adjnew(:,currN+1) = zeros(currN+1,1);
            adjnew(l,currN+1) = 50;
            adjnew(currN+1,m) = 50;
            adjnew(l,m)=0;
            
            currN=currN+1; % total number of nodes increases by 1
            
        % For a bottom-right-oriented bond, add 2 nodes in between. Just to
        % differentiate, the connections between the fictitious nodes and
        % the original nodes are represented using "60", instead of "6"
        elseif (adjnew(l,m)==6)
            adjnew(currN+1,:) = zeros(1,currN);
            adjnew(currN+2,:) = zeros(1,currN);
            adjnew(:,currN+1) = zeros(currN+2,1);
            adjnew(:,currN+2) = zeros(currN+2,1);
            
            adjnew(l,currN+1) = 60;
            adjnew(currN+1,currN+2) = 60;
            adjnew(currN+2,m) = 60;
            adjnew(l,m)=0;
            
            currN=currN+2; % total number of nodes increases by 2
        end
    end
end
end