function adj=add_weighing_nodes_in_between(adj)
N=size(adj,1);
[rows, cols]=find(adj==5);
for i=1:length(rows)
 % add 1 node in between
            adj(N+1,:) = zeros(1,N);
            adj(:,N+1) = zeros(N+1,1);
            adj(rows(i),N+1) = 50;
            adj(N+1,cols(i)) = 50;
            adj(rows(i),cols(i))=0;            
            N=N+1;
end

[rows, cols]=find(adj==6);
for i=1:length(rows)
% add 2 nodes in between
            adj(N+1,:) = zeros(1,N);
            adj(N+2,:) = zeros(1,N);
            adj(:,N+1) = zeros(N+2,1);
            adj(:,N+2) = zeros(N+2,1);           
            adj(rows(i),N+1) = 60;
            adj(N+1,N+2) = 60;
            adj(N+2,cols(i)) = 60;
            adj(rows(i),cols(i))=0;            
            N=N+2;
end
    
end