function feature_vector = generateFeatures(poreTxtFilename)

global N  lattice_const;

% Initialize time
t0=0;

% Grid size
% N=5002;
% num_dt = 10000;
% start_fitting = 5000;

pore_data = dlmread(poreTxtFilename);
N = pore_data(1,1);

% Determine number of each interaction
[C,~]=init(N,poreTxtFilename);
lattice_const = 1.4226*sqrt(3);

feature_vector = visualize(C,N,poreTxtFilename);
end


function y=pbc(x)
global N;
if (x==0 || x==N)
    y=N;
else
    y=mod(x,N);
end
end

function [nearest,next_nearest]=neighbor_sites(i,j,species)

if (species==1)
    
    nearest = [pbc(i-1),    pbc(j-1+(~mod(i,2))),        2;
               pbc(i-1),    pbc(j+(~mod(i,2))),          2;
               i,           j,                           2];
    
    next_nearest = [pbc(i-1),    pbc(j-1+(~mod(i,2))),        species;
        pbc(i-1),                pbc(j+(~mod(i,2))),          species;
        i,                       pbc(j-1),                    species;
        i,                       pbc(j+1),                    species;
        pbc(i+1),                pbc(j-1+(~mod(i,2))),        species;
        pbc(i+1),                pbc(j+(~mod(i,2))),          species];
    
elseif (species==2)

    nearest = [i,       j,                       1;
        pbc(i+1),       pbc(j+(~mod(i,2))),      1;
        pbc(i+1),       pbc(j-1+(~mod(i,2))),    1];
    
    next_nearest = [pbc(i-1),       pbc(j-1+(~mod(i,2))),        species;
        pbc(i-1),                   pbc(j+(~mod(i,2))),          species;
        i,                          pbc(j-1),                    species;
        i,                          pbc(j+1),                    species;
        pbc(i+1),                   pbc(j-1+(~mod(i,2))),        species;
        pbc(i+1),                   pbc(j+(~mod(i,2))),          species];
end
end


function y=count_neighbors(i,j,species,C)
[nearest,next_nearest]=neighbor_sites(i,j,species);
[n_near, ~] = size(nearest);
[n_nextnear, ~] = size(next_nearest);

y=[0;0];
for k=1:n_near
    posr = nearest(k,1);
    posc = nearest(k,2);
    neigh_spec = nearest(k,3);
    if (C(posr,posc,neigh_spec)==neigh_spec)
        y(1)=y(1)+1;
    end
end
for k=1:n_nextnear
    posr = next_nearest(k,1);
    posc = next_nearest(k,2);
    next_neigh_spec = next_nearest(k,3);
    if (C(posr,posc,next_neigh_spec)==next_neigh_spec)
        y(2)=y(2)+1;
    end
end
end

function [y,complete_members,incomplete_members]=is_complete_hexagon(i,j,C)
count_filled=0;
row_mod_two = ~mod(i,2);

incomplete_members=[];
complete_members  =[];
if (C(i,j,1)==0)   % 0 means vacancy not added so site is filled
    count_filled=count_filled+1;
    complete_members = [complete_members; i  j 1];
else
    incomplete_members = [incomplete_members; i  j 1];
end
if (C(i,j,2)==0)   % 0 means vacancy not added so site is filled
    count_filled=count_filled+1;
    complete_members = [complete_members; i  j 2];
else
    incomplete_members = [incomplete_members; i  j 2];
end
if (C(pbc(i-1),pbc(j+row_mod_two),2)==0)   % 0 means vacancy not added so site is filled
    count_filled=count_filled+1;
    complete_members = [complete_members; pbc(i-1),pbc(j+row_mod_two),2];
else
    incomplete_members = [incomplete_members; pbc(i-1),pbc(j+row_mod_two),2];
end
if (C(pbc(i+1),pbc(j+row_mod_two),1)==0)   % 0 means vacancy not added so site is filled
    count_filled=count_filled+1;
    complete_members = [complete_members; pbc(i+1),pbc(j+row_mod_two),1]; 
else
    incomplete_members = [incomplete_members; pbc(i+1),pbc(j+row_mod_two),1]; 
end
if (C(i,pbc(j+1),1)==0)   % 0 means vacancy not added so site is filled
    count_filled=count_filled+1;
    complete_members = [complete_members; i,pbc(j+1),1];     
else
    incomplete_members = [incomplete_members; i,pbc(j+1),1];     
end
if (C(i,pbc(j+1),2)==0)   % 0 means vacancy not added so site is filled
    count_filled=count_filled+1;
    complete_members = [complete_members; i,pbc(j+1),2];         
else
    incomplete_members = [incomplete_members; i,pbc(j+1),2];         
end

y=0;  % assuming hexagon is NOT complete
if (count_filled==6)
    y=1;  % hexagon is complete
end
end

function [C0,num_atoms]=init(N,poreTxtFilename)
% This function initializes a seed MoS2 crystal on the surface, and returns
% the surface state vector. The function currently returns a triangle
% 1 - SU
% 2 - Mo
% 3 - SL

% Initialize
C0 = zeros(N,N,2);

pore = dlmread(poreTxtFilename);
pore = pore(2:end,:);

num_atoms=[0,0,0];
for i=1:N
    for j=1:N
        for k=1:2
            [lia,~] =   ismember([k,i,j],pore,'rows');
            if (any(lia))
                C0(i,j,k) = k;
                num_atoms(k) = num_atoms(k)+1;
            end
        end
    end
end

end

function feature_vector = visualize(C,N,poreTxtFilename)
% This function plots the current state of the system in a hexagonal
% lattice, given the surface state matrix
% Inputs: Surface state matrix (C), Size of grid (N)
% Outputs: none

global lattice_const;

% Define box sizes in the x and y directions
box_x = N*2*(lattice_const/sqrt(3))*cos(pi/6);
box_y = N*(lattice_const/sqrt(3))*(2+2*cos(pi/3))/2;

[m,n]=size(C(:,:,1));
num_total_atoms = m*n*2;
num_left_atoms = sum(sum(C(:,:,1)==0))+sum(sum(C(:,:,2)==0));
num_removed_atoms = m*n*2-num_left_atoms;

rim_atoms=[];
% fourth column in rim_atoms is the type of edge site (0- not an edge site,
% 1- AC site, 2-ZZ site, 3-UA site)
num_5_membered_ring=0;
num_incomplete_hexagons=0;
for i=1:N
    for j=1:N        
            [y,complete_members,~]=is_complete_hexagon(i,j,C);
            if (y==0) % if hexagon is incomplete
                rim_atoms = [rim_atoms; complete_members];
                
                if (size(complete_members,1)==5)
                    num_5_membered_ring = num_5_membered_ring + 1;
                end
                num_incomplete_hexagons = num_incomplete_hexagons+1;
            end
   end
end

% remove duplicate entries from the rim atoms
rim_atoms = unique(rim_atoms,'rows');
num_rim_atoms = size(rim_atoms,1);

% add extra column at end to indicate kind of edge 
rim_atoms(:,end+1) = zeros(num_rim_atoms,1);

% count the number of AC, ZZ, and UA edge atoms
num_AC_atoms=0;
num_ZZ_atoms=0;
num_UA_atoms=0;
num_SB_atoms=0;
cont=0;

pos_H = [];
pos_curr_l = [];
for i=1:num_rim_atoms
    % rim atom is itself a filled site, so no vacancy, i.e. C(row_pos,col_pos,species)=0
    row_pos = rim_atoms(i,1);
    col_pos = rim_atoms(i,2);
    species = rim_atoms(i,3);
    
    [nearest,~]=neighbor_sites(row_pos,col_pos,species); 
    pos_curr = calculate_coordinates(row_pos,col_pos,species);
    pos_curr_l = [pos_curr_l; pos_curr];
    
    
    num_neigh_neigh = 0;
    for j=1:size(nearest,1)
        if (C(nearest(j,1), nearest(j,2), nearest(j,3)) == 0) % the neighboring site is occupied, so count neighbors of neighbors to assign ZZ, AC, etc.
            neigh_neigh = count_neighbors(nearest(j,1), nearest(j,2), nearest(j,3),C);
            num_neigh_neigh = num_neigh_neigh + neigh_neigh(1);
            
        else
            % neighboring site for the rim is unoccupied, so add functional group
            % find the vector connecting the current site to this empty
            % site, and scale it by scale_H for getting the H-atom location
            
            pos_curr = calculate_coordinates(row_pos,col_pos,species);
            pos_neigh = calculate_coordinates(nearest(j,1),nearest(j,2),nearest(j,3));
            
            len_CH = 1.09; % this is the ratio of C-H and C-C bond lengths
            len_CC = 1.42;          
            pos_H = [pos_H; (pos_curr*(len_CC-len_CH) + pos_neigh*len_CH)/len_CC ];
        end
    end
    
    y = count_neighbors(row_pos, col_pos, species,C);
    num_neighbors = y(1);
    
    % num_neighbors can be (a) 0, i.e. fully filled coordination sphere
    % (cannot be rim atom), (b) 1, i.e. 1 neighbor empty, (c) 2, i.e. 2
    % neighbors empty, (d) 3, i.e. no neighbors (cannot be rim atom)
    % 1 and 2 are plausible options
    
    if (num_neighbors ==1) % this means exactly one neighbor is a vacancy so it is an edge atom (not all rim atoms are edge atoms)
        % need to add functional group (termination)
        
        
        switch num_neigh_neigh
            case 0
                % this is a zigzag (ZZ) site
                num_ZZ_atoms = num_ZZ_atoms+1;
                rim_atoms(i,4) = 2;               
                
            case 1
                % this is an armchair (AC) site
                num_AC_atoms = num_AC_atoms+1;
                rim_atoms(i,4) = 1;
                
            case 2
                % this is an unassigned (UA) site
                num_UA_atoms = num_UA_atoms+1;
                rim_atoms(i,4) = 3;
        end
    elseif (num_neighbors ==2)
        num_SB_atoms = num_SB_atoms+1;
    end
    
end

dummy = unique(sort(pdist(pos_curr_l)));
max_dist = dummy(end);
min_dist = dummy(2);
equivalent_circular_dia = sqrt(4*num_incomplete_hexagons*(3*sqrt(3)*1.4226^2 / 2)/pi);
eccentricity = max_dist/equivalent_circular_dia;

num_edge_atoms = num_AC_atoms + num_ZZ_atoms + num_UA_atoms;
num_H_atoms = size(pos_H,1);
        
count_removed_atoms=0;
count_rim_atoms=0;
for i=1:N
    for j=1:N       
        if (C(i,j,1)==0)
            % atom is present, i.e. vacancy not added
            coord_SU = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2  ];
          
            % check if rim atom, and assign index
            [lia,~] = ismember([i,j,1],rim_atoms(:,1:3),'rows');
            if (any(lia))
               count_rim_atoms = count_rim_atoms+1;
               C(i,j,1) = num_removed_atoms + count_rim_atoms;
            end
        else          
            % assign index to removed atoms
%             fprintf(fid_isomerize,'%d %d %d\n',1,i,j);
            count_removed_atoms = count_removed_atoms+1;
            C(i,j,1) = count_removed_atoms;
            
            % print removed atoms with different name
             coord_SU = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2  ];           
  
            
        end
        
        if (C(i,j,2)==0)
            % atom is present, i.e. vacancy not added
            coord_Mo = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];
            
            % check if rim atom, and assign index
            [lia,~] = ismember([i,j,2],rim_atoms(:,1:3),'rows');
            if (any(lia))
               count_rim_atoms = count_rim_atoms+1;
               C(i,j,2) = num_removed_atoms + count_rim_atoms;
            end
        else
            
            % assign index to removed atoms
%            fprintf(fid_isomerize,'%d %d %d\n',2,i,j);
            count_removed_atoms = count_removed_atoms+1;
            C(i,j,2) = count_removed_atoms;
            
            % print removed atoms with different name
            coord_Mo = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];
            
        end
    end
end

rim_atoms_SU = rim_atoms(rim_atoms(:,3)==1,:);
rim_atoms_Mo = rim_atoms(rim_atoms(:,3)==2,:);

% print rim atoms with different element name
coord_rim_SU = [ ~mod(rim_atoms_SU(:,1),2)*(lattice_const/2) + (rim_atoms_SU(:,2)-1)*lattice_const, (rim_atoms_SU(:,1)-1)*sqrt(3)*lattice_const/2  ];
for i=1:size(coord_rim_SU,1)
%     fprintf(fid,['Ri ', num2str(coord_rim_SU(i,1)), '  ', num2str(coord_rim_SU(i,2)), '   0.000\n']);

    label = {'AC','ZZ','UA'};
    which_edge_atom = rim_atoms_SU(i,4);
    if (which_edge_atom)
%         fprintf(fid,[label{which_edge_atom},' ', num2str(coord_rim_SU(i,1)), '  ', num2str(coord_rim_SU(i,2)), '   0.000\n']);   
    end
end

coord_rim_Mo = [ ~mod(rim_atoms_Mo(:,1),2)*(lattice_const/2) + (rim_atoms_Mo(:,2)-1)*lattice_const, (rim_atoms_Mo(:,1)-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];
for i=1:size(coord_rim_Mo,1)
%     fprintf(fid,['Ri ', num2str(coord_rim_Mo(i,1)), '  ', num2str(coord_rim_Mo(i,2)), '   0.000\n']);

    label = {'AC','ZZ','UA'};
    which_edge_atom = rim_atoms_Mo(i,4);
    if (which_edge_atom)
%         fprintf(fid,[label{which_edge_atom},' ', num2str(coord_rim_Mo(i,1)), '  ', num2str(coord_rim_Mo(i,2)), '   0.000\n']);   
    end
end


% make adjacency matrix out of removed atoms, including bond orientations
adjmat_removed = zeros(num_removed_atoms, num_removed_atoms);
for i=1:size(C,1)
    for j=1:size(C,2)
        curr_index_1 = C(i,j,1);
        if (curr_index_1>0)  % current atom is removed atom or rim atom
            [nearest,~]=neighbor_sites(i,j,1);
            [n_near, ~] = size(nearest);
            coord_A = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2  ];

            for k=1:n_near
                posr = nearest(k,1);
                posc = nearest(k,2);
                neigh_spec = nearest(k,3);
                coord_B = [ ~mod(posr,2)*(lattice_const/2) + (posc-1)*lattice_const, (posr-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];

                direction = coord_B-coord_A; % find direction of bond
                
                if (direction(1)>box_x/2.0) % right
                    direction(1) = -(box_x-direction(1));
                end
                if (direction(1)<-box_x/2.0) % left
                    direction(1) = (box_x+direction(1));                   
                end
                if (direction(2)>box_y/2.0) % top
                    direction(2) = -(box_y-direction(2));
                end
                if (direction(2)<-box_y/2.0) % down
                    direction(2) = (box_y+direction(2));
                end                  
                
                slope_direction = atan(abs(direction(2)/direction(1)))*180/pi;
                
                if (slope_direction > 85 && slope_direction < 95)
                    if (direction(2)>0)
                        adj_direction = 1; % top
                    else
                        adj_direction = 2; % bottom
                    end
                elseif (slope_direction>25 && slope_direction<35)
                    if (direction(2)>0) % top
                        if (direction(1)<0)
                            adj_direction = 3; % top left
                        else
                            adj_direction = 4; % top right
                        end
                    else                % bottom
                        if (direction(1)<0)
                            adj_direction = 5; % bottom left
                        else
                            adj_direction = 6; % bottom right
                        end
                    end
                end
                                   
                neigh_index = C(posr,posc,neigh_spec); 
                if (neigh_index>0) 
                   if (curr_index_1 <= num_removed_atoms && neigh_index <= num_removed_atoms)  % both current and neighboring atoms are removed atom only
                        adjmat_removed(curr_index_1,neigh_index)=adj_direction;
                   end
                end
            end
        end
        
        curr_index_2 = C(i,j,2);
        if (curr_index_2>0) 
            [nearest,~]=neighbor_sites(i,j,2);
            [n_near, ~] = size(nearest);
            coord_B = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];

            for k=1:n_near
                posr = nearest(k,1);
                posc = nearest(k,2);
                neigh_spec = nearest(k,3);
                coord_A = [ ~mod(posr,2)*(lattice_const/2) + (posc-1)*lattice_const, (posr-1)*sqrt(3)*lattice_const/2  ];

                direction = coord_A-coord_B;
                
                if (direction(1)>box_x/2.0) % right
                    direction(1) = -(box_x-direction(1));
                end
                if (direction(1)<-box_x/2.0) % left
                    direction(1) = (box_x+direction(1));                   
                end
                if (direction(2)>box_y/2.0) % top
                    direction(2) = -(box_y-direction(2));
                end
                if (direction(2)<-box_y/2.0) % down
                    direction(2) = (box_y+direction(2));
                end                 
                                                 
                slope_direction = atan(abs(direction(2)/direction(1)))*180/pi;
                
                if (slope_direction > 85 && slope_direction < 95)
                    if (direction(2)>0)
                        adj_direction = 1; % top
                    else
                        adj_direction = 2; % bottom
                    end
                elseif (slope_direction>25 && slope_direction<35)
                    
                    if (direction(2)>0) % top
                        if (direction(1)<0)
                            adj_direction = 3; % top left
                        else
                            adj_direction = 4; % top right
                        end
                    else                % bottom
                        if (direction(1)<0)
                            adj_direction = 5; % bottom left
                        else
                            adj_direction = 6; % bottom right
                        end
                    end
                end
                                
                neigh_index = C(posr,posc,neigh_spec); 
                if (neigh_index>0)   % neighboring atom is removed or rim atom
                    if (curr_index_2 <= num_removed_atoms && neigh_index <= num_removed_atoms)  % both current and neighboring atoms are removed atom only
                        adjmat_removed(curr_index_2,neigh_index)=adj_direction;
                    end
                end
            end
        end
        
    end
end

shape_factor = 1/eccentricity;

feature_vector = [
    num_removed_atoms
    num_ZZ_atoms
    num_AC_atoms
    num_UA_atoms
    num_SB_atoms
    count_rim_atoms
    num_incomplete_hexagons
    num_5_membered_ring
    max_dist
    shape_factor
    which_symmetry_operation(adjmat_removed)
];

fclose('all');
end

function y=calculate_coordinates(posr, posc,species)
 global lattice_const;
 if (species==1)
    y = [ ~mod(posr,2)*(lattice_const/2) + (posc-1)*lattice_const, (posr-1)*sqrt(3)*lattice_const/2  ];
 else
    y = [ ~mod(posr,2)*(lattice_const/2) + (posc-1)*lattice_const, (posr-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];        

 end
end
