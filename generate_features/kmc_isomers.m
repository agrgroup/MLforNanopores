function [tf,tknock_timeseries,num_dangling_bonds_timeseries,num_dangling_bonds_CH_timeseries,num_dangling_bonds_CH2_timeseries,num_AC,num_ZZ,num_UA,num_5R]=kmc_isomers(isoindex,poresize,basedir)
% This function performs a single kinetic Monte Carlo simulation for the
% formation of a nanopore in a 2D hexagonal lattice and returns the final
% configuration of the nanopore, in the absence of edge diffusion effects

% Inputs: isoindex - index of the isomer, poresize - size of the nanopore
% in terms of the number of atoms removed, basedir - directory in which to
% store the generated files

% Outputs: tf - formation time of nanopore, tknock_timeseries - time
% required for atom to be knocked out at each timestep,
% num_dangling_bonds_timeseries - number of dangling bonds as a function of
% time, num_dangling_bonds_CH_timeseries - number of doubly bonded carbons, 
% num_dangling_bonds_CH2_timeseries - number of singly bonded carbons, 
% num_AC - number of armchair edge atoms, num_ZZ - number of zigzag edge atoms,
% num_UA - number of unassigned edge atoms (in this study, etching and diffusion
% barriers of UA = that of AC),num_5R - number of 5-membered rings 

global basedirectory nuperp N T R NA kb kb_eV area_site lattice_const num_dt isomer_index pore_size t_desired;

% Initialize time
t0=0;
t_desired=40;

% Set isomer index, pore size, and base directory global variables
isomer_index = isoindex;
pore_size = poresize;
basedirectory = basedir;

% Set grid size
N = 10;

% Number of KMC steps required to form pore of given size
num_dt = pore_size-2;

% Values of physical constants
R=8.314;        % Gas constant in J/mol-K
NA=6.022e23;    % Avogadro's number in mol^-1
kb=R/NA;        % Boltzmann constant in SI units
kb_eV=8.6173324e-5;% Boltzmann constant in eV/K

% Set attempt frequency for A and B sublattices
nuperp = 1e13*[1,1];

% Set temperature of system in Kelvin
T = 1000+273.15;

% Set lattice size of graphene
lattice_const = 1.422591*sqrt(3);
area_site = 0.50*sqrt(3)*(lattice_const^2)*1e-20;
      
% Initialize 2D material layer with a monoatom vacancy
[~,C,num_atoms]=init(N);

% Do the KMC algorithm
tic
[tf,tknock_timeseries,num_dangling_bonds_timeseries,num_dangling_bonds_CH_timeseries,num_dangling_bonds_CH2_timeseries,num_AC,num_ZZ,num_UA,num_5R]=kmc(C,num_atoms,t0,N);
toc
    
end

function [t,tknock_timeseries,num_dangling_bonds_timeseries,num_dangling_bonds_CH_timeseries,num_dangling_bonds_CH2_timeseries,num_AC,num_ZZ,num_UA,num_5R]=kmc(C0,num_atoms,t0,N)
% This function performs the KMC algorithm and returns final configuration
% of the system

% Inputs: initial surface configuration (C0), number of atoms on A and B sublattices (num_atoms), initial time (t0), size of
% grid (N)

% Outputs: tf - formation time of nanopore, tknock_timeseries - time
% required for atom to be knocked out at each timestep,
% num_dangling_bonds_timeseries - number of dangling bonds as a function of
% time, num_dangling_bonds_CH_timeseries - number of doubly bonded carbons, 
% num_dangling_bonds_CH2_timeseries - number of singly bonded carbons, 
% num_AC - number of armchair edge atoms, num_ZZ - number of zigzag edge atoms,
% num_UA - number of unassigned edge atoms (in this study, etching and diffusion
% barriers of UA atoms = that of AC atoms),num_5R - number of 5-membered rings 

global num_dt t_desired;

t=t0; % initial time
count=0; % counter for number of time steps
C=C0; % store initial state of system in variable C

% Determine the list of sites which can be etched
sites_list=sites(C,N);
r = sites_list(:,7); % vector of all rates
num=length(r);

% Initialize variables to store timeseries
num_dangling_bonds_timeseries =  zeros(num_dt+1,1);
num_dangling_bonds_CH_timeseries =  zeros(num_dt+1,1);
num_dangling_bonds_CH2_timeseries =  zeros(num_dt+1,1);
tknock_timeseries = zeros(num_dt+1,1);

while (num>0 && count<=num_dt)
        
    % Generate two uniformly-distributed random numbers u1 and u2
    random = rand(1,2);
    u1=random(1);
    u2=random(2);
      
    % Compute effective rate vector
    [num,~]=size(sites_list);
    rcum = cumsum(r);
    rcum=rcum/rcum(num); % normalize with total rate
    
    % Store number of dangling bonds vs time
    [data1,data2,data3] = count_dangling_bonds(sites_list);
    num_dangling_bonds_timeseries(count+1) = data1;
    num_dangling_bonds_CH_timeseries(count+1) = data2;
    num_dangling_bonds_CH2_timeseries(count+1) = data3;
    tknock_timeseries(count+1) = 1/sum(r);

    % Find which of the sites is going to have a reaction, based on
    % Gillespie's algorithm
    site_sel=0;
    for i=1:num       
        if (i==1 && u1<=rcum(i))
            site_sel=i;
        elseif (i>1 && u1>rcum(i-1) && u1<=rcum(i))
            site_sel=i;
        end     
    end
    
    % Increment the system time based on a Poisson process
    rnet=sum(r);
    if (rnet>0)
        dt = 1/rnet;
%         dt = log(1/u2)/rnet;
        t=t+dt;   
        count=count+1;
    end

    if (site_sel>0)
        % Find (r,c) position of the lattice site selected to be etched
        pos_sel_r = sites_list(site_sel,4); 
        pos_sel_c = sites_list(site_sel,5);
        
        % Find sublattice at which etching is going to occurr
        sublattice_sel = sites_list(site_sel,1);
        
        % Find if selected lattice site is already etched
        site_status = sites_list(site_sel,6);
        
        if (site_status==0)
           % If selected lattice site is 'not etched', then update site
           % status to 'etched'
           C(pos_sel_r,pos_sel_c,sublattice_sel) = sublattice_sel;
           num_atoms(sublattice_sel) = num_atoms(sublattice_sel)+1;         
        end           
        
        % Update list of available sites
        [sites_list, ~]=sites_change(C,sublattice_sel,[pos_sel_r,pos_sel_c],sites_list);
        
        % Compute the next step's rate vector
        r = sites_list(:,7);                    
    end
      
          
    if (mod(count,1000)==0)
        disp(['count = ', num2str(count)]);
    end
end

% Extract and save the number of dangling bonds, and the knocking time at
% each time step
[data1,data2,data3] = count_dangling_bonds(sites_list);
num_dangling_bonds_timeseries(end+1) = data1;
num_dangling_bonds_CH_timeseries(end+1) = data2;
num_dangling_bonds_CH2_timeseries(end+1) = data3;
tknock_timeseries(end+1) = 1/sum(r);

% Save XYZ file and antimolecule adjacency matrix at the last timestep
[num_AC,num_ZZ,num_UA,num_5R]=visualize(C,N,count+1);       
end

function [num_dangling_bonds,num_dangling_bonds_CH,num_dangling_bonds_CH2]=count_dangling_bonds(sites_list)
% This function determines the number of total, singly-bonded, and doubly-bonded
% dangling bonds in the system

% Input: list of atomic sites in the system
% Output: total number of dangling bonds (num_dangling_bonds),
% doubly-bonded dangling bond (num_dangling_bonds_CH), 
% singly-bonded dangling bond (num_dangling_bonds_CH2)

% Initialize the count for all dangling bonds to zero
num_dangling_bonds = 0;
num_dangling_bonds_CH = 0;
num_dangling_bonds_CH2 = 0;

% Loop through all the atomic sites
for i=1:size(sites_list,1)
    if (sites_list(i,6)==1) % If the site is a vacancy
        if (sites_list(i,2)>0) % If number of nearest neighboring vacancies is greater than zero
             num_dangling_bonds = num_dangling_bonds + (3-sites_list(i,2));
        end      
    else                    % If the site has an atom
        if (sites_list(i,2)==2)     % 2 vacancy neighbors, so a 'CH2' type dangling bond
            num_dangling_bonds_CH2 = num_dangling_bonds_CH2 + 1; 
        elseif (sites_list(i,2)==1) % 1 vacancy neighbor,  so a 'CH'  type dangling bond
            num_dangling_bonds_CH = num_dangling_bonds_CH + 1;
        end
    end
end
end

function y=atomtype(i,j,sublattice,C)
% This function classifies a given lattice site to be an armchair, zigzag,
% singly-bonded, lone, or fully-bonded carbon atom

% Inputs: location (i,j) of lattice site, sublattice, state of the system (C)
% Output: atomtype (AC, ZZ, SB, LC, FC)

% Find the list of nearest neighbor lattice sites to the given lattice site
[nearest,~]=neighbor_sites(i,j,sublattice);

% Initialize number of vacancy neighbors, number of filled neighbors, and
% number of filled neighbors of filled neighbors, all to be zero
num_vacancy_neigh = 0;
num_atomic_neigh  = 0;
num_atomic_neigh_atomic_neigh = zeros(2,1);

% Loop through all nearest neighbor lattice sites
for i=1:size(nearest,1)
    if (C(nearest(i,1),nearest(i,2),nearest(i,3))==0) % nearest site is not a vacancy, so it is filled with an atom
        num_atomic_neigh = num_atomic_neigh+1;
        how_many_vacancy_neigh_atomic_neigh = count_vacancy_neighbors(nearest(i,1),nearest(i,2),nearest(i,3),C); % this function counts the number of vacancy neighbors
        num_atomic_neigh_atomic_neigh(num_atomic_neigh) = 3-(how_many_vacancy_neigh_atomic_neigh(1));
    else
        num_vacancy_neigh = num_vacancy_neigh+1;
    end
end

if (num_vacancy_neigh==0) % 0 filled vacancy neighbors, i.e. 3 filled atom neighbors
    y = 'FC';  % fully coordinated
elseif (num_vacancy_neigh==1) % 1 filled vacancy neighbors, i.e. 2 filled atomic neighbors - armchair or zigzag
    if ((num_atomic_neigh_atomic_neigh(1)==2 && num_atomic_neigh_atomic_neigh(2)==3) || (num_atomic_neigh_atomic_neigh(1)==3 && num_atomic_neigh_atomic_neigh(2)==2))     %  armchair
        y = 'AC';
    elseif (num_atomic_neigh_atomic_neigh(1)==3 && num_atomic_neigh_atomic_neigh(2)==3)  % zigzag
        y = 'ZZ';
    else % unassigned edge, treated similar to armchair in our code
        y = 'UA';
    end    
elseif (num_vacancy_neigh==2) % 1 filled atomic neighbor - singly-coordinated atom
    y = 'KL';
elseif (num_vacancy_neigh==3) % 0 filled atomic neighbors, is possible sometimes (e.g., zigzag edge etching)
    y = 'LC'; % lone carbon
end
end

function y=barrier(atomtype)
% This function returns the etching barrier for a given atomtype
% Input: atomtype
% Output: etching barrier (in eV)

switch atomtype
    case 'AC'
        y = 2.28986495;
    case 'ZZ'
        y = 2.3010539487;
    case 'UA'
        y = 2.28986495;
    case 'KL'
        y = 1.02986589;
    case 'FC'
        y = inf;
    case 'LC'
        y = 0;
end
end

function sites_list=sites(C,N)
% Compute the positions of all possible relavent interactions
% The columns are 1. sublattice, 2. number of vacancy neighbors, 
% 3. number of vacancy next neighbors, 4. row-position (y),
% 5. column-position (x), 6. vacancy flag (1-occupied by vacancy,
% 0-occupied by atom), 7. etching rate

global nuperp T kb_eV;

% Initialize monoatom vacancy in lattice, if need other shape, change
% init() function
[initial_vacancy,~,~]=init(N);

[num,~]=size(initial_vacancy);
sites_list=zeros(num,7);
sites_list(:,[1,4,5])=initial_vacancy(:,[1,2,3]); % store sublattice, row location, column location in 1st, 4th, and 5th columns of sites_list

for i=1:num 
    % This has been made a loop so that if one wants to initialize with a
    % defect other than a monoatom vacancy, the code can be easily modified simply by changing `initial_vacancy`
    sublattice = initial_vacancy(i,1);
    pos = initial_vacancy(i,2:3);

    [nearest,next_nearest]=neighbor_sites(pos(1),pos(2),sublattice);
    [n_near, ~] = size(nearest);
    [n_nextnear, ~] = size(next_nearest);

    % Update vacancy flag (i.e., 6th column) to be 1
    [~,locb] =   ismember([sublattice,pos(1),pos(2)],sites_list(:,[1,4,5]),'rows');
    sites_list(locb,6) = 1;
    
    % Loop through nearest neighbor lattice sites
    for k=1:n_near
        posr = nearest(k,1);
        posc = nearest(k,2);
        neigh_sublattice = nearest(k,3);
        
        [lia,locb] =   ismember([neigh_sublattice,posr,posc],sites_list(:,[1,4,5]),'rows');      
        if (~any(lia)) % If lattice site is not in sites_list, add it there
            sites_list = [sites_list; neigh_sublattice,1,0,posr,posc,ismember([neigh_sublattice,posr,posc],initial_vacancy,'rows'),0];
        
        else % If lattice site is already there, update number of vacancy neighbors
            sites_list(locb,2) = sites_list(locb,2)+1;
        end
    end
    
    for k=1:n_nextnear
        posr = next_nearest(k,1);
        posc = next_nearest(k,2);
        neigh_sublattice = next_nearest(k,3);
            
        [lia,locb] =   ismember([neigh_sublattice,posr,posc],sites_list(:,[1,4,5]),'rows');
        if (~any(lia)) % If lattice site is not in sites_list, add it there
            sites_list = [sites_list; neigh_sublattice,0,1,posr,posc,ismember([neigh_sublattice,posr,posc],initial_vacancy,'rows'),0];
        else % If lattice site is already there, update number of vacancy next neighbors
            sites_list(locb,3) = sites_list(locb,3)+1;
        end
    end
end   

% Loop through all lattice sites in the sites list
for i=1:size(sites_list,1)
    sublattice = sites_list(i,1);
    num_vacancy_neigh = sites_list(i,2);
    vacancy_flag = sites_list(i,6);
    posr = sites_list(i,4);
    posc = sites_list(i,5);
    
    if (vacancy_flag==0) % site doesn't have vacancy
        if (num_vacancy_neigh>0) % site has at least one vacancy neighbor
            Ea = barrier(atomtype(posr,posc,sublattice,C));
            sites_list(i,7)=nuperp(sublattice)*exp(-Ea/(kb_eV*T));  % Rate constant for etching of atom
        else
            sites_list(i,7)=0; % site has no vacancy neighbor, so cannot be etched
        end
    else  % site already has vacancy, so etching rate is zero for lattice site
        sites_list(i,7) = 0;
    end
end
end

function [updated_sites_list,modified_sites_index]=sites_change(C,sublattice,pos,sites_list)
% This function updates the sites_list for the simulation and returns an
% updated list of potential sites for etching

% Inputs: state of the system (C), properties of currently etched site
% (sublattice, pos), current sites list (sites_list)
% Outputs: updated sites list and list of modified lattice sites

global nuperp T kb_eV;

updated_sites_list = sites_list;

% Update properties of etched lattice site
[~,locb] =   ismember([sublattice,pos(1),pos(2)],updated_sites_list(:,[1,4,5]),'rows'); % find etched site in sites_list
original_state = updated_sites_list(locb,6);
updated_sites_list(locb,6) = ~original_state;
modified_sites_index = locb;

[nearest,next_nearest]=neighbor_sites(pos(1),pos(2),sublattice);
[n_near, ~] = size(nearest);
[n_nextnear, ~] = size(next_nearest);

found_neighbor = 0;

% Loop over nearest neighbors of the etched site
for k=1:n_near
    posr = nearest(k,1);
    posc = nearest(k,2);
    neigh_sublattice = nearest(k,3);

    % Try to find neighoring site in updated_sites_list
    [lia,locb] =   ismember([neigh_sublattice,posr,posc],updated_sites_list(:,[1,4,5]),'rows');
    
    if (~any(lia)) % if not found, add it to the list
        updated_sites_list = [updated_sites_list; neigh_sublattice,1,0,posr,posc,C(posr,posc,neigh_sublattice)==neigh_sublattice,0];
        [num_sites,~] = size(updated_sites_list);
        modified_sites_index = [modified_sites_index; num_sites];
    
    else % if found, update the number of vacancy neighbors to the site
        updated_sites_list(locb,2) = updated_sites_list(locb,2)+1;
        modified_sites_index = [modified_sites_index; locb];
    end
end

% Loop over next nearest neighbors of the etched site
for k=1:n_nextnear
    posr = next_nearest(k,1);
    posc = next_nearest(k,2);
    neigh_sublattice = next_nearest(k,3);

    % Try to find next neighboring site in updated_sites_list
    [lia,locb] =   ismember([neigh_sublattice,posr,posc],updated_sites_list(:,[1,4,5]),'rows');
    
    if (~any(lia)) % if not found, add it to the list
        updated_sites_list = [updated_sites_list; neigh_sublattice,0,1,posr,posc,C(posr,posc,neigh_sublattice)==neigh_sublattice,0];
        [num_sites,~] = size(updated_sites_list);
        modified_sites_index = [modified_sites_index; num_sites];
    
    else % if found, update the number of next nearest vacancy neighbors to the site
        updated_sites_list(locb,3) = updated_sites_list(locb,3)+1;
        modified_sites_index = [modified_sites_index; locb];
    end
end

% Loop over updated_sites_list to perform housekeeping check
for i=1:size(updated_sites_list,1)
    sublattice = updated_sites_list(i,1);
    num_neigh = updated_sites_list(i,2);
    num_next_neigh = updated_sites_list(i,3);
    vacancy_flag = updated_sites_list(i,6);
    posr = updated_sites_list(i,4);
    posc = updated_sites_list(i,5);
    
    if (~all(count_vacancy_neighbors(posr,posc,sublattice,C)==[num_neigh;num_next_neigh]))
        display(posr)
        display(posc)
        display(vacancy_flag)
        stored_neighs=[num_neigh,num_next_neigh];
        actual_neighs=count_vacancy_neighbors(posr,posc,sublattice,C)';
        display(sublattice)
        display(stored_neighs)
        display(actual_neighs)
        error('Neighbor coordination number does not match');
    end
end

% Loop over modified sites to assign them rates, and to remove those modified sites without nearest
% neighbor or next-nearest neighbor vacancy sites
for_removal=[];
for i=1:length(modified_sites_index)
    index = modified_sites_index(i);
    sublattice = updated_sites_list(index,1);
    num_neigh = updated_sites_list(index,2);
    num_next_neigh = updated_sites_list(index,3);
    vacancy_flag = updated_sites_list(index,6);
    pos_r = updated_sites_list(index,4);
    pos_c = updated_sites_list(index,5);
    
    if (vacancy_flag==0) % site doesn't have vacancy yet
        if (num_neigh==0 && num_next_neigh==0)
            for_removal = [for_removal;index];
        else
            if (num_neigh>0) % needs to have atleast one vacancy as neighbor
                Ea = barrier(atomtype(pos_r,pos_c,sublattice,C));
                updated_sites_list(index,7)=nuperp(sublattice)*exp(-Ea/(kb_eV*T));  % Rate for etching
            else
                updated_sites_list(index,7)=0;
            end
        end
    else % site already has vacancy, so etching rate is set to zero
            updated_sites_list(index,7) = 0;
    end
end
updated_sites_list(for_removal,:) = [];


if (found_neighbor==0 && original_state==1)
    disp('Exception: Site with no neighbor / next-nearest neighbor removed');
end
end

function y=pbc(x)
% This function implements periodic boundary conditions for lattice
% indices

% Input: row/column position before wrapping
% Output: row/column position after periodic boundary wrapping

global N;
if (x==0 || x==N)
    y=N;
else
    y=mod(x,N);
end
end

function [nearest,next_nearest]=neighbor_sites(i,j,sublattice)
% This function returns the nearest and next-nearest lattice sites for a
% given site (i, j, sublattice)

% Inputs: row (i), column (j), and sublattice (1 or 2)
% Outputs: list of nearest neighbor and next-nearest neighbor lattice sites
% in the format (row, column, sublattice)

if (sublattice==1)
    
    nearest = [pbc(i-1),    pbc(j-1+(~mod(i,2))),        2;
               pbc(i-1),    pbc(j+(~mod(i,2))),          2;
               i,           j,                           2];
    
    next_nearest = [pbc(i-1),    pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i-1),                pbc(j+(~mod(i,2))),          sublattice;
        i,                       pbc(j-1),                    sublattice;
        i,                       pbc(j+1),                    sublattice;
        pbc(i+1),                pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i+1),                pbc(j+(~mod(i,2))),          sublattice];
    
elseif (sublattice==2)

    nearest = [i,       j,                       1;
        pbc(i+1),       pbc(j+(~mod(i,2))),      1;
        pbc(i+1),       pbc(j-1+(~mod(i,2))),    1];
    
    next_nearest = [pbc(i-1),       pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i-1),                   pbc(j+(~mod(i,2))),          sublattice;
        i,                          pbc(j-1),                    sublattice;
        i,                          pbc(j+1),                    sublattice;
        pbc(i+1),                   pbc(j-1+(~mod(i,2))),        sublattice;
        pbc(i+1),                   pbc(j+(~mod(i,2))),          sublattice];
end
end

function y=count_vacancy_neighbors(i,j,sublattice,C)
% This function counts the number of vacancy neighbors to a given lattice
% site

% Input: row (i), column (j), sublattice (1 or 2), state of the system (C)
% Output: number of vacancy neighbors in the format (nearest, next-nearest)

[nearest,next_nearest]=neighbor_sites(i,j,sublattice);
[n_near, ~] = size(nearest);
[n_nextnear, ~] = size(next_nearest);

y=[0;0];
for k=1:n_near
    posr = nearest(k,1);
    posc = nearest(k,2);
    neigh_sublattice = nearest(k,3);
    if (C(posr,posc,neigh_sublattice)==neigh_sublattice) % site contains vacancy
        y(1)=y(1)+1;
    end
end
for k=1:n_nextnear
    posr = next_nearest(k,1);
    posc = next_nearest(k,2);
    next_neigh_sublattice = next_nearest(k,3);
    if (C(posr,posc,next_neigh_sublattice)==next_neigh_sublattice) % site contains vacancy
        y(2)=y(2)+1;
    end
end
end

function [y,complete_members,incomplete_members]=is_complete_hexagon(i,j,C)
% This function determines if the (i,j)th hexagon in the graphene lattice
% has all 6 possible atomic sites filled with carbon atoms. This function
% is useful to determine the rim atoms of a nanopore, because the rim atoms
% are all those atoms which belong to an "incomplete" hexagon

% Inputs: row (i), column (j), state of the system (C)
% Output: flag (y), list of filled and vacant members of the hexagon

row_mod_two = ~mod(i,2);

count_filled=0;
incomplete_members=[];
complete_members  =[];

% We will sequentially go through all 6 possible lattice sites in the
% (i,j)th hexagon and determine the filled and unfilled lattice sites

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

function [initial_vacancy,C0,num_atoms]=init(N)
% This function initializes the state of the system

% Input: size of the lattice (N). Number of lattice sites is: N*N*2, the
% factor of 2 due to there being two sublattices in graphene
% Output: initial state of the system (initial_vacancy and C0), number of
% vacancies in each sublattice

% Initialize system with monoatom vacancy, can suitably tweak as required
C0 = zeros(N,N,3);
center=round(N/2);

initial_vacancy = [      1       center  center];
C0(center,center,1) = 1;
num_atoms = [1,0];
end

function [num_AC_atoms,num_ZZ_atoms,num_UA_atoms,num_5_membered_ring]=visualize(C,N,count)
% This function saves the current state of the system, given the surface
% state matrix, as an XYZ file and a directed adjacency matrix for the
% antimolecule of the formed nanopore

% Inputs: Surface state matrix (C), Size of grid (N), and current timestep (count)
% Outputs: number of armchair, zigzag, and unassigned atoms, and 5-membered
% rings

global lattice_const isomer_index pore_size basedirectory;

% Define box sizes in the x and y directions
box_x = N*2*(lattice_const/sqrt(3))*cos(pi/6);
box_y = N*(lattice_const/sqrt(3))*(2+2*cos(pi/3))/2;

if (count<pore_size && count>=3)
    
else
    fid_xyz = fopen(      [basedirectory,'/pore',num2str(pore_size),'/pore',            num2str(isomer_index)  ,'.xyz'],'wt');
    fid_txt = fopen(      [basedirectory,'/pore',num2str(pore_size),'/pore',            num2str(isomer_index)  ,'.txt'],'wt');
    fid_adjmat_removed =  [basedirectory,'/pore',num2str(pore_size),'/adjmat_antimol_', num2str(isomer_index)  ,'.txt'];
end

% Create variables for number of total atoms, number of remaining atoms,
% and number of removed atoms
[m,n]=size(C(:,:,1));
num_total_atoms = m*n*2;
num_remaining_atoms = sum(sum(C(:,:,1)==0))+sum(sum(C(:,:,2)==0));
num_removed_atoms = m*n*2-num_remaining_atoms;

rim_atoms=[];
% rim_atoms contains four columns: (row, col, sublattice, edge_atom_type)
% fourth column in rim_atoms is the type of edge site (0- not an edge site,
% 1- AC site, 2-ZZ site, 3-UA site)
num_5_membered_ring=0;
for i=1:N
    for j=1:N        
            [y,complete_members,~]=is_complete_hexagon(i,j,C);
            if (y==0) % if hexagon is incomplete
                rim_atoms = [rim_atoms; complete_members];
                
                if (size(complete_members,1)==5)
                    num_5_membered_ring = num_5_membered_ring + 1;
                end
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

for i=1:num_rim_atoms
    % rim atom is itself a filled site, so no vacancy, i.e. C(row_pos,col_pos,species)=0
    row_pos = rim_atoms(i,1);
    col_pos = rim_atoms(i,2);
    species = rim_atoms(i,3);
    
    [nearest,~]=neighbor_sites(row_pos,col_pos,species);    
    num_neigh_neigh = 0;
    for j=1:size(nearest,1)
        if (C(nearest(j,1), nearest(j,2), nearest(j,3)) == 0) % the neighboring site is occupied, so count neighbors of neighbors to assign ZZ, AC, etc.
            neigh_neigh = count_vacancy_neighbors(nearest(j,1), nearest(j,2), nearest(j,3),C);
            num_neigh_neigh = num_neigh_neigh + neigh_neigh(1);         
        end
    end
    
    y = count_vacancy_neighbors(row_pos, col_pos, species,C);
    num_neighbors = y(1);
    
    % num_neighbors can be (a) 0, i.e. fully filled coordination sphere
    % (cannot be rim atom), (b) 1, i.e. 1 neighbor empty, (c) 2, i.e. 2
    % neighbors empty, (d) 3, i.e. no neighbors (cannot be rim atom)
    % 1 and 2 are plausible options
    
    if (num_neighbors > 0) % this means atleast one neighbor is a vacancy so it is an edge atom (not all rim atoms are edge atoms)               
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
    end
end
num_edge_atoms = num_AC_atoms + num_ZZ_atoms + num_UA_atoms;

fprintf(fid_xyz,[num2str(num_total_atoms+num_rim_atoms+num_edge_atoms),'\n # of carbons removed: ', num2str(num_removed_atoms)  ...
            ,'; # of rim atoms: ', num2str(num_rim_atoms),  '; # of edge atoms: ', num2str(num_edge_atoms),  ' \n']);      
fprintf(fid_txt,'%d\n',N);

count_removed_atoms=0;
count_rim_atoms=0;
% Cycle through all locations on the 2D lattice
for i=1:N
    for j=1:N       
        if (C(i,j,1)==0)
            % atom is present, i.e. vacancy not added
            coord_sublattice_1 = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2  ];
            fprintf(fid_xyz,['C ', num2str(coord_sublattice_1(1)), '  ', num2str(coord_sublattice_1(2)), '   0.000\n']);           
            
            % check if rim atom, and assign index
            [lia,~] = ismember([i,j,1],rim_atoms(:,1:3),'rows');
            if (any(lia))
               count_rim_atoms = count_rim_atoms+1;
               C(i,j,1) = num_removed_atoms + count_rim_atoms;
            end
        else
            % print removed atoms with different name
            coord_sublattice_1 = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2  ];           
            fprintf(fid_xyz,['Re ', num2str(coord_sublattice_1(1)), '  ', num2str(coord_sublattice_1(2)), '   0.000\n']);
          
            % assign index to removed atoms
            fprintf(fid_txt,'%d %d %d\n',1,i,j);
            count_removed_atoms = count_removed_atoms+1;
            C(i,j,1) = count_removed_atoms;
        end
        
        if (C(i,j,2)==0)
            % atom is present, i.e. vacancy not added
            coord_sublattice_2 = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];
            fprintf(fid_xyz,['C ', num2str(coord_sublattice_2(1)), '  ', num2str(coord_sublattice_2(2)), '   0.000\n']);
            
            % check if rim atom, and assign index
            [lia,~] = ismember([i,j,2],rim_atoms(:,1:3),'rows');
            if (any(lia))
               count_rim_atoms = count_rim_atoms+1;
               C(i,j,2) = num_removed_atoms + count_rim_atoms;
            end
        else
            % print removed atoms with different name
            coord_sublattice_2 = [ ~mod(i,2)*(lattice_const/2) + (j-1)*lattice_const, (i-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];
            fprintf(fid_xyz,['Re ', num2str(coord_sublattice_2(1)), '  ', num2str(coord_sublattice_2(2)), '   0.000\n']);
            
            % assign index to removed atoms
            fprintf(fid_txt,'%d %d %d\n',2,i,j);
            count_removed_atoms = count_removed_atoms+1;
            C(i,j,2) = count_removed_atoms;
        end
    end
end
rim_atoms_sublattice_1 = rim_atoms(rim_atoms(:,3)==1,:);
rim_atoms_sublattice_2 = rim_atoms(rim_atoms(:,3)==2,:);

% print rim atoms with different element name
coord_rim_sublattice_1 = [ ~mod(rim_atoms_sublattice_1(:,1),2)*(lattice_const/2) + (rim_atoms_sublattice_1(:,2)-1)*lattice_const, (rim_atoms_sublattice_1(:,1)-1)*sqrt(3)*lattice_const/2  ];
for i=1:size(coord_rim_sublattice_1,1)
    fprintf(fid_xyz,['Ri ', num2str(coord_rim_sublattice_1(i,1)), '  ', num2str(coord_rim_sublattice_1(i,2)), '   0.000\n']);

    label = {'AC','ZZ','UA'};
    which_edge_atom = rim_atoms_sublattice_1(i,4);
    if (which_edge_atom)
        fprintf(fid_xyz,[label{which_edge_atom},' ', num2str(coord_rim_sublattice_1(i,1)), '  ', num2str(coord_rim_sublattice_1(i,2)), '   0.000\n']);   
    end
end

coord_rim_sublattice_2 = [ ~mod(rim_atoms_sublattice_2(:,1),2)*(lattice_const/2) + (rim_atoms_sublattice_2(:,2)-1)*lattice_const, (rim_atoms_sublattice_2(:,1)-1)*sqrt(3)*lattice_const/2 + (lattice_const/sqrt(3)) ];
for i=1:size(coord_rim_sublattice_2,1)
    fprintf(fid_xyz,['Ri ', num2str(coord_rim_sublattice_2(i,1)), '  ', num2str(coord_rim_sublattice_2(i,2)), '   0.000\n']);

    label = {'AC','ZZ','UA'};
    which_edge_atom = rim_atoms_sublattice_2(i,4);
    if (which_edge_atom)
        fprintf(fid_xyz,[label{which_edge_atom},' ', num2str(coord_rim_sublattice_2(i,1)), '  ', num2str(coord_rim_sublattice_2(i,2)), '   0.000\n']);   
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
dlmwrite(fid_adjmat_removed, adjmat_removed);
fclose('all');
end