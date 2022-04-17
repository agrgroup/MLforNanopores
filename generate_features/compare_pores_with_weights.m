function y=compare_pores_with_weights(pore_A,pore_B)
% This function compares the given directed adjacency matrices (i.e.,
% including information on orientation of bonds) of two pores to check if
% they are identical

% Input: directed adjacency matrices of two pores (pore_A, pore_B)
% Output: 1 if pores are identical and 0 otherwise


if (sum(sum(pore_A>0))~=sum(sum(pore_B>0)))
    y=0;
else
    
    pore_A_try = pore_A;

    % The bond orientations are numbered as follows (consistent with the file
    % generate_directed_graphs_antimolec
    
    % 1 - top
    % 2 - bottom
    % 3 - top left
    % 4 - top right
    % 5 - bottom left
    % 6 - bottom right
    
    % As explained in the paper, only 3 of the 6 orientations are unique.
    pore_A_try(pore_A_try>=2 & pore_A_try<=4) = 0;
    
    y=0;
    
    for j=1:12 % 12 symmetry operations in graphene
        pore_B_transformed = pore_B;
        
        if (j==2)
            % rotate by 60 clockwise(1->4, 6->2, 5->3 and 2->5,
            % 3->1, 4->6)
            pore_B_transformed = pore_B_transformed*10;
            
            pore_B_transformed(pore_B_transformed==10)=4;
            pore_B_transformed(pore_B_transformed==60)=2;
            pore_B_transformed(pore_B_transformed==50)=3;
            
            pore_B_transformed(pore_B_transformed==20)=5;
            pore_B_transformed(pore_B_transformed==30)=1;
            pore_B_transformed(pore_B_transformed==40)=6;
            
        elseif (j==3)
            % rotate by 120 clockwise(1->6, 6->5, 5->1  and  2->3, 3->4, 4->2)
            pore_B_transformed(pore_B_transformed==1)=10;
            pore_B_transformed(pore_B_transformed==5)=1;
            pore_B_transformed(pore_B_transformed==6)=5;
            pore_B_transformed(pore_B_transformed==10)=6;
            
            pore_B_transformed(pore_B_transformed==2)=20;
            pore_B_transformed(pore_B_transformed==4)=2;
            pore_B_transformed(pore_B_transformed==3)=4;
            pore_B_transformed(pore_B_transformed==20)=3;
            
        elseif (j==4)
            % rotate by 180 clockwise(1<->2, 3<->6, 4<->5)
            pore_B_transformed = pore_B_transformed*10;
            
            pore_B_transformed(pore_B_transformed==10)=2;
            pore_B_transformed(pore_B_transformed==20)=1;
            pore_B_transformed(pore_B_transformed==30)=6;
            
            pore_B_transformed(pore_B_transformed==60)=3;
            pore_B_transformed(pore_B_transformed==40)=5;
            pore_B_transformed(pore_B_transformed==50)=4;
        elseif (j==5)
            % rotate by 240 clockwise(1->5, 5->6, 6->1  and  2->4, 4->3, 3->2)
            pore_B_transformed(pore_B_transformed==1)=10;
            pore_B_transformed(pore_B_transformed==6)=1;
            pore_B_transformed(pore_B_transformed==5)=6;
            pore_B_transformed(pore_B_transformed==10)=5;
            
            pore_B_transformed(pore_B_transformed==2)=20;
            pore_B_transformed(pore_B_transformed==3)=2;
            pore_B_transformed(pore_B_transformed==4)=3;
            pore_B_transformed(pore_B_transformed==20)=4;
            
        elseif (j==6)
            % rotate by 300 clockwise(1->3, 2->6, 3->5, 4->1,
            % 5->2, 6->4)
            pore_B_transformed = pore_B_transformed*10;
            
            pore_B_transformed(pore_B_transformed==10)=3;
            pore_B_transformed(pore_B_transformed==20)=6;
            pore_B_transformed(pore_B_transformed==30)=5;
            
            pore_B_transformed(pore_B_transformed==40)=1;
            pore_B_transformed(pore_B_transformed==50)=2;
            pore_B_transformed(pore_B_transformed==60)=4;
            
        elseif  (j==7)
            % Mirror plane: 0 degree clockwise from vertical(3 <-> 4, 5 <-> 6)
            pore_B_transformed(pore_B_transformed==3)=30;
            pore_B_transformed(pore_B_transformed==4)=3;
            pore_B_transformed(pore_B_transformed==30)=4;
            
            pore_B_transformed(pore_B_transformed==5)=50;
            pore_B_transformed(pore_B_transformed==6)=5;
            pore_B_transformed(pore_B_transformed==50)=6;
            
        elseif (j==8)
            % Mirror plane: 90 degree clockwise from vertical (1<->2, 3<->5, 4<->6)
            pore_B_transformed(pore_B_transformed==1)=10;
            pore_B_transformed(pore_B_transformed==2)=1;
            pore_B_transformed(pore_B_transformed==10)=2;
            
            pore_B_transformed(pore_B_transformed==3)=30;
            pore_B_transformed(pore_B_transformed==5)=3;
            pore_B_transformed(pore_B_transformed==30)=5;
            
            pore_B_transformed(pore_B_transformed==4)=40;
            pore_B_transformed(pore_B_transformed==6)=4;
            pore_B_transformed(pore_B_transformed==40)=6;
            
        elseif (j==9)
            % Mirror plane: 150 degree clockwise from vertical (1<->3, 2<->6, 4<->5)
            pore_B_transformed(pore_B_transformed==1)=10;
            pore_B_transformed(pore_B_transformed==3)=1;
            pore_B_transformed(pore_B_transformed==10)=3;
            
            pore_B_transformed(pore_B_transformed==2)=20;
            pore_B_transformed(pore_B_transformed==6)=2;
            pore_B_transformed(pore_B_transformed==20)=6;
            
            pore_B_transformed(pore_B_transformed==4)=40;
            pore_B_transformed(pore_B_transformed==5)=4;
            pore_B_transformed(pore_B_transformed==40)=5;
        elseif (j==10)
            % Mirror plane: 30 degree clockwise from vertical (1<->4, 2<->5, 3<->6)
            pore_B_transformed(pore_B_transformed==1)=10;
            pore_B_transformed(pore_B_transformed==4)=1;
            pore_B_transformed(pore_B_transformed==10)=4;
            
            pore_B_transformed(pore_B_transformed==2)=20;
            pore_B_transformed(pore_B_transformed==5)=2;
            pore_B_transformed(pore_B_transformed==20)=5;
            
            pore_B_transformed(pore_B_transformed==3)=30;
            pore_B_transformed(pore_B_transformed==6)=3;
            pore_B_transformed(pore_B_transformed==30)=6;
            
        elseif (j==11)
            % Mirror plane: 60 degree clockwise from vertical (1 <-> 6, 3 <-> 2)
            pore_B_transformed(pore_B_transformed==1)=10;
            pore_B_transformed(pore_B_transformed==6)=1;
            pore_B_transformed(pore_B_transformed==10)=6;
            
            pore_B_transformed(pore_B_transformed==3)=30;
            pore_B_transformed(pore_B_transformed==2)=3;
            pore_B_transformed(pore_B_transformed==30)=2;
            
        elseif (j==12)
            % Mirror plane: 120 degree clockwise from vertical (1 <-> 5, 2 <-> 4)
            
            pore_B_transformed(pore_B_transformed==1)=10;
            pore_B_transformed(pore_B_transformed==5)=1;
            pore_B_transformed(pore_B_transformed==10)=5;
            
            pore_B_transformed(pore_B_transformed==2)=20;
            pore_B_transformed(pore_B_transformed==4)=2;
            pore_B_transformed(pore_B_transformed==20)=4;
            
        end
        
        pore_B_transformed(pore_B_transformed>=2 & pore_B_transformed<=4) = 0;
        
        %         pore_A_try
        weighted_A = add_weighing_nodes_in_between(pore_A_try);
        
        %         pore_B_rotate_flip
        weighted_B = add_weighing_nodes_in_between(pore_B_transformed);
        
        
        comparison = graphisomorphism(sparse(weighted_A), sparse(weighted_B));
        
        if (comparison==1)
            y=1;
            break;
        end
        
    end
       
end
end
