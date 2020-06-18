function [ weight] = get_weight(elem,NODE,VERT,VERT_ELE,PARAMS)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% gauss_six = [ 0.445948490915965   0.445948490915965   
%               0.445948490915965   0.108103018168070   
%               0.108103018168070   0.445948490915965   
%               0.091576213509771   0.091576213509771  
%               0.091576213509771   0.816847572980459   
%               0.816847572980459   0.091576213509771];
        
w_six = [0.111690794839006 0.111690794839006 0.111690794839006 0.054975871827661 0.054975871827661 0.054975871827661]; 

nlink = PARAMS.nlink;
ndof = PARAMS.ndof;

weight = zeros(6,1);

%if(~isempty(elem.vertices_ele))
    
    if(~isempty(elem.vertices))
        if(PARAMS.refinement)
            Id = find(elem.vertices_ele_virt(1:3)==1);
        else
            Id = find(elem.nodes(1:3)<=PARAMS.num_bg_nod);
        end
        
        if(length(Id)==1)          
            xe = zeros(3,1);
            ye = zeros(3,1);
            
            if(PARAMS.refinement)            
                xe(1) = [VERT_ELE(elem.vertices_ele(Id)).X];
                ye(1) = [VERT_ELE(elem.vertices_ele(Id)).Y];
            else
                 xe(1) = [NODE(elem.nodes(Id)).X];
                 ye(1) = [NODE(elem.nodes(Id)).Y];
            end
           
            xe(2:3) = [VERT(elem.vertices).X];
            ye(2:3) = [VERT(elem.vertices).Y];
            
            [gauss_x,gauss_y] = gauss_map(elem,NODE);
        
            weight = moments_fitting(xe,ye,gauss_x,gauss_y); 
            
            if(sum(weight)<0)
              weight = -weight;
            end
            
        else
             if(PARAMS.refinement)         
                 [quadCoord] = orderWedgeConn( VERT_ELE( elem.vertices_ele(Id) ), VERT(elem.vertices) );
             else
                 [quadCoord] = orderWedgeConn(NODE(elem.nodes(Id)),VERT(elem.vertices));
             end

        
            xe = quadCoord(1,:)';
            ye = quadCoord(2,:)';
           
            [gauss_x,gauss_y] = gauss_map(elem,NODE);
        
            weight = moments_fitting(xe,ye,gauss_x,gauss_y);
        end
    else
        if(PARAMS.refinement)         
            xe = [VERT_ELE(elem.vertices_ele).X];
            ye = [VERT_ELE(elem.vertices_ele).Y];
        else
            xe = [NODE(elem.nodes).X];
            ye = [NODE(elem.nodes).Y];
        end
        
        xe = xe'; ye = ye';
        
        [gauss_x,gauss_y] = gauss_map(elem,NODE);
        
        weight = moments_fitting(xe,ye,gauss_x,gauss_y);
    
    end
%else
%         xe = [NODE(elem.nodes).X];
%         ye = [NODE(elem.nodes).Y];
%         xe = xe';
%         ye = ye';
%         
%         for i = 1:6        
%             GN    =  [1 0 -1 ;
%                   0 1 -1];
%             Jacobian     = GN*[xe ye];        % Compute Jacobian matrix 
%             detJ  = det(Jacobian);
%             weight(i) = w_six(i)*detJ;
%         end
%end

end
       
function[quadCoord] = orderWedgeConn(pNod,vert)
if(isempty(pNod))
    quadCoord = [vert.X; vert.Y];
else
    quadCoord = [vert.X pNod.X; vert.Y pNod.Y];
end
%%%% Subroutine to order intersection points %%%%%%%%%%%%%%
%%%%% Prevents formation of a bow-tie configuration %%%%%%%%%%%

flag = 0;
count = 0;

while(flag==0)
    AB = [quadCoord(1,2)-quadCoord(1,1),quadCoord(2,2)-quadCoord(2,1),0];
    BC = [quadCoord(1,3)-quadCoord(1,2),quadCoord(2,3)-quadCoord(2,2),0];
    crossp1 = cross(AB,BC); 
    
    CD = [quadCoord(1,4)-quadCoord(1,3),quadCoord(2,4)-quadCoord(2,3),0];
    DA = [quadCoord(1,1)-quadCoord(1,4),quadCoord(2,1)-quadCoord(2,4),0];
    crossp2 = cross(CD,DA);
    
    dotnormals = dot(crossp1,crossp2);    
    
    if(dotnormals<0)
        count = count + 1;
        if(count==1)
            quadCoord(:,[3 4]) = quadCoord(:,[4 3]);
        elseif(count==2)
            quadCoord(:,[2 3]) = quadCoord(:,[3 2]);
        end
        
    else
        flag = 1;
        checknormal = dot(crossp1,[0 0 -1]);
        if (checknormal>0)
            quadCoord(:,[2 4]) = quadCoord(:,[4 2]);
        end
    end
end
end

