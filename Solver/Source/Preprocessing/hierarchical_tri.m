function [ node, x, y, ls, num_sub_elements] = hierarchical_tri( e, node, x, y , ls, PARAMS)

cell_p(1,1:3) = [1 2 3];
cell_node_p = [x(node(1:3,e)) y(node(1:3,e))];
index_node = size(node,2);
index_x   = size(x,1);

for p = 1:PARAMS.refine_levels
    cell_node_temp_ind = 0;
    cell_temp_ind = 0;
    cell_node_temp = [];
    cell_temp = [];
    for i = 1:size(cell_p,1)
         %lse = get_level_set(cell_node_p(cell_p(i,:),1), cell_node_p(cell_p(i,:),2),PARAMS); 
         cp1 = [cell_node_p(cell_p(i,1),1) cell_node_p(cell_p(i,1),2)];
         cp2 = [cell_node_p(cell_p(i,2),1) cell_node_p(cell_p(i,2),2)];
         cp3 = [cell_node_p(cell_p(i,3),1) cell_node_p(cell_p(i,3),2)];
         cps = [cp1;cp2;cp3;0.5*cp1+0.5*cp2;0.5*cp2+0.5*cp3;0.5*cp3+0.5*cp1];
         lse = get_level_set(cps(:,1),cps(:,2),PARAMS); 
        if (max(lse)*min(lse)<0.0)
            pt1_temp = cell_node_p(cell_p(i,1),:);
            pt2_temp = cell_node_p(cell_p(i,2),:);
            pt3_temp = cell_node_p(cell_p(i,3),:);
            
            cell_node_temp(cell_node_temp_ind+1,1:2) = pt1_temp;
            cell_node_temp(cell_node_temp_ind+2,1:2) = pt2_temp;
            cell_node_temp(cell_node_temp_ind+3,1:2) = pt3_temp;
            cell_node_temp(cell_node_temp_ind+4,1:2) = 0.5*(pt2_temp+pt3_temp);
            cell_node_temp(cell_node_temp_ind+5,1:2) = 0.5*(pt1_temp+pt3_temp);
            cell_node_temp(cell_node_temp_ind+6,1:2) = 0.5*(pt1_temp+pt2_temp);
            
            cell_temp(cell_temp_ind+1,1:3) = [cell_node_temp_ind+1 cell_node_temp_ind+6 cell_node_temp_ind+5];
            cell_temp(cell_temp_ind+2,1:3) = [cell_node_temp_ind+2 cell_node_temp_ind+4 cell_node_temp_ind+6];
            cell_temp(cell_temp_ind+3,1:3) = [cell_node_temp_ind+4 cell_node_temp_ind+5 cell_node_temp_ind+6];
            cell_temp(cell_temp_ind+4,1:3) = [cell_node_temp_ind+3 cell_node_temp_ind+5 cell_node_temp_ind+4];
            
            cell_node_temp_ind = cell_node_temp_ind+6;
            cell_temp_ind = cell_temp_ind+4;
        else
            pt1_temp = cell_node_p(cell_p(i,1),:);
            pt2_temp = cell_node_p(cell_p(i,2),:);
            pt3_temp = cell_node_p(cell_p(i,3),:);
            
            cell_node_temp(cell_node_temp_ind+1,1:2) = pt1_temp;
            cell_node_temp(cell_node_temp_ind+2,1:2) = pt2_temp;
            cell_node_temp(cell_node_temp_ind+3,1:2) = pt3_temp;
            
            cell_temp(cell_temp_ind+1,1:3) = [cell_node_temp_ind+1 cell_node_temp_ind+2 cell_node_temp_ind+3];
            
            cell_node_temp_ind = cell_node_temp_ind+3;
            cell_temp_ind = cell_temp_ind+1;
        end
    end
    
    cell_p = cell_temp;
    cell_node_p = cell_node_temp;
end

if(PARAMS.nlink==6)
   cell_p = [cell_p cell_p];
end

node = [node (cell_p+index_x)'];
x = [x;cell_node_p(:,1)];
y = [y;cell_node_p(:,2)];


% update level set function
%ls = [ls; cell_node_p(:,1) - (PARAMS.xint*(PARAMS.Xmax-PARAMS.Xmin)+PARAMS.Xmin)];
%ls = [ls; cell_node_p(:,1).^2+cell_node_p(:,2).^2-1/4];
ls = [ls; get_level_set(cell_node_p(:,1),cell_node_p(:,2),PARAMS)];


num_sub_elements = size(cell_p,1);

% hold on
% plot_mesh(cell_node_p,cell_p,'T3','g.-',1.2);

end

