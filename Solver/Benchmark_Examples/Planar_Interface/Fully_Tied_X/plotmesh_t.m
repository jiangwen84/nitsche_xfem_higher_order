%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heat conduction in 2D                     %
% Haim Waisman, Columbia University         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotmesh_t;
include_flags;
 
figure(1);
 
if strcmpi(plot_mesh,'yes')==1;  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % plot Natural BC-edge ED
%     for i=1:nbe1
%         node1 = n_bc1(1,i);        % first node
%         node2 = n_bc1(2,i);        % second node
%         x1 = x(node1); y1=y(node1);    % coordinate of the first node
%         x2 = x(node2); y2=y(node2);    % coordinate of the second node
%         plot([x1 x2],[y1 y2],'g','LineWidth',3); 
%         hold on
%     end
%    
%     xlabel('X (m)');
%     ylabel('y (m)');
%     title('Mesh Plot with Natural B.C.-edge');
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % plot Natural BC-edge EF
%     for i=1:nbe2
%         node1 = n_bc2(1,i);        % first node
%         node2 = n_bc2(2,i);        % second node
%         x1 = x(node1); y1=y(node1);    % coordinate of the first node
%         x2 = x(node2); y2=y(node2);    % coordinate of the second node
%         plot([x1 x2],[y1 y2],'r','LineWidth',3); 
%         hold on
%     end
%     
%        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % plot Natural BC-edge GF
%     for i=1:nbe3
%         node1 = n_bc3(1,i);        % first node
%         node2 = n_bc3(2,i);        % second node
%         x1 = x(node1); y1=y(node1);    % coordinate of the first node
%         x2 = x(node2); y2=y(node2);    % coordinate of the second node
%         plot([x1 x2],[y1 y2],'y','LineWidth',3); 
%         hold on
%     end
    
    %legend('Natural B.C. (Flux on Line ED)','Natural B.C. (Flux on Line EF)','Natural B.C. (Flux on Line FG)');
 
    % plot Mesh
    for i = 1:nel
        x1 = [x(IEN(1,i)) x(IEN(2,i)) x(IEN(3,i)) x(IEN(1,i))];
        y1 = [y(IEN(1,i)) y(IEN(2,i)) y(IEN(3,i)) y(IEN(1,i))];
        plot(x1,y1,'g');hold on;
 
%         if strcmpi(plot_nod,'yes')==1 
%             text(x1(1),y1(1),sprintf('%0.5g',IEN(1,i)));
%             text(x1(2),y1(2),sprintf('%0.5g',IEN(2,i)));
%             text(x1(3),y1(3),sprintf('%0.5g',IEN(3,i)));
%             text(x1(4),y1(4),sprintf('%0.5g',IEN(4,i)));
%         end
        
    end
end
 
fprintf(1,'  Mesh Params \n');
fprintf(1,'No. of Elements  %d \n',nel);
fprintf(1,'No. of Nodes     %d \n',nnp);
fprintf(1,'No. of Equations %d \n\n',neq);
