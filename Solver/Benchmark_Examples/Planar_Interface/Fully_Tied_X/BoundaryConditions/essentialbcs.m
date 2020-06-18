function [ifixu, u2] = essentialbcs(NODE,PARAMS,u1)

 tol = 1.0e-13;
 
num_dof=PARAMS.ndof*length(NODE);

ifixu = zeros(num_dof,1); u2 = u1;

% ifixu(1) = 1; ifixu(2) = 1; ifixu(4) = 1;
% for nod=1:PARAMS.num_bg_nod    
%     loc_nod = PARAMS.ndof*(nod-1)+1;
%     %%%% Top Surface %%%%%%%%%%%
%     if(abs(NODE(nod).X-PARAMS.Xmax)<tol)
%         [uex] = exactsolution(NODE(nod).X,NODE(nod).Y,NODE(nod).domain,PARAMS);
%         if(PARAMS.ndof==2)
%             u(loc_nod:loc_nod+1) = timeStepping(PARAMS)*uex;
%             ifixu(loc_nod:loc_nod+1) = [1;0];
%         else
%             u(loc_nod) = uex(2);
%             ifixu(loc_nod) = 1;
%         end
%     end
%     %%%% Bottom Surface %%%%%%%%%
%     if(abs(NODE(nod).X-PARAMS.Xmin)<tol)
%         [uex] = exactsolution(NODE(nod).X,NODE(nod).Y,NODE(nod).domain,PARAMS);
%         if(PARAMS.ndof==2)
%             u(loc_nod:loc_nod+1) = timeStepping(PARAMS)*uex;
%             ifixu(loc_nod:loc_nod+1) = [1;0];
%         else
%             u(loc_nod) = uex(2);
%             ifixu(loc_nod) = 1;
%         end
%     end
% end
% ifixu(1)=1;
end