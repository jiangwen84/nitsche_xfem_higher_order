function[Amat] = interpolationmatrix(elem,NODE,PARAMS)

nlink = PARAMS.nlink;

xe = [NODE(elem.nodes).X]; ye = [NODE(elem.nodes).Y];
% if(nlink==4)
%     ze = [NODE(elem.nodes).Z];
%     Amat = [...
%         1     1     1     1    ;...
%         xe(1) xe(2) xe(3) xe(4);...
%         ye(1) ye(2) ye(3) ye(4);...
%         ze(1) ze(2) ze(3) ze(4)];
% elseif(nlink==3)
%     Amat = [...
%         1     1     1    ;...
%         xe(1) xe(2) xe(3);...
%         ye(1) ye(2) ye(3)];
% end

%Amat = zeros(3,nlink);

Amat = [ones(1,nlink);
        xe;
        ye];
end