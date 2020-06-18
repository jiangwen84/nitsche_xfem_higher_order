function[quadNod] = getDegenQuadCoord(elem,NODE,VERT,VERT_ELE, PARAMS)

%%  Local coordinate system for a cube element:
%
%           A
%           |
%           5----------------8
%          /|               /|
%         / |              / |
%        /  |             /  |
%       /   |            /   |
%      6----+-----------7    |
%      |    |           |    |
%      |    |           |    |
%      |    |           |    |
%      |    |           |    |
%      |    1-----------+----4-->
%      |   /            |   /
%      |  /             |  /
%      | /              | /
%      |/               |/
%      2----------------3
%     /
%    V
%%

if (~PARAMS.refinement)
    quadNod(4) = struct('X',[],'Y',[]);
    if(isempty(elem.siblings))
        %%%% Global coordinates of the uncut tetrahedral element
        xe = [NODE(elem.nodes(1:3)).X];
        ye = [NODE(elem.nodes(1:3)).Y];
        %%%% Degenerated coordinates of the uncut tetrahedral element
        quadNod(1,1).X = xe(1); quadNod(1,1).Y = ye(1);
        quadNod(2,1).X = xe(2); quadNod(2,1).Y = ye(2);
        quadNod(3,1).X = xe(3); quadNod(3,1).Y = ye(3);
        quadNod(4,1).X = xe(3); quadNod(4,1).Y = ye(3);
    elseif(length(elem.siblings)==1)
        %%%% For cut elements with just two subdomains
        %%%% Index of the physical node of the cut tetrahedral element
        physNode = get_physical_nodes(elem,NODE);
        if(length(elem.vertices)==2)
            if(length(physNode)==1)
                %%%% Global coordinates of the sub-element
                xp = [NODE(physNode).X];   xv = [VERT(elem.vertices).X];
                yp = [NODE(physNode).Y];   yv = [VERT(elem.vertices).Y];
                %%%% Degenerated coordinates of the sub-triangular element
                quadNod(1,1).X = xp(1); quadNod(1,1).Y = yp(1);
                quadNod(2,1).X = xp(1); quadNod(2,1).Y = yp(1);
                quadNod(3,1).X = xv(1); quadNod(3,1).Y = yv(1);
                quadNod(4,1).X = xv(2); quadNod(4,1).Y = yv(2);
            elseif(length(physNode)==2)
                [quadCoord] = orderWedgeConn(NODE(physNode),VERT(elem.vertices));
                for iCoord=1:4
                    quadNod(iCoord,1).X = quadCoord(1,iCoord);
                    quadNod(iCoord,1).Y = quadCoord(2,iCoord);
                end
            else
                err = MException('tetCutChk:physNodeChk',...
                    'in fourthOrderGaussQuad 3 vertices and %d physical nodes',...
                    length(physNode));
                throw(err);
            end
        else
            err = MException('tetCutChk:physNodeChk',...
                'in fourthOrderGaussQuad, wrong number of vertices: %d',...
                length(elem.vertices));
            throw(err);
        end
    elseif(length(elem.siblings)>1)
        %%%% For cut elements with more than two subdomains
        %%%% Index of the physical node of the cut tetrahedral element
        physNode = get_physical_nodes(elem,NODE);
        locVert = unique([elem.vertices(1,:) elem.vertices(2,:)]);
        if(length(locVert)==3)
            if(~isempty(physNode))
                if(length(physNode)==1)
                    [quadCoord] = orderWedgeConn(NODE(physNode),VERT(locVert));
                    for iCoord=1:4
                        quadNod(iCoord,1).X = quadCoord(1,iCoord);
                        quadNod(iCoord,1).Y = quadCoord(2,iCoord);
                    end
                else
                    % Locally three vertices and two physical nodes
                    %   --> sub-element is neither a triangle nor a quad
                    %   --> Break the sub-element into a triangle and a quad
                    % There are many ways of doing this. The strategy followed here
                    % is (a) The tri sub-element contains all three vertices and
                    % (b) the quad contains both physical nodes and the two
                    % vertices which lie on the element edge
                    
                    % Search for vertices which lie on element edges
                    xe = [NODE(elem.nodes).X]; ye = [NODE(elem.nodes).Y];
                    xVert = [VERT(locVert).X]; yVert = [VERT(locVert).Y];
                    [~,ON] = inpolygon(xVert,yVert,xe,ye);
                    quadVertId = locVert(ON==1);
                    % Coordinates for the triangular sub-element
                    xv = [VERT(locVert).X];
                    yv = [VERT(locVert).Y];
                    %%%% Degenerated coordinates of the sub-triangular element
                    quadNod(1,1).X = xv(1); quadNod(1,1).Y = yv(1);
                    quadNod(2,1).X = xv(1); quadNod(2,1).Y = yv(1);
                    quadNod(3,1).X = xv(2); quadNod(3,1).Y = yv(2);
                    quadNod(4,1).X = xv(3); quadNod(4,1).Y = yv(3);
                    % Coordinates for the quad sub-element
                    [quadCoord] = orderWedgeConn(NODE(physNode),VERT(quadVertId));
                    for iCoord=1:4
                        quadNod(iCoord,2).X = quadCoord(1,iCoord);
                        quadNod(iCoord,2).Y = quadCoord(2,iCoord);
                    end
                    
                end
            elseif(isempty(physNode))
                %%%% Global coordinates of the sub-element
                xv = [VERT(locVert).X];
                yv = [VERT(locVert).Y];
                %%%% Degenerated coordinates of the sub-triangular element
                quadNod(1,1).X = xv(1); quadNod(1,1).Y = yv(1);
                quadNod(2,1).X = xv(1); quadNod(2,1).Y = yv(1);
                quadNod(3,1).X = xv(2); quadNod(3,1).Y = yv(2);
                quadNod(4,1).X = xv(3); quadNod(4,1).Y = yv(3);
            else
                err = MException('tetCutChk:physNodeChk',...
                    'in fourthOrderGaussQuad 3 vertices and %d physical nodes',...
                    length(physNode));
                throw(err);
            end
        elseif(length(locVert)==4)
            if(isempty(physNode))
                % No physical nodes and locally four vertices
                %   --> sub-element is a quad
                %       Coordinates for the quad sub-element
                [quadCoord] = orderWedgeConn(NODE(physNode),VERT(locVert));
                for iCoord=1:4
                    quadNod(iCoord,1).X = quadCoord(1,iCoord);
                    quadNod(iCoord,1).Y = quadCoord(2,iCoord);
                end
            elseif(length(physNode)==1)
                % Locally four vertices and one physical node
                %   --> sub-element is neither a triangle nor a quad
                %   --> Break the sub-element into a triangle and a quad such that
                %      (a) The quad sub-element contains all four vertices and
                %      (b) the tri contains physical node and two vertices which lie
                %          on element edges which meet at the physical node.
                
                % Search for vertices which lie on element edges that meet at
                % the physical node.
                nonPhysNode = elem.nodes(elem.nodes~=physNode); vertId = [];
                for iNode=1:length(nonPhysNode)
                    for iVert=1:length(locVert)
                        xCheck = [NODE(physNode).X NODE(nonPhysNode(iNode)).X VERT(locVert(iVert)).X];
                        yCheck = [NODE(physNode).Y NODE(nonPhysNode(iNode)).Y VERT(locVert(iVert)).Y];
                        AreaCheck = polyarea(xCheck, yCheck);
                        if(AreaCheck<=1e-11)
                            vertId(end+1) = iVert;  %#ok<AGROW> array is small: pre-allocation not essential
                        end
                    end
                end
                % Coordinates for the triangular sub-element
                xp = [NODE(physNode).X];   xv = [VERT(locVert(vertId)).X];
                yp = [NODE(physNode).Y];   yv = [VERT(locVert(vertId)).Y];
                %%%% Degenerated coordinates of the sub-triangular element
                quadNod(1,1).X = xp(1); quadNod(1,1).Y = yp(1);
                quadNod(2,1).X = xp(1); quadNod(2,1).Y = yp(1);
                quadNod(3,1).X = xv(1); quadNod(3,1).Y = yv(1);
                quadNod(4,1).X = xv(2); quadNod(4,1).Y = yv(2);
                % Coordinates for the quad sub-element
                [quadCoord] = orderWedgeConn(NODE(physNode),VERT(locVert));
                for iCoord=1:4
                    quadNod(iCoord,2).X = quadCoord(1,iCoord);
                    quadNod(iCoord,2).Y = quadCoord(2,iCoord);
                end
            else
                err = MException('tetCutChk:physNodeChk',...
                    'in fourthOrderGaussQuad 4 vertices and %d physical nodes',...
                    length(physNode));
                throw(err);
            end
        else
            err = MException('tetCutChk:physNodeChk',...
                'in fourthOrderGaussQuad, wrong number of vertices: %d',...
                length(locVert));
            throw(err);
        end
    end
    
else
    quadNod(4) = struct('X',[],'Y',[]);
    if(isempty(elem.siblings))
        %%%% Global coordinates of the uncut tetrahedral element
        xe = [VERT_ELE(elem.vertices_ele(1:3)).X];
        ye = [VERT_ELE(elem.vertices_ele(1:3)).Y];
        
        %%%% Degenerated coordinates of the uncut tetrahedral element
        quadNod(1,1).X = xe(1); quadNod(1,1).Y = ye(1);
        quadNod(2,1).X = xe(2); quadNod(2,1).Y = ye(2);
        quadNod(3,1).X = xe(3); quadNod(3,1).Y = ye(3);
        quadNod(4,1).X = xe(3); quadNod(4,1).Y = ye(3);
    elseif(length(elem.siblings)==1)
        %%%% For cut elements with just two subdomains
        %%%% Index of the physical node of the cut tetrahedral element
        physVERT = get_physical_vertices(elem,VERT_ELE);
        if(length(elem.vertices)==2)
            if(length(physVERT)==1)
                %%%% Global coordinates of the sub-element
                xp = [VERT_ELE(physVERT).X];   xv = [VERT(elem.vertices).X];
                yp = [VERT_ELE(physVERT).Y];   yv = [VERT(elem.vertices).Y];
                %%%% Degenerated coordinates of the sub-triangular element
                quadNod(1,1).X = xp(1); quadNod(1,1).Y = yp(1);
                quadNod(2,1).X = xp(1); quadNod(2,1).Y = yp(1);
                quadNod(3,1).X = xv(1); quadNod(3,1).Y = yv(1);
                quadNod(4,1).X = xv(2); quadNod(4,1).Y = yv(2);
            elseif(length(physVERT)==2)
                [quadCoord] = orderWedgeConn(VERT_ELE(physVERT),VERT(elem.vertices));
                for iCoord=1:4
                    quadNod(iCoord,1).X = quadCoord(1,iCoord);
                    quadNod(iCoord,1).Y = quadCoord(2,iCoord);
                end
            else
                err = MException('tetCutChk:physNodeChk',...
                    'in fourthOrderGaussQuad 3 vertices and %d physical nodes',...
                    length(physNode));
                throw(err);
            end
        else
            err = MException('tetCutChk:physNodeChk',...
                'in fourthOrderGaussQuad, wrong number of vertices: %d',...
                length(elem.vertices));
            throw(err);
        end
    elseif(length(elem.siblings)>1)
        %%%% For cut elements with more than two subdomains
        error('For cut elements with more than two subdomains: Not implemented yet');
    end
end
end

function[physNode] = get_physical_nodes(elem,NODE)

physNode = zeros(3,1);
pNcount = 0;
for in=1:3
    if(~isempty(NODE(elem.nodes(in)).childrenId))
        pNcount=pNcount+1;
        physNode(pNcount) = elem.nodes(in);
    end
end
physNode(pNcount+1:end) = [];
end

function[physVERT] = get_physical_vertices(elem,VERT_ELE)

physVERT = zeros(3,1);
pNcount = 0;
domain_ele = elem.domain;
for in=1:3
    ls_vert = VERT_ELE(elem.vertices_ele(in)).ls;
    if(ls_vert>0.0)
        domain_vert = 2;
    else
        domain_vert = 1;
    end
    if(domain_ele==domain_vert)
        pNcount=pNcount+1;
        physVERT(pNcount) = elem.vertices_ele(in);
    end
end
physVERT(pNcount+1:end) = [];
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

