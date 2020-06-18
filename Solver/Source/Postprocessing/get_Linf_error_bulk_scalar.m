function[LinfErr] = get_Linf_error_bulk_scalar(NODE,ELEM,PARAMS,u)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
disp = u(1:length(NODE)*PARAMS.ndof); mat_param;

LinfEx = 0.0; LinfErr = 0.0;

for e=1:length(ELEM)
    
    
    
    [ue] = get_nodal_solution(ELEM(e),PARAMS,disp);
    if(isempty(ELEM(e).siblings))
        physNode = ELEM(e).nodes;
    else
        physNode = get_physical_nodes(ELEM(e),NODE);
    end
    
    len = length(physNode);
    
    for i = 1:len
        [uex,~] = exactsolution(NODE(physNode(i)).X,NODE(physNode(i)).Y,ELEM(e).domain,PARAMS);
        uh = ue(find(ELEM(e).nodes==physNode(i)));
        LExElm = abs(uh);
        LErrElm = abs(uex-uh);
        
        if(LinfEx<LExElm)
            LinfEx = LExElm;
        end
        
        if(LinfErr<LErrElm)
            LinfErr = LErrElm;
        end
    end
end
    LinfErr = LinfErr/LinfEx;

   
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





