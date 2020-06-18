function [uLoc] = get_local_solution(elem,PARAMS,u)

for in=1:length(elem.nodes)
    uLoc(PARAMS.ndof*(in-1)+1:PARAMS.ndof*(in-1)+PARAMS.ndof,1)=...
        u(PARAMS.ndof*(elem.nodes(in)-1)+1:PARAMS.ndof*(elem.nodes(in)-1)+PARAMS.ndof);
end

end