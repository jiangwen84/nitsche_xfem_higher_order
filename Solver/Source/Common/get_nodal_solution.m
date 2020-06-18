function[ue] = get_nodal_solution(elem,PARAMS,disp)

ndof = PARAMS.ndof; nlink = PARAMS.nlink;

ue = zeros(nlink*ndof,1);
for jnod=1:nlink
    cg = ndof*(elem.nodes(jnod)-1) + 1;
    ce = ndof*(jnod-1)+1;
    ue(ce:(ce+ndof-1)) = disp(cg:(cg+ndof-1));
end
end
