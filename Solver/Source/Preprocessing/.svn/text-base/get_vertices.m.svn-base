function[VERT,ELEM,cut_edge_conn,vsize,n_cut_edge] = get_vertices(VERT,ELEM,e,vsize,csize,node,cut_edge_conn,edge_id,nb_int,xint,yint,n_cut_edge)

vertices_e = zeros(1,nb_int);

% Get the vertices
if(n_cut_edge==0)
    for i=1:nb_int
        edge_info_e = get_local_cut_edges(i,e,node,edge_id);
        cut_edge_conn(1:2,n_cut_edge+1) = edge_info_e(1:2);
        VERT(vsize+1).X = xint(i);
        VERT(vsize+1).Y = yint(i);
        vertices_e(i) = vsize+1;
        vsize=vsize+1;
        n_cut_edge=n_cut_edge+1;
    end
else
    for i=1:nb_int
        edge_info_e = get_local_cut_edges(i,e,node,edge_id);
        set=find(cut_edge_conn(1,:) == edge_info_e(1) & cut_edge_conn(2,:) == edge_info_e(2));
        if(isempty(set))
            cut_edge_conn(:,(n_cut_edge+1)) = edge_info_e(1:2);
            VERT(vsize+1).X = xint(i);
            VERT(vsize+1).Y = yint(i);
            vertices_e(i) = vsize+1;
            vsize=vsize+1;
            n_cut_edge=n_cut_edge+1;
        else
            vertices_e(i) = set;
        end
    end
end
vertices_e = sort(vertices_e);
ELEM(csize+1).vertices = vertices_e;
ELEM(csize+2).vertices = vertices_e;

function[edge_info_e] = get_local_cut_edges(i_int,e,node,edge_id)

if(edge_id(i_int)==21)
    edge_info_e(1:2) = sort([node(2,e);node(1,e)]);
elseif(edge_id(i_int)==31)
    edge_info_e(1:2) = sort([node(3,e);node(1,e)]);
elseif(edge_id(i_int)==32)
    edge_info_e(1:2) = sort([node(3,e);node(2,e)]);
end
