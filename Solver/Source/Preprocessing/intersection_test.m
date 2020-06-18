function [int_flag] = intersection_test(ls,node,e,n)
%%%%%%% Flags indicating intersection with elements %%%%%
int_flag = 0;    %%%% Flag for a cut element

for i=2:n
    schange=ls(node(i-1,e))*ls(node(i,e));   %%%%% Sign change in the level set
    if (schange<=0)
        int_flag = 1;
    end
end
