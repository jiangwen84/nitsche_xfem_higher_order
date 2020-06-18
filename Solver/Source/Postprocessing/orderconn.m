function [reordered_vert_conn] = orderconn(vert,vert_conn,normal_ls)

flag = 0;
count = 0;
xint = vert(vert_conn,1); yint = vert(vert_conn,2); zint = vert(vert_conn,3);
reordered_vert_conn = vert_conn;
while(flag==0)
    AB = [xint(2)-xint(1),yint(2)-yint(1),zint(2)-zint(1)];
    BC = [xint(3)-xint(2),yint(3)-yint(2),zint(3)-zint(2)];
    crossp1 = cross(AB,BC);

    CD = [xint(4)-xint(3),yint(4)-yint(3),zint(4)-zint(3)];
    DA = [xint(1)-xint(4),yint(1)-yint(4),zint(1)-zint(4)];
    crossp2 = cross(CD,DA);

    dotnormals = dot(crossp1,crossp2);

    if(dotnormals<0)
        count = count + 1;
        if(count==1)
            vert_tmp = reordered_vert_conn(3);
            reordered_vert_conn(3) = reordered_vert_conn(4);
            reordered_vert_conn(4) = vert_tmp;
            xint = vert(reordered_vert_conn,1);
            yint = vert(reordered_vert_conn,2);
            zint = vert(reordered_vert_conn,3);
%             xtemp = xint(4); ytemp = yint(4); ztemp = zint(4);
%             xint(4) = xint(3); yint(4) = yint(3); zint(4) = zint(3);
%             xint(3) = xtemp; yint(3) = ytemp; zint(3) = ztemp;
        elseif(count==2)
            vert_tmp = reordered_vert_conn(2);
            reordered_vert_conn(2) = reordered_vert_conn(3);
            reordered_vert_conn(3) = vert_tmp;
            xint = vert(reordered_vert_conn,1);
            yint = vert(reordered_vert_conn,2);
            zint = vert(reordered_vert_conn,3);
%             xtemp = xint(3); ytemp = yint(3); ztemp = zint(3);
%             xint(3) = xint(2); yint(3) = yint(2); zint(3) = zint(2);
%             xint(2) = xtemp; yint(2) = ytemp; zint(2) = ztemp;
        end

    else
        flag = 1;
        checknormal = dot(crossp1,normal_ls);
        if (checknormal>0)
            vert_tmp = reordered_vert_conn(2);
            reordered_vert_conn(2) = reordered_vert_conn(4);
            reordered_vert_conn(4) = vert_tmp;
            xint = vert(reordered_vert_conn,1);
            yint = vert(reordered_vert_conn,2);
            zint = vert(reordered_vert_conn,3);
%             xtemp = xint(2); ytemp = yint(2); ztemp = zint(2);
%             xint(2) = xint(4); yint(2) = yint(4); zint(2) = zint(4);
%             xint(4) = xtemp; yint(4) = ytemp; zint(4) = ztemp;
        end
    end
end
