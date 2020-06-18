function [ls] = calc_level_set(x,y,z,numnod,IntGeom)

switch lower(IntGeom)
    case {'plane'}  %%%%%% Planar level set parallel to xy-plane %%%%%%%%%
        ls = z - 0.4856;
    case {'incplane'} %%%%% Inclined plane %%%%%%%%
        ls = 0.2*x - 0.2*y + z - 0.4856;
    case {'spherori'}
        rad = 0.45;
        Rad = sqrt(x.*x+y.*y+z.*z);
        ls = Rad-rad;
    case {'spherfedkiw'}
        rad = 0.2486;
        Rad = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) + (z-0.5)*(z-0.5));
        ls = Rad-rad;
    case {'popcorn'}
        r0 = 0.25; A = 4; sigma = 0.1;
        for k=1:5
            xk(k) = (r0/sqrt(5))*(2*cos(2*(k-1)*pi/5));
            yk(k) = (r0/sqrt(5))*(2*sin(2*(k-1)*pi/5));
            zk(k) = r0/sqrt(5);
        end
        for k=6:10
            xk(k) = (r0/sqrt(5))*(2*cos((2*(k-6)-1)*pi/5));
            yk(k) = (r0/sqrt(5))*(2*sin((2*(k-6)-1)*pi/5));
            zk(k) = -r0/sqrt(5);
        end
        xk(11) = 0; yk(11) = 0; zk(11) = r0;
        xk(12) = 0; yk(12) = 0; zk(12) = -r0;
        
        for i=1:numnod
            ls(i) = sqrt((x(i)-0.5)^2+(y(i)-0.5)^2+(z(i)-0.5)^2)-r0...
                -A*exp(-((x(i)-0.5-xk(1))^2+(y(i)-0.5-yk(1))^2+(z(i)-0.5-zk(1))^2)/sigma^2)...
                -A*exp(-((x(i)-0.5-xk(2))^2+(y(i)-0.5-yk(2))^2+(z(i)-0.5-zk(2))^2)/sigma^2)...
                -A*exp(-((x(i)-0.5-xk(3))^2+(y(i)-0.5-yk(3))^2+(z(i)-0.5-zk(3))^2)/sigma^2)...
                -A*exp(-((x(i)-0.5-xk(4))^2+(y(i)-0.5-yk(4))^2+(z(i)-0.5-zk(4))^2)/sigma^2)...
                -A*exp(-((x(i)-0.5-xk(5))^2+(y(i)-0.5-yk(5))^2+(z(i)-0.5-zk(5))^2)/sigma^2)...
                -A*exp(-((x(i)-0.5-xk(6))^2+(y(i)-0.5-yk(6))^2+(z(i)-0.5-zk(6))^2)/sigma^2)...
                -A*exp(-((x(i)-0.5-xk(7))^2+(y(i)-0.5-yk(7))^2+(z(i)-0.5-zk(7))^2)/sigma^2)...
                -A*exp(-((x(i)-0.5-xk(8))^2+(y(i)-0.5-yk(8))^2+(z(i)-0.5-zk(8))^2)/sigma^2)...
                -A*exp(-((x(i)-0.5-xk(9))^2+(y(i)-0.5-yk(9))^2+(z(i)-0.5-zk(9))^2)/sigma^2)...
                -A*exp(-((x(i)-0.5-xk(10))^2+(y(i)-0.5-yk(10))^2+(z(i)-0.5-zk(10))^2)/sigma^2)...
                -A*exp(-((x(i)-0.5-xk(11))^2+(y(i)-0.5-yk(11))^2+(z(i)-0.5-zk(11))^2)/sigma^2)...
                -A*exp(-((x(i)-0.5-xk(12))^2+(y(i)-0.5-yk(12))^2+(z(i)-0.5-zk(12))^2)/sigma^2);
        end
end