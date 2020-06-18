function[fBar] = bodyforce(x,y,z,domain,PARAMS)

if(domain == 1)
    fBar = 2*sin(x)*cos(y);
else
    fBar = 0;
end

end