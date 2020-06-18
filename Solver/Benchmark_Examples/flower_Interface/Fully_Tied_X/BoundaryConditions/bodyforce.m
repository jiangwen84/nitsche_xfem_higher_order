function[fBar] = bodyforce(x,y,domain,PARAMS)

if(domain == 1)
    fBar = -4;
else
    fBar = -16*(x^2 + y^2);
end

end