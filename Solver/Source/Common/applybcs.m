function [bigk,fext] = applybcs(ifixu,u,bigk,fext)

dof = length(ifixu);

for n=1:dof
    if (ifixu(n) == 1)
        for b=1:dof
            fext(b) = fext(b) - bigk(b,n)*u(n);
        end
        bigk(n,1:dof) = zeros(1,dof);
        bigk(1:dof,n) = zeros(1,dof);
        bigk(n,n) = 1.0;
        fext(n) = u(n);
    end
end
