function [normal_ls] = normal(coeff)

%%%%%%%%%%% Calculate normal to the interface %%%%%%%%%%%%%%%%%

norm_grad = sqrt (coeff(2)*coeff(2) + coeff(3)*coeff(3));
normal_ls(1) = -coeff(2)/norm_grad;
normal_ls(2) = -coeff(3)/norm_grad;

normal_ls = [normal_ls(1), normal_ls(2)];