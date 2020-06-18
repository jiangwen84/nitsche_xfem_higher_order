function [ls ] = get_level_set(x,y,PARAMS)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ls = x - (PARAMS.xint*(PARAMS.Xmax-PARAMS.Xmin)+PARAMS.Xmin);


end