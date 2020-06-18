function [ls ] = get_level_set(x,y,PARAMS)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%ls = x.^2 + y.^2 - 1/4; 

ls1 = (x-0.25).^2+ (y-0.25).^2- 0.1^2;
ls2 = (x-0.75).^2+ (y-0.75).^2- 0.1^2;
ls3 = (x-0.25).^2+ (y-0.75).^2- 0.1^2;
ls4 = (x-0.75).^2+ (y-0.25).^2- 0.1^2;

ls5 = (x+0.25).^2+ (y-0.25).^2- 0.1^2;
ls6 = (x+0.75).^2+ (y-0.75).^2- 0.1^2;
ls7 = (x+0.25).^2+ (y-0.75).^2- 0.1^2;
ls8 = (x+0.75).^2+ (y-0.25).^2- 0.1^2;

ls9 = (x-0.25).^2+ (y+0.25).^2- 0.1^2;
ls10 = (x-0.75).^2+ (y+0.75).^2- 0.1^2;
ls11 = (x-0.25).^2+ (y+0.75).^2- 0.1^2;
ls12 = (x-0.75).^2+ (y+0.25).^2- 0.1^2;

ls13 = (x+0.25).^2+ (y+0.25).^2- 0.1^2;
ls14 = (x+0.75).^2+ (y+0.75).^2- 0.1^2;
ls15 = (x+0.25).^2+ (y+0.75).^2- 0.1^2;
ls16 = (x+0.75).^2+ (y+0.25).^2- 0.1^2;

ls = min([ls1 ls2 ls3 ls4 ls5 ls6 ls7 ls8 ls9 ls10 ls11 ls12 ls13 ls14 ls15 ls16],[],2);

end

