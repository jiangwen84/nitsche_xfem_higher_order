function [x,y,node] = readmesh(mesh_file)

fid = fopen(mesh_file, 'r');
[~] = fscanf(fid,'%s',1);
[~] = fscanf(fid,'%s',3);
[~] = fscanf(fid,'%s',1);
[~] = fscanf(fid,'%s',1);
numnod = fscanf(fid,'%d',1);
x = zeros(numnod,1);
y = zeros(numnod,1);

for n=1:numnod
    nodeno = fscanf(fid,'%d',1);
    coord = fscanf(fid, '%lf',3);
    x(nodeno) = coord(1);
    y(nodeno) = coord(2);
end

[~] = fscanf(fid,'%s',1);
[~] = fscanf(fid,'%s',1);

numsideele = fscanf(fid,'%d',1);
init_alloc = 100;
block_size = floor(init_alloc + 0.1*init_alloc);
alloc_size = block_size;
numele = 0;
node=zeros(3,alloc_size);
for e=1:(numsideele+1)
    c1 = str2num(fgetl(fid));
    if(length(c1)>2)
        if (c1(2) == 2)
            numele = numele +1;
            node(1:3,numele) = c1(6:8);
        end
    end
    if(numele + 0.1*alloc_size>alloc_size)
        alloc_size = alloc_size + block_size;
        node = [node zeros(3,block_size)];
    end
end
node(:,numele+1:end)=[];

fclose(fid);