function vtk_write(fname,title,xyz,t,p,datatype)

fid=fopen(fname,'w');
Nnode=size(xyz,1);
Nele=size(t,1);
fprintf(fid,'# vtk DataFile Version 2.0  \n');
fprintf(fid,'%s\n',title);
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid,'POINTS %d double\n',Nnode);
for node=1:Nnode
    fprintf(fid,'%f\t%f\t%f\n',xyz(node,:));
end
ele_order=4;
cell_size=Nele*(ele_order+1);
fprintf(fid,'CELLS %d %d\n',Nele,cell_size);
for ele=1:Nele
    fprintf(fid,'\t%d',ele_order);
    for order=1:ele_order
        fprintf(fid,'\t%d',t(ele,order)-1);
    end
    fprintf(fid,'\n');
end
fprintf(fid,'\n');
fprintf(fid,'CELL_TYPES %d\n',Nele);
for ele=1:Nele
    fprintf(fid,'10\n');
end
fprintf(fid,'\n');
if(nargin>=5)
    if(datatype==0)
        fprintf(fid,'POINT_DATA %d\n',Nnode);
        fprintf(fid,'SCALARS value double\n');
        fprintf(fid,'LOOKUP_TABLE default\n');
        for node=1:Nnode
            fprintf(fid,'%f\n',p(node));
        end
    else
        fprintf(fid,'CELL_DATA %d\n',Nele);
        fprintf(fid,'SCALARS value double\n');
        fprintf(fid,'LOOKUP_TABLE default\n');
        for ele=1:Nele
            fprintf(fid,'%f\n',p(ele));
        end
    end
end
fclose(fid);
return;
end