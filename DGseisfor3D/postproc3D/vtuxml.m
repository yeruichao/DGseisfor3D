function vtuxml(filename,p,t,s,v)
t=t-1;
Nele=size(t,2);
Nnode=size(p,2);
Nscal=length(s)/Nnode;
Nvect=size(v,2)/Nnode;
fid=fopen([filename,'.vtu'],'w');
fprintf(fid,'<?xml version="1.0"?>\n');
fprintf(fid,'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
fprintf(fid,'  <UnstructuredGrid>\n');
fprintf(fid,'    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n',Nnode,Nele);
fprintf(fid,'      <Points>\n');
fprintf(fid,'        <DataArray type="Float32" NumberOfComponents="3" Name="Points" format="appended" offset="0"/>\n');
fprintf(fid,'      </Points>\n');
fprintf(fid,'      <Cells>\n');
offset=3*Nnode*4+4;
fprintf(fid,'        <DataArray type="Int32" Name="connectivity" format="appended" offset="%d"/>\n',offset);
offset=offset+4*Nele*4+4;
fprintf(fid,'        <DataArray type="Int32" Name="offsets" format="appended" offset="%d"/>\n',offset);
offset=offset+4*Nele+4;
fprintf(fid,'        <DataArray type="Int8" Name="types" format="appended" offset="%d"/>\n',offset);
offset=offset+Nele+4;
fprintf(fid,'      </Cells>\n');
fprintf(fid,'      <PointData>\n');
for i=1:Nscal
    fprintf(fid,['        <DataArray type="Float32" Name="scalars',num2str(i),'" NumberOfComponents="1" format="appended" offset="%d"/>\n'],offset);
    offset=offset+Nnode*4+4;
end
for i=1:Nvect
    fprintf(fid,['        <DataArray type="Float32" Name="vector',num2str(i),'" NumberOfComponents="3" format="appended" offset="%d"/>\n'],offset);
    offset=offset+3*Nnode*4+4;
end
fprintf(fid,'      </PointData>\n');
fprintf(fid,'    </Piece>\n');
fprintf(fid,'  </UnstructuredGrid>\n');
fprintf(fid,'  <AppendedData encoding="raw">\n');
fwrite(fid,'_','char*1')
fwrite(fid,3*Nnode*4,'int32');
fwrite(fid,p(:),'float32');
fwrite(fid,4*Nele*4,'int32');
fwrite(fid,t(:),'int32');
fwrite(fid,Nele*4,'int32');
fwrite(fid,(1:Nele)*4,'int32');
fwrite(fid,Nele,'int32');
fwrite(fid,ones(Nele,1)*10,'int8');
for i=1:Nscal
    fwrite(fid,Nnode*4,'int32');
    fwrite(fid,s((i-1)*Nnode+1:i*Nnode),'float32');
end
for i=1:Nvect
    fwrite(fid,3*Nnode*4,'int32');
    vtmp=v(:,(i-1)*Nnode+1:i*Nnode);
    fwrite(fid,vtmp(:),'float32');
end
fwrite(fid,char(10),'char*1');
fprintf(fid,'  </AppendedData>\n');
fprintf(fid,'</VTKFile>');
fclose(fid);

end
