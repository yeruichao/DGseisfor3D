function fid=vtuxml_head_tri(filename,Nnode,Nele,Scalname,Vectname)

Nscal=size(Scalname,1);Nvect=size(Vectname,1);

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
offset=offset+4*Nele*3+4;
fprintf(fid,'        <DataArray type="Int32" Name="offsets" format="appended" offset="%d"/>\n',offset);
offset=offset+4*Nele+4;
fprintf(fid,'        <DataArray type="Int8" Name="types" format="appended" offset="%d"/>\n',offset);
offset=offset+Nele+4;
fprintf(fid,'      </Cells>\n');
fprintf(fid,'      <PointData>\n');
for i=1:Nscal
    fprintf(fid,['        <DataArray type="Float32" Name="',Scalname{i},'" NumberOfComponents="1" format="appended" offset="%d"/>\n'],offset);
    offset=offset+Nnode*4+4;
end
for i=1:Nvect
    fprintf(fid,['        <DataArray type="Float32" Name="',Vectname{i},'" NumberOfComponents="3" format="appended" offset="%d"/>\n'],offset);
    offset=offset+3*Nnode*4+4;
end
fprintf(fid,'      </PointData>\n');
fprintf(fid,'    </Piece>\n');
fprintf(fid,'  </UnstructuredGrid>\n');
fprintf(fid,'  <AppendedData encoding="raw">\n');
fwrite(fid,'_','char*1');

return
end
