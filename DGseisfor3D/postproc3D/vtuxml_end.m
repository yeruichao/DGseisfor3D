function vtuxml_end(fid)

fwrite(fid,char(10),'char*1');
fprintf(fid,'  </AppendedData>\n');
fprintf(fid,'</VTKFile>');
fclose(fid);

return
end
