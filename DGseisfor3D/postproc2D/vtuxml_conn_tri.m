function vtuxml_conn_tri(fid,t)

if(size(t,2)==3)
    t=t';
end

Nele=size(t,2);
fwrite(fid,4*Nele*3,'int32');
fwrite(fid,t(:),'int32');

fwrite(fid,Nele*4,'int32');
fwrite(fid,(1:Nele)*3,'int32');
fwrite(fid,Nele,'int32');
fwrite(fid,ones(Nele,1)*5,'int8');

return
end
