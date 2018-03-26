function vtuxml_conn(fid,t)

if(size(t,2)==4)
    t=t';
end

Nele=size(t,2);
fwrite(fid,4*Nele*4,'int32');
fwrite(fid,t(:),'int32');

fwrite(fid,Nele*4,'int32');
fwrite(fid,(1:Nele)*4,'int32');
fwrite(fid,Nele,'int32');
fwrite(fid,ones(Nele,1)*10,'int8');

return
end
