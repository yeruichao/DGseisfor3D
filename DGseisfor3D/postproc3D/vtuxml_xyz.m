function vtuxml_xyz(fid,p)

if(size(p,2)==3)
    p=p';
end

Nnode=size(p,2);
fwrite(fid,3*Nnode*4,'int32');
fwrite(fid,p(:),'float32');

return
end
