function vtuxml_var(fid,var)

%size(var)

if(size(var,2)==3)
    var=var';
end

var=var(:);
%[min(var),max(var)]

fwrite(fid,length(var)*4,'int32');
fwrite(fid,var,'float32');

return
end
