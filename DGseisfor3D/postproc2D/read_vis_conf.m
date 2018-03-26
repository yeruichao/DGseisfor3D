function [var,infile,outfile,Np]=read_vis_conf(fname)
    fid=fopen(fname,'r');
    tline=read_line(fid);
    Np=str2num(tline);
    tline=read_line(fid);
    a=str2num(tline);
    Nvar=a(1);
    var=cell(Nvar,1);
    for i=1:Nvar
        tline=strtrim(read_line(fid));
        k=find(tline==' ',1,'first');
        var{i}.varname=tline(1:k-1);
        tline(1:k-1)=' ';
        a=str2num(tline);
        var{i}.start=a(1);
        var{i}.final=a(2);
	if(length(a)>=3)
            var{i}.intv=a(3);
        else
            var{i}.intv=1;
        end
    end
    tline=read_line(fid);
    infile=strtrim(tline);
    tline=read_line(fid);
    outfile=strtrim(tline);
end

function tline=read_line(fid)
while(1)
    tline=fgetl(fid);
    if(isempty(tline) || all(tline==' '))
        continue
    end
    if(tline(1) ~= '#')
        k=find(tline=='=',1,'last');
        if(not(isempty(tline)))
            tline=strtrim(tline(k+1:end));
        end
        return
    end
end
end
