function [cond,regions,vars,infile,outfile,tmpdir]=...
    read_visual_para(fname)

    fid=fopen(fname,'r');
    % read number of regions
    tline=read_line(fid);
    a=str2num(tline);
    Nregion=a(1);
    if(Nregion==0)
        regions=[];
    else
        tline=read_line(fid);
        a=str2num(tline);
        regions=a(1:Nregion);
    end
    % read number of Conditions
    tline=read_line(fid);
    a=str2num(tline);
    Ncond=a(1);
    cond=cell(Ncond,1);
    for i=1:Ncond
        % read number of subConditions
        tline=read_line(fid);
        a=str2num(tline);
        Nsubcond=a(1);
        cond{i}=cell(Nsubcond,1);
        for j=1:Nsubcond
            % read type
            tline=read_line(fid);
            tline=strtrim(tline);
            k=find(tline==' ',1,'first');
            a=tline(1:k-1);tline(1:k-1)=' ';
            cond{i}{j}.type=str2double(a);
            % read normal/center
            cond{i}{j}.value=zeros(3,1);
            for m=1:3
                tline=strtrim(tline);
                k=find(tline==' ',1,'first');
                a=tline(1:k-1);tline(1:k-1)=' ';
                cond{i}{j}.value(m)=str2double(a);
            end
            % read offset/radius
            a=str2num(tline);
            cond{i}{j}.offset=a(1);
        end
    end
    % read time invarient scalar
    vars.scalvar=read_varname(fid);
    % read time invarient vector
    vars.vectvar=read_varname(fid);
    % read snap times
    tline=read_line(fid);
    a=str2num(tline);
    vars.snap=(a(1):a(3):a(2));
    % read time varient scalar
    vars.scalvarT=read_varname(fid);
    % read time varient vector
    vars.vectvarT=read_varname(fid);

    tline=read_line(fid);
    infile=strtrim(tline);
    tline=read_line(fid);
    outfile=strtrim(tline);
    tline=read_line(fid);
    tmpdir=strtrim(tline);
    if(tmpdir(end)~='/')
        tmpdir=[tmpdir,'/'];
    end
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

function var=read_varname(fid)
    tline=read_line(fid);
    a=str2num(tline);
    Nvar=a(1);
    var=cell(Nvar,1);
    for i=1:Nvar
        tline=strtrim(read_line(fid));
        k=find(tline==' ',1,'first');
        var{i}.varname=tline(1:k-1);
        tline(1:k-1)=' ';
        a=strtrim(tline);
        var{i}.flag=strcmp(a,'on');
    end
end 
