function tt=tet_subelement(N)

nodal=load('porder2nodal.mat','porder');
fmask3=find( abs(1+nodal.porder{N}.r ...
                  +nodal.porder{N}.s ...
                  +nodal.porder{N}.t) < 1e-7)';

nsele=N*(N+1)/2;
if(N >= 2);nsele=nsele+5*N*(N-1)/2;end
if(N >= 3);nsele=nsele+5*(N*(N+1)*(N-4)/6+N);end

tt=zeros(nsele,4);
sk=0;skm=0;
for i=N:-1:1
    tt(sk+1:sk+i,2)=fmask3(skm+1:skm+i);
    tt(sk+1:sk+i,1)=tt(sk+1:sk+i,2)-1;
    tt(sk+1:sk+i,3)=fmask3(skm+2:skm+i+1);
    tt(sk+1:sk+i,4)=fmask3(skm+i+2:skm+2*i+1);
    sk=sk+i;
    skm=skm+i+1;
end

z0=zeros(nsele,1);

if(N >= 2)
    v1=z0;v2=z0;v3=z0;v4=z0;v5=z0;v6=z0;v7=z0;
    skm=0;skv=0;
    for i=N-1:-1:1
       v1(skv+1:skv+i)=tt(skm+1:skm+i,1)-1;
       v2(skv+1:skv+i)=tt(skm+1:skm+i,1);
       v3(skv+1:skv+i)=tt(skm+1:skm+i,3);
       v4(skv+1:skv+i)=tt(skm+2:skm+i+1,1);
       v5(skv+1:skv+i)=tt(skm+1:skm+i,4)-1;
       v6(skv+1:skv+i)=tt(skm+1:skm+i,4);
       v7(skv+1:skv+i)=tt(skm+2:skm+i+1,4);
       skv=skv+i;skm=skm+i+1;
    end
    tt(sk+1:sk+skv,1)=v1(1:skv);
    tt(sk+1:sk+skv,2)=v2(1:skv);
    tt(sk+1:sk+skv,3)=v4(1:skv);
    tt(sk+1:sk+skv,4)=v5(1:skv);
    sk=sk+skv;
    tt(sk+1:sk+skv,1)=v5(1:skv);
    tt(sk+1:sk+skv,2)=v2(1:skv);
    tt(sk+1:sk+skv,3)=v4(1:skv);
    tt(sk+1:sk+skv,4)=v6(1:skv);
    sk=sk+skv;
    tt(sk+1:sk+skv,1)=v4(1:skv);
    tt(sk+1:sk+skv,2)=v2(1:skv);
    tt(sk+1:sk+skv,3)=v3(1:skv);
    tt(sk+1:sk+skv,4)=v6(1:skv);
    sk=sk+skv;
    tt(sk+1:sk+skv,1)=v4(1:skv);
    tt(sk+1:sk+skv,2)=v6(1:skv);
    tt(sk+1:sk+skv,3)=v3(1:skv);
    tt(sk+1:sk+skv,4)=v7(1:skv);
    sk=sk+skv;
    tt(sk+1:sk+skv,1)=v5(1:skv);
    tt(sk+1:sk+skv,2)=v6(1:skv);
    tt(sk+1:sk+skv,3)=v4(1:skv);
    tt(sk+1:sk+skv,4)=v7(1:skv);
    sk=sk+skv;
end

if(N >= 3)
    v1=z0;v2=z0;v3=z0;v4=z0;v5=z0;v6=z0;v7=z0;v8=z0;
    sp=0;skv=0;
    ei=N+1;
    for i=1:ei
        ej=ei-i+1;
        for j=1:ej
            ek=ej-j+1;
            for k=1:ek
                sp=sp+1;
                if((i <= ei-3) && (j <= ej-3) && (k <= ek-3))
                    skv=skv+1;
                    v1(skv)=sp;
                    v2(skv)=v1(skv)+1;
                    v3(skv)=v2(skv)+ek;
                    v4(skv)=v1(skv)+ek;
                    v5(skv)=v1(skv)+ej*(ej+1)/2+1-j;
                    v6(skv)=v5(skv)+1;
                    v7(skv)=v6(skv)+ek-1;
                    v8(skv)=v5(skv)+ek-1;
                end
            end
        end
    end
    tt(sk+1:sk+skv,1)=v1(1:skv);
    tt(sk+1:sk+skv,2)=v2(1:skv);
    tt(sk+1:sk+skv,3)=v4(1:skv);
    tt(sk+1:sk+skv,4)=v5(1:skv);
    sk=sk+skv;
    tt(sk+1:sk+skv,1)=v3(1:skv);
    tt(sk+1:sk+skv,2)=v4(1:skv);
    tt(sk+1:sk+skv,3)=v2(1:skv);
    tt(sk+1:sk+skv,4)=v7(1:skv);
    sk=sk+skv;
    tt(sk+1:sk+skv,1)=v5(1:skv);
    tt(sk+1:sk+skv,2)=v2(1:skv);
    tt(sk+1:sk+skv,3)=v4(1:skv);
    tt(sk+1:sk+skv,4)=v7(1:skv);
    sk=sk+skv;
    tt(sk+1:sk+skv,1)=v6(1:skv);
    tt(sk+1:sk+skv,2)=v7(1:skv);
    tt(sk+1:sk+skv,3)=v2(1:skv);
    tt(sk+1:sk+skv,4)=v5(1:skv);
    sk=sk+skv;
    tt(sk+1:sk+skv,1)=v8(1:skv);
    tt(sk+1:sk+skv,2)=v5(1:skv);
    tt(sk+1:sk+skv,3)=v4(1:skv);
    tt(sk+1:sk+skv,4)=v7(1:skv);
    sk=sk+skv;
end

end
