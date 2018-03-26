function tt=tetsubele(N)

if(N==1)
    tt=1:4;
    return
end

if(N == 2) 
    Nsele = N*(N+1)*(N+2)/6 + (N-1)*N*(N+1)/6*4;
elseif(N >= 3) 
    Nsele = N*(N+1)*(N+2)/6 + (N-1)*N*(N+1)/6*4 + (N-2)*(N-1)*N/6;
end
tt=zeros(Nsele,4);
M=zeros(N+1,N+1,N+1);

sk=0;
for k=1:N+1
    for j=1:N+2-k
        for i=1:N+3-k-j
            sk=sk+1;
            M(i,j,k)=sk;
        end
    end
end
sk=0;
for k=1:N
    for j=1:N+1-k
        for i=1:N+2-k-j
            sk=sk+1;
            tt(sk,1)=M(i,j,k);
            tt(sk,2)=M(i+1,j,k);
            tt(sk,3)=M(i,j+1,k);
            tt(sk,4)=M(i,j,k+1);
        end
    end
end
for k=1:N-1
    for j=1:N-k
        for i=1:N+1-k-j
            t1=M(i+1,j,k);t2=M(i+1,j+1,k);t3=M(i,j+1,k);
            t4=M(i,j,k+1);t5=M(i+1,j,k+1);t6=M(i,j+1,k+1);
            sk=sk+1;
            tt(sk,1)=t4;tt(sk,2)=t2;tt(sk,3)=t5;tt(sk,4)=t1;
            sk=sk+1;
            tt(sk,1)=t4;tt(sk,2)=t2;tt(sk,3)=t1;tt(sk,4)=t3;
            sk=sk+1;
            tt(sk,1)=t4;tt(sk,2)=t2;tt(sk,3)=t3;tt(sk,4)=t6;
            sk=sk+1;
            tt(sk,1)=t4;tt(sk,2)=t2;tt(sk,3)=t6;tt(sk,4)=t5;
        end
    end
end
if(N==2)
    return
end
for k=1:N-2
    for j=1:N-k-1
        for i=1:N-k-j
            sk=sk+1;
            tt(sk,1)=M(i,j+1,k+1);
            tt(sk,2)=M(i+1,j+1,k+1);
            tt(sk,3)=M(i+1,j,k+1);
            tt(sk,4)=M(i+1,j+1,k);
        end
    end
end

end 
