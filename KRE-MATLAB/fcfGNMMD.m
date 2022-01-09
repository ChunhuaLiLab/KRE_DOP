
function[]=fcfGNMMD()
[ce1]=pdb_xyz();
posall=ce1;
md_cof=xlsread('md_cof.xlsx');
md_msf=xlsread('md_msf.xlsx');
md_dccm=csvread('md_cij.csv');
[U,M]=eig(md_cof);
U1=U';
n=303;
rr=[];
Coll=[];
for cutoff=13 
    for s=5
        for p=10^-7*(10^s) 
            for i=1:n
                mk(i,i)=2/(M(i,i)+sqrt(M(i,i)^2+8*p));
            end
            K=U*mk*U1; 
            k=-K; 
            netmat=zeros(n);
            for i=1:n
                for j=i:n
                    if k(i,j)<0
                        k(i,j)=0;
                        k(j,i)=0;
                    end
                end
            end
            for i=1:n
                for j=i:n
                    if i==j
                        continue;
                    else
                        dis=sqrt(sum((posall(i,:)-posall(j,:)).^2));
                        if dis<=cutoff
                            netmat(i,j)=-k(i,j);
                            netmat(j,i)=-k(i,j);
                        end
                    end
                end
            end
        end
    

for i=1:n
    netmat(i,i)=-1*sum(netmat(i,:));
end
[V,D]=eig(netmat);
d=diag(D);
for i=1:n
    fluu(i)=0;
    for k=2:n
        fluu(i)=fluu(i)+(V(i,k)*V(i,k)/D(k,k));
    end
end
fcfgnmmd_msf=fluu';
%Pearson Correlation Coefficient
[S1]=corrcoef(md_msf,fcfgnmmd_msf);
PCC_msf=S1(1,2);

for i=1:n
    for j=1:n
        flu(i,j)=0;
        for k=2:n
            flu(i,j)=flu(i,j)+(V(i,k)*V(j,k)/D(k,k));
    end
end
fcfgnmmd_cof=flu;
[S2]=corrcoef(md_cof,fcfgnmmd_cof);
PCC_cof=S3(1,2);

fcfgnmmd_dccm=zeros(n);
for i=1:n
    for j=1:n
        fcfgnmmd_dccm(i,j)=flu(i,j)/sqrt(flu(i,i)*flu(j,j));
    end
end
[S3]=corrcoef(fcfgnmmd_dccm,md_dccm);
PCC_dccm=S4(1,2);

CC_msf_cof=(PCC_msf+PCC_cof)/2;
result=[PCC_msf PCC_cof CC_msf_cof PCC_dccm];
rr=[rr;result];
arr=roundn(rr,-4);
end
end
end

  




