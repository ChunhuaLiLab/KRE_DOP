function[]=fcfGNMMD()
[ce1]=pdb_xyz();
posall=ce1;
md_cof=xlsread('md_cof.xlsx');
md_msf=xlsread('md_msf.xlsx');
md_dccm=csvread('md_cij.csv');
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
    end
end
for i=1:n
    netmat(i,i)=-1*sum(netmat(i,:));
end
[V,D]=eig(netmat);
%%the dissipated work
      d=diag(D);
      WW=[];
      W1=0;
     for w=0:40 
          W1=0;
         for p=1:n
             W(p)=0;
           for k=2:n
            J2(k)=w/((d(k)*d(k))+(w*w));
             W(p)=W(p)+V(p,k)*V(p,k)*J2(k);
           end
          W1=W1+W(p);
         end
         WW=[WW W1];
      end
     W=WW';
 
%the synchronous and asynchronous correlations
p=60;%95-35=60
d=diag(D);
 w=d(2);
for k=2:n
    for m=2:n
        Q(k,m)=(d(k)*d(m)+w*w)/(2*(d(k)*d(k)+w*w)*(d(m)*d(m)+w*w));
    end
end
for i=1:n
    for j=1:n
        A1(i,j)=0; 
        for k=2:n
            for m=2:n
                A1(i,j)=A1(i,j)+Q(k,m)*V(p,k)*V(p,m)*V(i,k)*V(j,m);
            end
        end
    end
end

for k=2:n
    for m=2:n
        H(k,m)=((d(k)-d(m))*w)/(2*(d(k)*d(k)+w*w)*(d(m)*d(m)+w*w));
    end
end

for i=1:n
    for j=1:n
        A2(i,j)=0; 
        for k=2:n
            for m=2:n
                A2(i,j)=A2(i,j)+H(k,m)*V(p,k)*V(p,m)*V(i,k)*V(j,m);
            end
        end
    end
end 
a=A1(p,:)';
figure(1);
x=36:338;
plot(x,a);
hold on ;
b=A2(p,:)';
x=36:338;
plot(x,b);
axis([36,338,-0.03,0.62]);
end