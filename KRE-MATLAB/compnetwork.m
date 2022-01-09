function[]=compnetwork()
 [ce1]=pdb_xyz()
 n1=303;
 posall=ce1;
 netmat=zeros(n1);
 for i=1:n1
     for j=i:n1
         if i==j
             continue;
         else 
		     x1=posall(i,:);
			 x2 = posall(j,:);
		     dis=sqrt((x1(1)-x2(1))^2+(x1(2)-x2(2))^2+(x1(3)-x2(3))^2);
			 cutoff =13;
             if dis<cutoff
                 netmat1(i,j)=1;
                 netmat1(j,i)=1;
             end
         end
     end
end
for i=1:n1
         netmat1(i,i)=0;
end
[NUM]=load('fcfgnmmd_dccm.dat');
for i=1:n1
    for j=1:n1
       w(i,j)=-log(abs(NUM(i,j))+eps);
    end
end
for i=1:n1
    for j=1:n1
        netmat(i,j)=w(i,j).*netmat1(i,j);
    end
end 
A=netmat;

N=size(A,2)-1;  
aver_D=[];   
 for m=1:n1
     D=A;
     D(m,:)=[];
     D(:,m)=[];
	 D(find(D==0))=inf;    
	 for i=1:N           
		 D(i,i)=0;       
	 end   
     for k=1:N
         for i=1:N
             for j=1:N
                 if D(i,j)>D(i,k)+D(k,j)
                     D(i,j)=D(i,k)+D(k,j);
                 end
             end
         end
     end
	 
 aver_D1=sum(sum(D))/(N*(N-1)); 
 aver_D=[aver_D aver_D1];
 end

N0=size(netmat,2);     
D0=netmat;
D0(find(D0==0))=inf;    
for i=1:N0           
 D0(i,i)=0;       
end   
for k=1:N0          
 for i=1:N0
	 for j=1:N0
		 if D0(i,j)>D0(i,k)+D0(k,j)
			D0(i,j)=D0(i,k)+D0(k,j);
		 end
	 end
 end
aver_D0=sum(sum(D0))/(N0*(N0-1)); 

CPL1=aver_D;
CPL0=aver_D0;
N1 = size(CPL1,1);
dif_CPL= [];
for i=1:N1
    dif_CPL(i)= CPL1(i)-CPL0(1);
end
%%%Z-score%%%
aver_dif = mean(dif_CPL);
std_dif = std(dif_CPL);
for i = 1:N1
    Z_score(i) = (dif_CPL(i)-aver_dif)/std_dif;
end
%%%degree%%%
for i=1:n1
    k(i)=0;
    for j=1:n1
        k(i)=k(i)+A(i,j);
    end
end
 
end