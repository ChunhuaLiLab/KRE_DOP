function [ce1]=pdb_xyz()
filename='4n6h_A.pdb';
fid1=fopen(filename);
MAX0=20000;
resall=[];
posall=[];
ca_b1=[];
resseq=[];
c=[];
c1=[];
ce=[];
for i=1:MAX0
    pdb=fgetl(fid1);
    if pdb(1)=='E'&pdb(2)=='N'&pdb(3)=='D'
        break;
    elseif length(pdb)<3
        continue;
    elseif (pdb(1:4)=='ATOM'& pdb(14)=='C' & pdb(15)=='A')  %% ¶ÁÈ¡AÁ´Ô­×Ó
        position=str2num(pdb(30:54));
        posall=[posall;position];
    end
end
ce1=posall;
fclose('all');
end