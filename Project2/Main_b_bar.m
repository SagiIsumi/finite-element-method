addpath('D:\isumi\大四\有限元素法\matlab_FEM\Project2\mydata')

Myfile=fopen("23mm.txt");
[ndime,nnode,nelem,nelnd,npres,ntrac,mate,coor,conn,pres,trac]=Read_nl(Myfile);
wglob=zeros(nnode*ndime,1);
rglob=zeros(ndime*nnode,1);
fclose(Myfile);
%rglob = GlobTrac(ndime,nnode,nelem,nelnd,ntrac,mate,coor,conn,trac);
fcnode=a;
r_num=length(fcnode);
for i = 1:r_num
    y=coor(2,fcnode(i));
    z=coor(3,fcnode(i));
    idx=ndime*(fcnode(i)-1);
    F=15998.64*0.05*0.05;
    rglob(idx+2)=F*y*sqrt(y*y+z*z);
    rglob(idx+3)=F*z*sqrt(y*y+z*z);
end
rglob=-rglob;
%rglob=zeros(ndime*nnode,1);
nstpho = mate(13);
nprtho=mate(14);
filename='D:\isumi\大四\有限元素法\matlab_FEM\Project2\output\sin0.vtk';
%MakeVTK(wglob,ndime,nnode,nelem,nelnd,coor,conn,filename);
iprt=1;
count=0;
kglob = GlobStif(ndime,nnode,nelem,nelnd,mate,coor,conn,wglob);
for j = 1:npres
    ir = ndime*(pres(1,j)-1)+pres(2,j);
    for ic = 1:ndime*nnode
        kglob(ir,ic) = 0;
    end
    kglob(ir,ir) = 1;
    rglob(ir) =  pres(3,j);
end
check =rank(kglob);
for i = 2:nstpho+1 % 根據時間步數變換
    faci = (i-1)/nstpho;
    rpress=faci*rglob;
    wglob=kglob\rpress;
    if (iprt==nprtho)
        iprt=0;
        count=count+1;
        filename=['D:\isumi\大四\有限元素法\matlab_FEM\Project2\output\sin',num2str(count),'.vtk'];
        %MakeVTK(wglob,ndime,nnode,nelem,nelnd,coor,conn,filename);
    end
    iprt=iprt+1;

end






