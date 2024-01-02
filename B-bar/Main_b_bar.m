addpath('D:\isumi\大四\有限元素法\matlab_FEM\B-bar\test_data')

Myfile=fopen("7-1.txt");
[ndime,nnode,nelem,nelnd,npres,ntrac,mate,coor,conn,pres,trac]=Read_nl(Myfile);
wglob=zeros(nnode*ndime,1);
rglob = GlobTrac(ndime,nnode,nelem,nelnd,ntrac,mate,coor,conn,trac);

%rglob=zeros(ndime*nnode,1);
nstpho = mate(13);
nprtho=mate(14);
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
wglob=kglob\rglob;
pos=coor;
pos_xin=coor(1,1:100)+wglob(1:2:199)';
pos_yin=coor(2,1:100)+wglob(2:2:200)';
pos_xout=coor(1,101:200)+wglob(201:2:399)';
pos_yout=coor(2,101:200)+wglob(202:2:400)';
pos_xin(1,101)=pos_xin(1,1);
pos_yin(1,101)=pos_yin(1,1);
pos_xout(1,101)=pos_xout(1,1);
pos_yout(1,101)=pos_yout(1,1);
cc=nexttile;
plot(cc,pos_xin,pos_yin,'b-x',pos_xout,pos_yout,'b-x');
strain=postproc(1,ndime,nelnd,conn,coor,wglob);





