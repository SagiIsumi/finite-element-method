addpath('D:\isumi\大四\有限元素法\matlab_FEM\Dynamic\5-3')

Myfile=fopen("5-3.txt");
[ndime,nnode,nelem,nelnd,npres,ntrac,mate,coor,conn,pres,trac]=Read_dyn(Myfile);
mglob = GlobMass(ndime,nnode,nelem,nelnd,mate,coor,conn);
kglob = GlobStif(ndime,nnode,nelem,nelnd,mate,coor,conn);
rglob = GlobTrac(ndime,nnode,nelem,nelnd,ntrac,mate,coor,conn,trac);
mpres = mglob;
kpres = kglob;
rpres = rglob;
for i = 1:npres
    ir = ndime*(pres(1,i)-1)+pres(2,i);    
    for ic = 1:ndime*nnode
        mpres(ir,ic) = 0;
        mpres(ic,ir) = 0;
        kpres(ir,ic) = 0;
        kpres(ic,ir) = 0;
    end
    mpres(ir,ir) = 1.;
    rpres(ir) = pres(3,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dynamic motion calculate
beta1 = mate(5);
beta2 = mate(6);
dt = mate(7);
nstp = mate(8);
nprt = mate(9);
vc = zeros(nnode*ndime,1);
uc = zeros(nnode*ndime,1);
mkpres = mpres+0.5*beta2*dt*dt*kpres;
ac = mkpres\(-kpres*uc+rpres);
stpnum=nstp/nprt;
stpnum=floor(stpnum);
u_output=zeros(nnode*ndime,1);
count=1;
coorc=coor;
filename='D:\isumi\大四\有限元素法\matlab_FEM\Dynamic\5-3\5-3_output0.vtk';
iprt=1;
MakeVTK(uc,ndime,nnode,nelem,nelnd,coor,conn,filename);
for i = 1:nstp         % output 每nstp的動態變化
    rpres=0;
    an = mkpres\(rpres-kpres*(uc+dt*vc+0.5*(1.-beta2)*dt*dt*ac));
    vn = (vc+dt*(1.-beta1)*ac+dt*beta1*an);
    un = (uc+dt*vc+(1.-beta2)*0.5*dt*dt*ac+0.5*beta2*dt*dt*an);
    ac = an;
    vc = vn;
    uc = un;
    if(iprt == nprt)
        iprt = 0;
        u_output(:,count)=uc;
        filename=['D:\isumi\大四\有限元素法\matlab_FEM\Dynamic\5-3\5-3_output',num2str(count),'.vtk'];
        MakeVTK(uc,ndime,nnode,nelem,nelnd,coor,conn,filename);
        count=count+1;
    end
    iprt = iprt+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Natural Frequency & Mode shapes
mpresroot = sqrtm(mpres);
mpresrootinv = inv(mpresroot);
hpres = mpresrootinv*(kpres*mpresrootinv);
[q,lambda,~] = svd(hpres);
lambdasort = zeros(ndime*nnode,1);
for i = 1:ndime*nnode
    lambdasort(i) = lambda(i,i);
end
lambdasort = sort(lambdasort);
check1=q*lambda;
check2=hpres/q;
check3=q*q.';
nmod = mate(10);
if(ndime == 2)
    nmodrbm = 3;
else
    nmodrbm = 6;
end
u = zeros(nnode*ndime,nmod);
for k = 1:nmod
    for i = 1:nnode*ndime
        if((lambda(i,i)-lambdasort(nmodrbm+k))^2 < 1e-8)
            ipick = i;
            break;
        end
    end
    u(:,k) = mpresrootinv*q(:,ipick);
end





