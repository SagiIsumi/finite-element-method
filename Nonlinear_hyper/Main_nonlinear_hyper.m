addpath('D:\isumi\大四\有限元素法\matlab_FEM\Nonlinear_hyper\test_data')

Myfile=fopen("6-3.txt");
[ndime,nnode,nelem,nelnd,npres,ntrac,mate,coor,conn,pres,trac]=Read_nl(Myfile);
rglob = GlobTrac(ndime,nnode,nelem,nelnd,ntrac,coor,conn,trac);
fclose(Myfile);
nstpho = mate(13);
nprtho=mate(14);
tol = 1e-4;%誤差
maxit = 100;%最大迭代數
wglob = zeros(nnode*ndime,1);%假設的u值
filename='D:\isumi\大四\有限元素法\matlab_FEM\Nonlinear\test_data\6-2_output0.vtk';
iprt=1;
count=1;
plotdata=zeros(2,nprtho+1);
%MakeVTK(wglob,ndime,nnode,nelem,nelnd,coor,conn,filename);
for i = 2:nstpho+1 % 根據時間步數變換
    faci = (i-1)/nstpho;
    erri = tol*100;
    niti = 0;%現在迭代次數
    disp(['Step: ' num2str(i) ', factor: ' num2str(faci)]);
    while(erri > tol && niti < maxit) %迭代至小於誤差
        niti = niti+1;
        kglob = GlobStif(ndime,nnode,nelem,nelnd,mate,coor,conn,wglob);
        fglob = GlobResi(ndime,nnode,nelem,nelnd,mate,coor,conn,wglob);
        bglob = faci*rglob-fglob;
        for j = 1:npres
            ir = ndime*(pres(1,j)-1)+pres(2,j);
            for ic = 1:ndime*nnode
                kglob(ir,ic) = 0;
            end
            kglob(ir,ir) = 1;
            bglob(ir) = -wglob(ir) + faci*pres(3,j);
        end
        dwglob = kglob\bglob;
        dwglobsq = dwglob.'*dwglob;
        wglob = wglob + dwglob;
        wglobsq = wglob.'*wglob;
        erri = sqrt(dwglobsq/wglobsq);
        disp(['Iter: ' num2str(niti) ', err: ' num2str(erri)]);
    end
     if (iprt==nprtho)
        iprt=0;
        count=count+1;
        filename=['D:\isumi\大四\有限元素法\matlab_FEM\Nonlinear\test_data\6-2_output',num2str(count),'.vtk'];
        %MakeVTK(wglob,ndime,nnode,nelem,nelnd,coor,conn,filename);
        ee=abs(wglob(84)/10);
        sig=Stress(20,ndime,nelnd,coor,conn,mate,wglob);
        sigsq=sig.'*sig;
        plotdata(1,count)=sqrt(sigsq(2,2));
        plotdata(2,count)=ee;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%繪圖%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iprt=iprt+1;

end
a1=nexttile;
plot(a1,plotdata(2,:),plotdata(1,:),'b-x');
ylabel(a1,'strees');
xlabel(a1,'strain');






