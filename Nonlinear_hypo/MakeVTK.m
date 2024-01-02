function MakeVTK(u,ndime,nnode,nelem,nelnd,coor,conn,filename)
    Myfile=fopen(filename,'w');
    fprintf(Myfile,'# vtk DataFile Version 2.0\r\nGenerated Volume Mesh\r\n');
    fprintf(Myfile,'ASCII\r\nDATASET UNSTRUCTURED_GRID\r\n');
    fprintf(Myfile,'POINTS %d float\r\n',nnode);
    uarr=zeros(1,nnode);
    
    if(ndime==2)
        uarr(1,:)=u(1:2:2*nnode-1);
        uarr(2,:)=u(2:2:2*nnode);
        coor(3,:)=0;
        coor(1,:)=coor(1,:)+uarr(1,:);
        coor(2,:)=coor(2,:)+uarr(2,:);
        fprintf(Myfile,'%f %f %f\r\n',coor);
    else
        uarr(1,:)=u(1:3:3*nnode-2);
        uarr(2,:)=u(2:3:3*nnode-1);
        uarr(3,:)=u(3:3:3*nnode);
        coor(1,:)=coor(1,:)+uarr(1,:);
        coor(2,:)=coor(2,:)+uarr(2,:);
        coor(3,:)=coor(3,:)+uarr(3,:);
        fprintf(Myfile,'%f %f %f\r\n',coor);
    end
    a1=nelem*(nelnd+1);
    fprintf(Myfile,'CELLS %d %d\r\n',nelem,a1);
    conn=conn(:,:)-1;
    elem=zeros(1,nelem);
    elem(:)=nelnd;
    for i = 1 : nelnd
        elem(i+1,:)=conn(i,:);
    end
    fprintf(Myfile,'%d %d %d %d %d\r\n',elem);
    fprintf(Myfile,'CELL_TYPES %d\r\n',nelem);
    arr1=zeros(1,nelem);
    arr1(:)=24;%%%%%%% Paraview elements type
    fprintf(Myfile,'%d\r\n',arr1);
    fprintf(Myfile,'POINT_DATA %d\r\n',nnode);
    fprintf(Myfile,'SCALARS var1 float\r\nLOOKUP_TABLE default\r\n');
    arr2=0:0.1:49.8;
    fprintf(Myfile,'%f\r\n',arr2);
end

