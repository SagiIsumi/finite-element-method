function rglob = GlobTrac(ndime,nnode,nelem,nelnd,ntrac,coor,conn,trac)
    rglob = zeros(ndime*nnode,1);
    nfcnd = numFaceNd(ndime,nelnd);
    coorif = zeros(ndime,nfcnd);
    for j = 1:ntrac
        iel = trac(1,j);
        ifc = trac(2,j);
        fcnd = FaceNd(ndime,nelnd,ifc);
        for a = 1:nfcnd
            for i = 1:ndime
                coorif(i,a) = coor(i,conn(fcnd(a),iel));
            end
        end
        rel = ElemLoad(ndime,nelnd,nfcnd,coorif,trac(:,j));
        for a = 1:nfcnd
            for i = 1:ndime
                ir = (conn(fcnd(a),iel)-1)*ndime+i;
                rglob(ir) = rglob(ir)+rel((a-1)*ndime+i);
            end
        end
    end
end