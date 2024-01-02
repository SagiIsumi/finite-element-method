function strain = postproc(iel,ndime,nelnd,conn,coor,uglob)
    strain=zeros(3,1);
    
    u=zeros(6,1);
    M = numIntegPt(ndime,nelnd);
    coorie = zeros(ndime,nelnd);
    xii = zeros(ndime,1);
    xi = IntegPt(ndime,nelnd,M);
    for a = 1:nelnd
        for i = 1:ndime
            coorie(i,a) = coor(i,conn(a,iel));
        end
    end
        for i = 1:ndime
            xii(i) = xi(i,1);
        end
        dNdxi = ShpFuncDeri(nelnd,ndime,xii);
        dxdxi = zeros(ndime,ndime);
        for i = 1:ndime
            for j = 1:ndime
                for a = 1:nelnd
                    dxdxi(i,j) = dxdxi(i,j)+coorie(i,a)*dNdxi(a,j);
                end
            end
        end
        dxidx = inv(dxdxi);
        dNdx = zeros(nelnd,ndime);
        for a = 1:nelnd
            for i = 1:ndime
                for j = 1:ndime
                    dNdx(a,i) = dNdx(a,i)+dNdxi(a,j)*dxidx(j,i);
                end
            end
        end
        B = zeros(3,8);
        B(1,1) = dNdx(1,1);
        B(1,3) = dNdx(2,1);
        B(1,5) =  dNdx(3,1);
        B(1,7) =  dNdx(4,1);
        B(2,2) =dNdx(1,2);
        B(2,4) =  dNdx(2,2);
        B(2,6) = dNdx(3,2);
        B(2,8) =  dNdx(4,2);
        B(3,1) =dNdx(1,2);
        B(3,2) =dNdx(1,1);
        B(3,3) =dNdx(2,2);
        B(3,4) = dNdx(2,1);
        B(3,5) = dNdx(3,2);
        B(3,6) = dNdx(3,1);
        B(3,7) = dNdx(4,2);
        B(3,8) = dNdx(4,1);
        u(1,1) = uglob(2*conn(1,iel)-1,1);
        u(2,1) = uglob(2*conn(1,iel),1);
        u(3,1) = uglob(2*conn(2,iel)-1,1);
        u(4,1) = uglob(2*conn(2,iel),1);
        u(5,1) = uglob(2*conn(3,iel)-1,1);
        u(6,1) = uglob(2*conn(3,iel),1);
        u(7,1) = uglob(2*conn(4,iel)-1,1);
        u(8,1) = uglob(2*conn(4,iel),1);
        temp=B*u;
        strain(1,1)=temp(1,1);
        strain(2,1)=temp(2,1);
        strain(3,1)=temp(3,1);

end

