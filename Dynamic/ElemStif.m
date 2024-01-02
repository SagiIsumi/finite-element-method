function kel = ElemStif(iel,ndime,nelnd,coor,conn,mate)%Kaibk勁度矩陣
    coorie = zeros(ndime,nelnd);%一個暫存的元素
    xii = zeros(ndime,1);
    dxdxi = zeros(ndime,ndime);
    dNdx = zeros(nelnd,ndime);
    M = numIntegPt(ndime,nelnd);
    xi = IntegPt(ndime,nelnd,M);
    w = integWt(ndime,nelnd,M);
    kuu = zeros(ndime*nelnd,ndime*nelnd);
    kau = zeros(ndime*ndime,ndime*nelnd);
    kua = zeros(ndime*nelnd,ndime*ndime);
    kaa = zeros(ndime*ndime,ndime*ndime);
    for a = 1:nelnd
        for i = 1:ndime
            coorie(i,a) = coor(i,conn(a,iel));%conn(a,iel)a是元素中的第a個節點 iel是第幾號元素 coor依節點編號取x,y值
            %建立該元素節點的所有座標位置
        end
    end
    dNdxi = ShpFuncDeri(nelnd,ndime,xii);
    for i = 1:ndime
        for j = 1:ndime
            for a = 1:nelnd
                dxdxi(i,j) = dxdxi(i,j)+coorie(i,a)*dNdxi(a,j);
            end
        end
    end
    jcb0 = det(dxdxi);
    for im = 1:M % M個積分點,依序對所有積分點運算計算對整體元素影響
        for i = 1:ndime
            xii(i) = xi(i,im);% epsilon值;空間維度決定數量
        end
        dNdxi = ShpFuncDeri(nelnd,ndime,xii);
        dxdxi(:)=0.;
        for i = 1:ndime%i表示對x所有的偏微分，此迴圈表示所有節點對該元素在不同維度的積分
            for j = 1:ndime
                for a = 1:nelnd
                    dxdxi(i,j) = dxdxi(i,j)+coorie(i,a)*dNdxi(a,j);%j表示維度,a表示所有積分點的作用
                end
            end
        end
        dxidx = inv(dxdxi);
        jcb = det(dxdxi);
        dNdx(:)= 0.;
        for a = 1:nelnd
            for i = 1:ndime
                for j = 1:ndime
                    dNdx(a,i) = dNdx(a,i)+dNdxi(a,j)*dxidx(j,i);
                end
            end
        end
        cmat = MatStif(ndime,mate);
        for a = 1:nelnd %在第幾號節點施加作用力
            for i = 1:ndime %決定作用力的面
                for b = 1:nelnd%在第幾號節點產生的形變
                    for k = 1:ndime  %決定形變的面
                        ir = ndime*(a-1)+i;
                        ic = ndime*(b-1)+k;
                        for j = 1:ndime%和i共同決定作用力方向
                            for l = 1:ndime%和k共同決定形變方向
                                kuu(ir,ic)=kuu(ir,ic)+cmat(i,j,k,l)*dNdx(a,j)*dNdx(b,l)*w(im)*jcb;
                            end
                        end
                    end
                 end
            end
        end
        for a = 1:nelnd
            for i = 1:ndime
                for m = 1:ndime
                    for k = 1:ndime
                        ic = ndime*(a-1)+i;
                        ir = ndime*(m-1)+k;
                        for j = 1:ndime
                            for l = 1:ndime
                                tmp = cmat(i,j,k,l)*dNdx(a,j)*(xii(m)*dxidx(m,l)*jcb0/jcb)*w(im)*jcb;
                                kau(ir,ic) = kau(ir,ic)+tmp;
                                kua(ic,ir) = kua(ic,ir)+tmp;
                            end
                        end
                    end
                end
            end
        end
        for m = 1:ndime
            for i = 1:ndime
                for n = 1:ndime
                    for k = 1:ndime
                        ir = ndime*(m-1)+i;
                        ic = ndime*(n-1)+k;
                        for j = 1:ndime
                            for l = 1:ndime
                                temp=cmat(i,j,k,l)*xii(m)*dxidx(m,j)*xii(n)*dxidx(n,l)*jcb0/jcb*jcb0/jcb*w(im)*jcb;
                                kaa(ir,ic) = kaa(ir,ic)+temp;
                            end
                        end
                    end
                end
            end
        end
    end

    kel=kuu-kua*(kaa\kau);
end

