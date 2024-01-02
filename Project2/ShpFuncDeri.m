function dNdxi = ShpFuncDeri(nelnd,ndime,xii)
    dNdxi = zeros(nelnd,ndime);
    if (ndime == 1)
        dNdxi(1,1) = 0.5;
        dNdxi(2,1) = -0.5;
    elseif (ndime == 2)
        if (nelnd == 3)
            dNdxi(1,1) = 1;
            dNdxi(2,2) = 1;
            dNdxi(3,1) = -1;
            dNdxi(3,2) = -1;
        end
        if (nelnd == 6)
            dNdxi(1,1) = 4*xii(1)-1;
            dNdxi(2,2) = 4*xii(2)-1;
            dNdxi(3,1) = -(4*(1-xii(1)-xii(2))-1);
            dNdxi(3,2) = -(4*(1-xii(1)-xii(2))-1);
            dNdxi(4,1) = 4*xii(2);
            dNdxi(4,2) = 4*xii(1);
            dNdxi(5,1) = -4*xii(2);
            dNdxi(5,2) = 4*(1-xii(1)-xii(2))-4*xii(2);
            dNdxi(6,1) = 4*(1-xii(1)-xii(2))-4*xii(1);
            dNdxi(6,2) = -4*xii(1);
        end
        if (nelnd == 4)
            dNdxi(1,1) = -0.25*(1-xii(2));
            dNdxi(1,2) = -0.25*(1-xii(1));
            dNdxi(2,1) = 0.25*(1-xii(2));
            dNdxi(2,2) = -0.25*(1+xii(1));
            dNdxi(3,1) = 0.25*(1+xii(2));
            dNdxi(3,2) = 0.25*(1+xii(1));
            dNdxi(4,1) = -0.25*(1+xii(2));
            dNdxi(4,2) = 0.25*(1-xii(1));
        end
        if (nelnd == 8)
            dNdxi(1,1) = 0.25*(1-xii(2))*(2*xii(1)+xii(2));
            dNdxi(1,2) = 0.25*(1-xii(1))*(xii(1)+2*xii(2));
            dNdxi(2,1) = 0.25*(1-xii(2))*(2*xii(1)-xii(2));
            dNdxi(2,2) = 0.25*(1+xii(1))*(2*xii(2)-xii(1));
            dNdxi(3,1) = 0.25*(1+xii(2))*(2*xii(1)+xii(2));
            dNdxi(3,2) = 0.25*(1+xii(1))*(2*xii(2)+xii(1));
            dNdxi(4,1) = 0.25*(1+xii(2))*(2*xii(1)-xii(2));
            dNdxi(4,2) = 0.25*(1-xii(1))*(2*xii(2)-xii(1));
            dNdxi(5,1) = -xii(1)*(1-xii(2));
            dNdxi(5,2) = -0.5*(1-xii(1)*xii(1));
            dNdxi(6,1) = 0.5*(1-xii(2)*xii(2));
            dNdxi(6,2) = -(1+xii(1))*xii(2);
            dNdxi(7,1) = -xii(1)*(1+xii(2));
            dNdxi(7,2) = 0.5*(1-xii(1)*xii(1));
            dNdxi(8,1) = -0.5*(1-xii(2)*xii(2));
            dNdxi(8,2) = -(1-xii(1))*xii(2);
        end
    elseif (ndime == 3)
        if (nelnd == 4)
            dNdxi(1,1) = 1;
            dNdxi(2,2) = 1;
            dNdxi(3,3) = 1;
            dNdxi(4,1) = -1;
            dNdxi(4,2) = -1;
            dNdxi(4,3) = -1;
        end
        if (nelnd == 10)
            xii4=1-xii(1)-xii(2)-xii(3);
            dNdxi(1,1) = 4*xii(1)-1;
            dNdxi(2,2) = 4*xii(2)-1;
            dNdxi(3,3) = 4*xii(3)-1;
            dNdxi(4,1) = -4*xii4+1;
            dNdxi(4,2) = -4*xii4+1;
            dNdxi(4,3) = -4*xii4+1;
            dNdxi(5,1) = 4*xii(2);
            dNdxi(5,2) = 4*xii(1);
            dNdxi(6,2) = 4*xii(3);
            dNdxi(6,3) = 4*xii(2);
            dNdxi(7,1) = 4*xii(3);
            dNdxi(7,3) = 4*xii(1);
            dNdxi(8,1) = 4*(xii4-xii(1));
            dNdxi(8,2) = -4*xii(1);
            dNdxi(8,3) = -4*xii(1);
            dNdxi(9,1) = -4*xii(2);
            dNdxi(9,2) = 4*(xii4-xii(2));
            dNdxi(9,3) = -4*xii(2);
            dNdxi(10,1) = -4*xii(3);
            dNdxi(10,2) = -4*xii(3);
            dNdxi(10,3) = 4*(xii4-xii(3));
           
        end
        if (nelnd == 8)
            dNdxi(1,1) = -1*(1-xii(2))*(1-xii(3))/8;
            dNdxi(1,2) = -1*(1-xii(1))*(1-xii(3))/8;
            dNdxi(1,3) = -1*(1-xii(2))*(1-xii(1))/8;
            dNdxi(2,1) = (1-xii(2))*(1-xii(3))/8;
            dNdxi(2,2) = -1*(1+xii(1))*(1-xii(3))/8;
            dNdxi(2,3) = -1*(1+xii(1))*(1-xii(2))/8;
            dNdxi(3,1) = (1+xii(2))*(1-xii(3))/8;
            dNdxi(3,2) = (1+xii(1))*(1-xii(3))/8;
            dNdxi(3,3) = -1*(1+xii(1))*(1+xii(2))/8;
            dNdxi(4,1) = -1*(1+xii(2))*(1-xii(3))/8;
            dNdxi(4,2) = (1-xii(1))*(1-xii(3))/8;
            dNdxi(4,3) = -1*(1+xii(2))*(1-xii(1))/8;
            dNdxi(5,1) = -1*(1+xii(3))*(1-xii(2))/8;
            dNdxi(5,2) = -1*(1+xii(3))*(1-xii(1))/8;
            dNdxi(5,3) = (1-xii(1))*(1-xii(2))/8;
            dNdxi(6,1) = (1-xii(2))*(1+xii(3))/8;
            dNdxi(6,2) = -1*(1+xii(1))*(1+xii(3))/8;
            dNdxi(6,3) = (1+xii(1))*(1-xii(2))/8;
            dNdxi(7,1) = (1+xii(3))*(1+xii(2))/8;
            dNdxi(7,2) = (1+xii(1))*(1+xii(3))/8;
            dNdxi(7,3) = (1+xii(1))*(1+xii(2))/8;
            dNdxi(8,1) = -1*(1+xii(2))*(1+xii(3))/8;
            dNdxi(8,2) = (1-xii(1))*(1+xii(3))/8;
            dNdxi(8,3) = (1+xii(2))*(1-xii(1))/8;
        end
        if (nelnd == 20)
            a=0;
        end
    end
end

