function [ndime,nnode,nelem,nelnd,npres,ntrac,mate,coor,conn,pres,trac] = Read_nl(infile)
cellarray = textscan(infile,'%s');
ndime = str2num(cellarray{1}{3});
mate = zeros(3,1);
mate(1) = str2num(cellarray{1}{6});
mate(2) = str2num(cellarray{1}{8});
mate(3) = str2num(cellarray{1}{10});
mate(4) = str2num(cellarray{1}{12});
mate(5) = str2num(cellarray{1}{14});
mate(6) = str2num(cellarray{1}{16});
mate(7) = str2num(cellarray{1}{18});
mate(8) = str2num(cellarray{1}{20});
mate(9) = str2num(cellarray{1}{22});
mate(10) = str2num(cellarray{1}{24});
mate(11) = str2num(cellarray{1}{26});
mate(12) = str2num(cellarray{1}{28});
mate(13) = str2num(cellarray{1}{30});
mate(14) = str2num(cellarray{1}{32});
mate(15) = str2num(cellarray{1}{34});
mate(16) = str2num(cellarray{1}{36});
nnode = str2num(cellarray{1}{39});
% #node 
ind = 40;
coor = zeros(ndime,nnode);
%全節點位置
for i = 1:nnode
    for j = 1:ndime
        ind = ind+1;
        coor(j,i) = str2num(cellarray{1}{ind});
    end
end
ind = ind+3;
nelem = str2num(cellarray{1}{ind});
%#元素
ind = ind+2;
nelnd = str2num(cellarray{1}{ind});
%一元素的節點數
conn = zeros(nelnd,nelem);
%一元素內所有的節點編號
ind= ind+1;
for i = 1:nelem
    for j = 1:nelnd
        ind = ind+1;
        conn(j,i) = str2num(cellarray{1}{ind});
    end
end
ind = ind+3;
npres = str2num(cellarray{1}{ind});
%位移邊界條件
ind= ind+1;
pres = zeros(3,npres);
for i = 1:npres
    for j = 1:3
        ind = ind+1;
        pres(j,i) = str2num(cellarray{1}{ind});
    end
end
ind = ind+2;
ntrac = str2num(cellarray{1}{ind});
%外力邊界條件
ind= ind+1;
trac = zeros(2+ndime,ntrac);
for i = 1:ntrac
    for j = 1:2+ndime
        ind = ind+1;
        trac(j,i) = str2num(cellarray{1}{ind});
    end
end