cd ~/IHEPBox/Projects/GitHub/21cmFAST/src/py21cmfast/_data/HaloProfile/
clear

reload=0
FileName = 'ProfileTable.h'
lm1=0
lm2=16.69
z1=0
z2=35

clc
if reload
    run Data.m
    save tmp Tab
else
    load tmp
end
n=size(Tab);
nz=n(1);
nm=n(2);

delete(FileName)
FileID=fopen(FileName,'w');

fprintf(FileID,'#define Z_axis_Size ');
fprintf(FileID,'%i\n',nz);
fprintf(FileID,'#define M_axis_Size ');
fprintf(FileID,'%i\n',nm);
fprintf(FileID,'#define LgM_axis_min ');
fprintf(FileID,'%f\n',lm1);
fprintf(FileID,'#define LgM_axis_max ');
fprintf(FileID,'%f\n',lm2);
fprintf(FileID,'#define Z_axis_min ');
fprintf(FileID,'%f\n',z1);
fprintf(FileID,'#define Z_axis_max ');
fprintf(FileID,'%f\n',z2);

fprintf(FileID,'double ProfileTable[');
fprintf(FileID,'%i',n(1));
fprintf(FileID,']');
fprintf(FileID,'[');
fprintf(FileID,'%i',n(2));
fprintf(FileID,']={\n');

for zid=1:nz
    for mid=1:nm-1
        fprintf(FileID,'%E, ',Tab(zid,mid));
    end
    if zid==nz
        fprintf(FileID,'%E};\n',Tab(zid,nm));
    else
        fprintf(FileID,'%E,\n',Tab(zid,nm));
    end
    
end

size(Tab)
x=Tab(:,1);

clf
semilogy(x)
max(x)

