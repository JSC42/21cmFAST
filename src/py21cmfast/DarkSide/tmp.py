from PyLab import np, plt
from main import Load_EFF_Data, PrintData_Kernel

n, f = Load_EFF_Data()
file = '/Users/cangtao/Desktop/tmp.c'
PrintData_Kernel(n,f,file)
