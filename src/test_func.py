from utils import unit_symbol, plot_picture,plot_picture_trans
from calculation_module_0129 import ac_analysis,Dc_analysis
from parse import parse
import numpy as np
import sympy as sp
from dc_bias import Dc_bias
from scipy import linalg
#from mpmath import mp

import matplotlib.pyplot as plt


if __name__ == "__main__":#當此在主程式可用
    file_name = input("file name: ")
    #file_name="D:\project_spice_new\netlist_put_in_here\\"+file_name
    circuit_name, components, namelist = parse(file_name)
    typelist=[c.type_ for c in components]
    nodelist = [c.nodes for c in components]
    valuelist = [c.value for c in components]
element=len(namelist)
templist=[]
#print(components)
#將nodelist中元素加入templist以方便後續計算node數
#templist:唯一暫存容器
for i in range(element):
    for j in range(2):
        templist.append(nodelist[i][j])
#刪除重複內容 
#print(templist)
templist=list(dict.fromkeys(templist))
#print(templist)
#矩陣大小變數:matrix_size
matrix_size=0
#先算node構成的部分
for i in range(len(templist)):
    if(templist[i]=="0"):
        continue
    elif(templist[i]!=""):
        matrix_size=matrix_size+1
    else:
        continue
#設定一個變數用以儲存目前(僅包含第一類元件)的數量
element2=matrix_size
#print(element2)
#print(matrix_size)
#接著計算第二類元素的部分1
for i in range(element):
    if(namelist[i][0]=="v"):
        matrix_size=matrix_size+1
    elif(namelist[i][0]=="l"):
        matrix_size=matrix_size+1
    else:
        continue
dic_node={}
#用來建立dic引數的一個變數:counter
counter=1
for i in range(len(templist)):
    #nodename等於0就給dic的值為0
    if(templist[i]=='0'):
        dic_node["0"]=0
    else:
        temp_string=templist[i]
        dic_node[temp_string]=counter
        counter=counter+1
#print(dic_node)
print("open circiut:",circuit_name)
#設定模擬指令
cmd=str(input("請輸入模擬內容:"))
vds_his=[]
vgs_his=[]
null=[]
#分割指令
cmd=cmd.split()
#讀取指令
cmd_length=len(cmd)
#print(cmd_length)
if(cmd[0].lower()=="ac"):
#確認起始頻率、步徑、終止頻率
    start_freq=cmd[cmd_length-2].lower()
    end_freq=cmd[cmd_length-1].lower()
    freq_step=cmd[cmd_length-3].lower()
    #將起始頻率、步徑、終止頻率從字串轉成值
    #處理start_freq
    length=len(start_freq)
    start_freq=unit_symbol(start_freq)
    length=len(end_freq)
    end_freq=unit_symbol(end_freq)
    length=len(freq_step)
    freq_step=unit_symbol(freq_step)
    pi=np.pi
    object_source="null"
# =============================================================================
#     null,vds_his,vgs_his=Dc_bias(
#         templist,
#         element2,
#         element,
#         dic_node,
#         matrix_size,
#         cmd,
#         namelist,
#         nodelist,
#         circuit_name,
#         valuelist,
#         typelist,
#         components,
#         object_source,
#         vgs_his,
#         vds_his,
#         null,
#         0
#     )
# =============================================================================
    print(start_freq, " ",end_freq," ",freq_step)
    print("kkk: ",vds_his," ",vgs_his)
    #k=str(input("請輸入模擬內容:"))
    result=[]
    place=[]
    dic_result={}
    times=ac_analysis(
        templist,
        place,
        element2,
        element,
        dic_node,
        matrix_size,
        cmd,
        result,
        namelist,
        nodelist,
        circuit_name,
        valuelist,
        typelist,
        components,
        dic_result,
        freq_step,
        start_freq,
        end_freq,
        vds_his,
        vgs_his
    )
    #times=ac_analysis(templist,element2,element,dic_node,matrix_size,cmd,namelist,nodelist,circuit_name,valuelist,typelist,components,freq_step,start_freq,end_freq)
elif(cmd[0].lower()=="dc"):
    cmd_cond=True
    while cmd_cond==True:
        #讀取指令
        cmd_length=len(cmd)
        #print(cmd_length)
        #確認起始頻率、步徑、終止頻率
        start_V=cmd[cmd_length-3].lower()
        end_V=cmd[cmd_length-2].lower()
        V_step=cmd[cmd_length-1].lower()
        #將起始頻率、步徑、終止頻率從字串轉成值
        #處理start_freq
        start_V=unit_symbol(start_V)
        end_V=unit_symbol(end_V)
        V_step=unit_symbol(V_step)
        result=[]
        times=0
        object_source=cmd[1].lower()
        if(object_source not in namelist):
            print("no source called ",object_source)
            #print(("\n再輸入一次模擬內容:"))
        else:
            cmd_cond=False
            break
    #dc_bias
    C_pass=[[]]
    L_pass=[[]]
    nonlin_his,vds_his,vgs_his,C_pass,L_pass,X0=Dc_bias(
        templist,
        element2,
        element,
        dic_node,
        matrix_size,
        cmd,
        namelist,
        nodelist,
        circuit_name,
        valuelist,
        typelist,
        components,
        object_source,
        vgs_his,
        vds_his,
        null,
        start_V,
        C_pass,
        L_pass
    )
    #sweeping
    
    Dc_analysis(templist,element2,element,dic_node,matrix_size,cmd,namelist,nodelist,circuit_name,valuelist,typelist,components,V_step,start_V,end_V,object_source,vds_his,vgs_his)
else:
    print("開發中!!!!")
#===========================================

