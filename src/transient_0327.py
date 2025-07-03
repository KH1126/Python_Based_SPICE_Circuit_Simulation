from utils import unit_symbol, plot_picture,plot_picture_trans
# from new_parser import parse
from new_parser import parse
from transientfunc_0201 import step_euler,BDF2
from dc_bias import Dc_bias
#from dc_biasing import Dc_bias
import numpy as np
import sympy as sp
from scipy import linalg
from mpmath import mp
import matplotlib.pyplot as plt
#
#讀檔
if __name__ == "__main__":#當此在主程式可用
    file_name = input("file name: ")
    circuit_name, components, namelist = parse(file_name)
    typelist=[c.type_ for c in components]
    nodelist = [c.nodes for c in components]
    valuelist = [c.value for c in components]
    print(circuit_name)
    print(namelist)
    print(nodelist)
    print(valuelist)

#===========================================
#構成矩陣
#step1:找除了0以外的有幾個node
#數有幾個元件，一樣用element當變數
#print("few ",namelist)
element=len(namelist)
templist=[]
timefault=[]
time_neg=[]
revive_step=0
fault=0

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
#儲存非線性元件節點
element1=matrix_size
#設定一個變數用以儲存目前(僅包含第一類元件)的數量
element2=matrix_size
#print(element2)
#print(matrix_size)
#接著計算第二類元素的部分
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
end=0
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

C_total=0
L_total=0
for i in range(len(typelist)):
    if((typelist[i]=='C') ):
        C_total+=1
    elif((typelist[i]=='L') ):
        L_total+=1
    else:
        pass
C_pass=[]
L_pass=[]
for i in range(C_total):
    C_pass.append([])
for i in range(L_total):
    L_pass.append([])
#創建component的dic
dic_component={}
tinv_namelist=[]
tinv_nodelist=[]
tv_namelist=[]
tv_nodelist=[]
nonlin_namelist=[]
nonlin_nodelist=[]
nonlin_his=[]
nonlin_his_check=[]
mos_device_list=[]
mos_device_node=[]
for i in range(len(namelist)):
    if(namelist[i][0].lower()=='r'):
        dic_component.update({namelist[i]:i})
        tinv_namelist.append(namelist[i])
        tinv_nodelist.append(nodelist[i])
    elif(namelist[i][0].lower()!='r' 
         and namelist[i][0].lower()!='d' 
         and namelist[i][0].lower()!='m'):
        dic_component.update({namelist[i]:i})
        tv_namelist.append(namelist[i])
        tv_nodelist.append(nodelist[i])
    elif(namelist[i][0].lower()=='d'):
        nonlin_namelist.append(namelist[i])
        nonlin_nodelist.append(nodelist[i])
        nonlin_his.append([0.0])
    elif(namelist[i][0].lower()=='m'):
        dic_component.update({namelist[i]:i})
        mos_device_list.append(namelist[i])
        mos_device_node.append(nodelist[i])
#print(tinv_namelist)
#print(tv_namelist)
#print(nonlin_namelist)
#print(mos_device_list)
#print(dic_node)
print((templist))
print(len(templist))
while True:
    #設定模擬指令
    cmd=str(input("請輸入模擬內容:"))
    #分割指令
    cmd=cmd.split()
    #讀取指令
    cmd_length=len(cmd)
    #確認起始頻率、步徑、終止頻率
    increment=cmd[cmd_length-2].lower()
    end_time=cmd[cmd_length-1].lower()
    #將起始頻率、步徑、終止頻率從字串轉成值
    #處理start_freq
    increment=unit_symbol(increment)
    end_time=unit_symbol(end_time)
    pi=np.pi
    result=[]
    times=0
    x=[]
    if(cmd[0].lower()=="tran"):
        break
    else :
        print("指令無效")
#填入非時變內容:
#生成TIV線性元件矩陣A
matrixA=np.zeros((matrix_size,matrix_size),dtype=float)
#生成TIV線性元件矩陣B
vectorB=np.zeros(matrix_size,dtype=float)#float maybe danger
for i in range(len(tinv_namelist)):
    #用兩個str容器裝nodename
    str1=tinv_nodelist[i][0]
    str2=tinv_nodelist[i][1]
    #用place去找dic中nodename對應值
    place1=dic_node[str1]
    place2=dic_node[str2]
    comp_number=dic_component[tinv_namelist[i]]
    value=components[comp_number].value
    #print(value)
    if(str1=="0"):
        matrixA[place2-1][place2-1]+=value
    elif(str2=="0"):
        matrixA[place1-1][place1-1]+=(1/value)
    else:
        matrixA[place1-1][place1-1]+=(1/value)
        matrixA[place1-1][place2-1]-=(1/value)
        matrixA[place2-1][place1-1]-=(1/value)
        matrixA[place2-1][place2-1]+=(1/value)
#print(matrixA)
#print(dic_component)
#print(CL_pass)
t=0
sample_point=0
time_step=increment
time_step_new=increment
x=[]
#生成線性矩陣
radiant=2*np.pi*50
offset=0
#print(nonlin_namelist)
#print(nonlin_nodelist)
count=0
C_last=[]
L_last=[]
vgs_his=[]
vds_his=[]
iti_fail=0
object_source="non"
X0=[]
start_V=0
if(len(mos_device_list)!=0):
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
        nonlin_his,
        start_V,
        C_pass,
        L_pass
    )
    print("time ",0,"voltage on cap: ",C_pass)
#print(vds_his,vgs_his,nonlin_his)
#t+=time_step
while (t <= end_time):
    C_pass_temp=C_pass
    L_pass_temp=L_pass
    X=[]
    X_BDF=[]
    checknode=[]
    check_current=[]
    X_diff=[]
    vgs_his_temp=vgs_his
    vds_his_temp=vds_his
    iti_fail=0
    #print(t/end_time*100)
    if(t==0):
        
        #print("time ",t,"voltage on cap: ",C_pass)
        listtest=[L_pass,element1,nonlin_his,matrixA,vectorB,nodelist,C_pass,dic_node
        ,matrix_size,typelist,dic_component,tv_namelist,
        tv_nodelist,element2,components,t,time_step,nonlin_namelist,nonlin_nodelist,
        mos_device_list,mos_device_node,vgs_his,vds_his]
        print(len(listtest))
        for u in range(len(listtest)):
            if listtest[u] is not None:
                print(listtest[u]," ",u)
        X,check_current,checknode,nonlin_his,vgs_his,vds_his,iti_fail=step_euler(
            L_pass,element1,nonlin_his,matrixA,vectorB,nodelist,C_pass,dic_node
            ,matrix_size,typelist,dic_component,tv_namelist,
            tv_nodelist,element2,components,t,time_step,nonlin_namelist,nonlin_nodelist,
            mos_device_list,mos_device_node,vgs_his,vds_his)
        #X=X0
        print(X)
    else:
        #print("III")
        #print("time ",t,"voltage on cap: ",C_pass)
        nonlin_his_check=nonlin_his
        vgs_his_check=vgs_his
        vds_his_check=vds_his
        X,check_current,checknode,nonlin_his,vgs_his,vds_his,iti_fail=step_euler(
            L_pass,element1,nonlin_his,matrixA,vectorB,nodelist,C_pass,dic_node
            ,matrix_size,typelist,dic_component,tv_namelist,
            tv_nodelist,element2,components,t,time_step,nonlin_namelist,nonlin_nodelist,
            mos_device_list,mos_device_node,vgs_his,vds_his)
        if(t<=20*10**-6):
            print(X[-1]," ",t)
        X_BDF=BDF2(L_pass,element1,nonlin_his_check,matrixA,vectorB,nodelist,
                   C_pass,dic_node,matrix_size,typelist,dic_component,
                   tv_namelist,tv_nodelist,element2,components,t,
                   time_step,nonlin_namelist,nonlin_nodelist
                   ,mos_device_list,mos_device_node,vgs_his_check,vds_his_check)
        #X_BDF,check_current,checknode,nonlin_his=BDF2(L_pass,element1,nonlin_his,matrixA,vectorB,nodelist,C_pass,dic_node,matrix_size,typelist,dic_component,tv_namelist,tv_nodelist,element2,components,t,time_step,nonlin_namelist,nonlin_nodelist)
    #迭帶運算
#===============================================
    #存C
    
    for i in range(len(checknode)):
        str1=checknode[i][0]
        str2=checknode[i][1]
        #用place去找dic中nodename對應值
        place1=dic_node[str1]
        place2=dic_node[str2]
        #check_current=[element2-1]
        #填入matrixA
        
        
        if(str1=="0"):
            temp=C_pass[i][-1]
            C_pass[i][0]=temp
            C_pass[i][1]=((-X[place2-1]))
        elif(str2=="0"):
            temp=C_pass[i][-1]
            C_pass[i][0]=temp
            C_pass[i][1]=((X[place1-1]))
        else:
            temp=C_pass[i][-1]
            C_pass[i][0]=temp
            C_pass[i][1]=((X[place1-1]-X[place2-1]))
    #存L
    #print("time ",t,"voltage on cap: ",C_pass)
    
    for i in range(len(check_current)):
    

        L_pass[i][0]=L_pass[i][-1]
        L_pass[i][1]=(X[check_current[i]])
    #print("step: ",time_step)
    #time_step=new_time_step
    print("progress: ",t/end_time*100,"%"," iti times: ",iti_fail)
    false_out=0
    if(t>0):
        i_max_diff=0;
        for u in range(len(X_BDF)):
            #print(X_BDF[u]-X[u])
            X_diff.append(X_BDF[u]-X[u])
            X_diff[-1]=abs(X_diff[-1])
        #print(element2," ",matrix_size)
        temp2=element2
        while(element2<matrix_size):
            #print(element2)
            if(i_max_diff<=abs(X_BDF[element2]-X[element2])):
                i_max_diff=abs(X_BDF[element2]-X[element2])
            element2+=1
        element2=temp2
        #print(X_diff)
        #print(max(X_diff)," ",i_max)
        x_norm=linalg.norm(X)
        i_fail_flag=0
        v_fail_flag=0
        V_max=0
        i_max=0
        temp2=element2
        while(element2<matrix_size):
            #print(element2)
            if(i_max_diff==abs(X_BDF[element2]-X[element2])):
                i_max=max(abs(X_BDF[element2]),abs((X[element2])))
            element2+=1
        element2=temp2
        for u in range(len(X_BDF)):
            if(max(X_diff)==abs(X_BDF[u]-X[u])):
                V_max=max(abs(X_BDF[u]),abs(X_BDF[u]))
        if((((max(X_diff)>V_max*10**-3 and max(X_diff)>10**-6) or i_max_diff>=i_max) and t>0 ) and  iti_fail==101 and time_step/10>10**-15):
            print("fault out",t, " ",time_step," ",max(X_diff)," ",V_max*10**-3)#and iti_fail==100
            t-=time_step
            time_step/=10
            C_pass=C_pass_temp
            L_pass=L_pass_temp
            vgs_his=vgs_his_temp
            vds_his=vds_his_temp
            #new_time_step/=10
            #X=[]
            #X_BDF=[]
            #max(X_diff)<=V_max*10**-3 and i_max_diff<i_max and t>0
        elif( (((max(X_diff)>V_max*10**-3 and max(X_diff)>10**-6)and i_max_diff<i_max and t>0) and iti_fail<=1) and time_step*10<=increment):
            print("revive time step")
            time_step*=10
    #print(t,"success out")
    if(len(X)==0):
        continue
    #畫圖用的x矩陣加入當下頻率值
    #存入上一個值
    
        #print("L ",L_pass[-1])
    #解的矩陣:result
    #t=round(t,4)
    if(t>=sample_point or t==end_time):
        
        result.append(X)
        x.append(t)
        sample_point+=increment
    else:
        pass
    count+=1
    if(t+time_step>end_time and end==0):
        t=end_time
        end=1
    elif(end==1):
        t+=time_step
    else:
        t+=time_step
#print(len(x))
#將結點名稱對上result矩陣
#np.savetxt('out.txt',result,delimiter=',')
# =============================================================================
#畫圖function define in utils
#print(time_neg)
#print(timefault[0])
#print(len(timefault))
#print(result[-1][1])
path='rea.txt'
print(len(result))
f=open(path,"w")

for i in range(len(result)):
    #print("pop")
    string=str(x[i])+" "+str(result[i][1]-result[i][2])+"\n"
    f.writelines(string)
    #f.write("")
f.close()
plot_picture_trans(templist,namelist,result,cmd,x,circuit_name)
print("~bye")
