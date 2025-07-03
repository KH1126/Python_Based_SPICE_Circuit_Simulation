#被我改過了喔
#transient functions
from utils import unit_symbol, plot_picture,plot_picture_trans,plot_picture_dc
import sympy as sp
from mpmath import mp
import numpy as np
import matplotlib.pyplot as plt
import re
from scipy import linalg
def step_euler(
    L_pass,element1,nonlin_his,matrixA,vectorB,nodelist,C_pass,dic_node
    ,matrix_size,typelist,dic_component,tv_namelist,
    tv_nodelist,element2,components,t,time_step,nonlin_namelist,nonlin_nodelist,
    mos_device_list,mos_device_node,vgs_his,vds_his):    
    checknode=[]
    check_current=[]
    C_pointer=0
    L_pointer=0
    pi=np.pi
    radiant=2*np.pi*50
    #mos====================
    Vt = 0.4
    # W_L = 1/0.18
    iti_report=0
    Lambda = 0.02  # 為1/VA
    va_n = 1 / Lambda
    va_p = 1 / Lambda
    u_n = 600e-4  # ngspice預設(600cm^2/V sec)
    u_p=600e-4
    t_ox = 10**-7  # ngspice預設
    epsilon_0 = 8.854e-12  # 真空介電常數
    epsilon_SiO2 = 3.9  # SiO2相對介電常數
    C_ox = epsilon_0 * epsilon_SiO2 / t_ox
    k_n = u_n * C_ox
    k_p=u_p * C_ox
    #mos====================
    #生成TV線性元件矩陣A
    matrixA_TV=np.zeros((matrix_size,matrix_size),dtype=float)
    #生成TV線性元件矩陣B
    vectorB_TV=np.zeros(matrix_size,dtype=float)#float maybe danger
    #生成Non線性元件矩陣A
    matrixA_Non=np.zeros((matrix_size,matrix_size),dtype=float)
    #生成Non線性元件矩陣B
    vectorB_Non=np.zeros(matrix_size,dtype=float)#float maybe danger
    #生成Mos元件矩陣B
    vectorB_mos = np.zeros(matrix_size)
    #生成Mos元件矩陣A
    matrixA_mos = np.zeros((matrix_size, matrix_size))
    temp=element2
    for i in range(len(tv_namelist)):
        comp_number=dic_component[tv_namelist[i]]
        ##print("KKK")
        #填入電阻的stamp
        #填入矩陣A的部分
        #填入電容的stamp
        #填入矩陣A的部分
        if(typelist[comp_number]=='C'):
           ##print("KKR")
            #用兩個str容器裝nodename
            str1=tv_nodelist[i][0]
            str2=tv_nodelist[i][1]
            #用place去找dic中nodename對應值
            place1=dic_node[str1]
            place2=dic_node[str2]
            comp_number=dic_component[tv_namelist[i]]
            value=components[comp_number].value
            checknode.append([str1,str2])
            if(str1=="0"):
                matrixA_TV[place2-1][place2-1]+=value/time_step
            elif(str2=="0"):
                matrixA_TV[place1-1][place1-1]+=value/time_step
            else:
                matrixA_TV[place1-1][place1-1]+=value/time_step
                matrixA_TV[place1-1][place2-1]-=value/time_step
                matrixA_TV[place2-1][place1-1]-=value/time_step
                matrixA_TV[place2-1][place2-1]+=value/time_step
            if(str1=="0"):
                vectorB_TV[place2-1]-=C_pass[C_pointer][-1]*value/time_step
            elif(str2=="0"):
                #print(C_pass[C_pointer])
                vectorB_TV[place1-1]+=C_pass[C_pointer][-1]*value/time_step
            else:
                vectorB_TV[place1-1]+=C_pass[C_pointer][-1]*value/time_step
                vectorB_TV[place2-1]-=C_pass[C_pointer][-1]*value/time_step
            C_pointer+=1
        #填入電流源的stamp
        #無須填入矩陣A的部分
        #僅填入vectorB
        elif(typelist[comp_number]=='I'):#暫時先不動，之後要改
            #用兩個str容器裝nodename
            str1=nodelist[i][0]
            str2=nodelist[i][1]
            #用place去找dic中nodename對應值
            place1=dic_node[str1]
            place2=dic_node[str2]
            if(components[i].sin!=0):
                if(str1=="0"):
                    vectorB[place2-1]+=components[i].sin*np.sin(t*radiant)
                elif(str2=="0"):
                    vectorB[place1-1]-=components[i].sin*np.sin(t*radiant)  
                else:
                    vectorB[place2-1]+=components[i].sin*np.sin(t*radiant)
                    vectorB[place1-1]-=components[i].sin*np.sin(t*radiant)
            else:
                if(str1=="0"):
                    vectorB[place2-1]+=components[i].dc
                elif(str2=="0"):
                    vectorB[place1-1]-=components[i].dc   
                else:
                    vectorB[place2-1]+=components[i].dc
                    vectorB[place1-1]-=components[i].dc
        #填入電壓源的stamp
        #須填入矩陣A的部分
        #同時填入vectorB
        elif(typelist[comp_number]=='V'):
            element2=element2+1#作為一個指標用來定出元件所放置的位置
            #用兩個str容器裝nodename
            str1=tv_nodelist[i][0]
            str2=tv_nodelist[i][1]
            #用place去找dic中nodename對應值
            place1=dic_node[str1]
            place2=dic_node[str2]
            comp_number=dic_component[tv_namelist[i]]
            #填入vectorB
            if(components[comp_number].sin[1]!=0):
                #print("EUer ",components[comp_number].sin," ",components[comp_number])
                voff=components[comp_number].sin[0]
                amp=components[comp_number].sin[1]
                rads=components[comp_number].sin[2]*2*pi
                delay=components[comp_number].sin[3]
                damp=components[comp_number].sin[4]
                phase=components[comp_number].sin[5]
                vectorB_TV[element2-1]+=amp*np.sin(t*rads+phase)+voff
                #填入matrixA_TV
                if(str1=="0"):
                    matrixA_TV[element2-1][place2-1]-=1
                    matrixA_TV[place2-1][element2-1]-=1
                elif(str2=="0"):
                    matrixA_TV[element2-1][place1-1]+=1
                    matrixA_TV[place1-1][element2-1]+=1  
                else:
                    matrixA_TV[element2-1][place1-1]+=1
                    matrixA_TV[element2-1][place2-1]-=1
                    matrixA_TV[place1-1][element2-1]+=1
                    matrixA_TV[place2-1][element2-1]-=1
            elif(components[comp_number].pulse[6]!=0):
                #print("EUer ",components[comp_number].sin," ",components[comp_number])
                V1=components[comp_number].pulse[0]
                V2=components[comp_number].pulse[1]
                TD=components[comp_number].pulse[2]
                TR=components[comp_number].pulse[3]
                TF=components[comp_number].pulse[4]
                PW=components[comp_number].pulse[5]
                PER=components[comp_number].pulse[6]
                NP=components[comp_number].pulse[7]
                t_place=t-TD
                T_in_PER=t_place%PER
                if(t<=TD or T_in_PER==0):
                    vectorB_TV[element2-1]+=V1
                elif(t>TD and T_in_PER>0 and T_in_PER<(TR)):
                    vectorB_TV[element2-1]+=V1+(V2-V1)*(T_in_PER)/TR
                elif(t>TD and T_in_PER>=TR and T_in_PER<(TR+PW)):
                    vectorB_TV[element2-1]+=V2
                elif(t>TD and T_in_PER>=TR+PW and T_in_PER<(TR+PW+TF)):
                    vectorB_TV[element2-1]+=V2+(V1-V2)/TF*(T_in_PER-TR-PW)
                elif(t>TD and T_in_PER>=TR+PW+TF and T_in_PER<PER):
                    vectorB_TV[element2-1]+=V1
                
                #填入matrixA_TV
                if(str1=="0"):
                    matrixA_TV[element2-1][place2-1]-=1
                    matrixA_TV[place2-1][element2-1]-=1
                elif(str2=="0"):
                    matrixA_TV[element2-1][place1-1]+=1
                    matrixA_TV[place1-1][element2-1]+=1  
                else:
                    matrixA_TV[element2-1][place1-1]+=1
                    matrixA_TV[element2-1][place2-1]-=1
                    matrixA_TV[place1-1][element2-1]+=1
                    matrixA_TV[place2-1][element2-1]-=1
            
            elif(len(components[comp_number].pwl)>0):
                time_table=[]
                V_table=[]
                for u in range(len(components[comp_number].pwl)):
                        V_table.append(components[comp_number].pwl[u][1])
                        
                        time_table.append(components[comp_number].pwl[u][0])
                #填入vector TV
                #print(V_table)
                #print(time_table)
                voltage_ptr=0
                for u in range(len(time_table)):
                    if(t<time_table[u]):
                        voltage_ptr=u-1#指向第n個電壓
                        break
                if(t>=time_table[-1]):
                    vectorB_TV[element2-1]+=V_table[-1]
                elif(V_table[voltage_ptr]==V_table[voltage_ptr+1]):
                    #print(V_table[voltage_ptr+1])
                    vectorB_TV[element2-1]+=V_table[voltage_ptr]
                
                else:
                    voltage_gap=(V_table[voltage_ptr+1]-V_table[voltage_ptr])
                    time_gap=time_table[voltage_ptr+1]-time_table[voltage_ptr]
                    voltage_in_variation=V_table[voltage_ptr]+(t-time_table[voltage_ptr])*voltage_gap/time_gap
                    print(voltage_in_variation)
                    vectorB_TV[element2-1]+=voltage_in_variation
                #填入matrixA_TV
                if(str1=="0"):
                    matrixA_TV[element2-1][place2-1]-=1
                    matrixA_TV[place2-1][element2-1]-=1
                elif(str2=="0"):
                    matrixA_TV[element2-1][place1-1]+=1
                    matrixA_TV[place1-1][element2-1]+=1  
                else:
                    matrixA_TV[element2-1][place1-1]+=1
                    matrixA_TV[element2-1][place2-1]-=1
                    matrixA_TV[place1-1][element2-1]+=1
                    matrixA_TV[place2-1][element2-1]-=1
            else:#需要修改，給DC電壓源在外面填
                vectorB[element2-1]+=components[comp_number].dc
                #print("ewfew ",t)
            #填入matrixA
            
                if(str1=="0"):
                    matrixA[element2-1][place2-1]-=1
                    matrixA[place2-1][element2-1]-=1
                elif(str2=="0"):
                    matrixA[element2-1][place1-1]+=1
                    matrixA[place1-1][element2-1]+=1  
                else:
                    matrixA[element2-1][place1-1]+=1
                    matrixA[element2-1][place2-1]-=1
                    matrixA[place1-1][element2-1]+=1
                    matrixA[place2-1][element2-1]-=1
        #填入電感的stamp
        #須填入矩陣A的部分
        #無須填入vectorB
        elif(typelist[comp_number]=='L'):
            element2=element2+1#作為一個指標用來定出元件所放置的位置
            #用兩個str容器裝nodename
            str1=tv_nodelist[i][0]
            str2=tv_nodelist[i][1]
            #用place去找dic中nodename對應值
            place1=dic_node[str1]
            place2=dic_node[str2]
            comp_number=dic_component[tv_namelist[i]]
            value=components[comp_number].value
            check_current.append(element2-1)
            #填入matrixA
            if(str1=="0"):
                matrixA_TV[element2-1][place2-1]-=1
                matrixA_TV[place2-1][element2-1]-=1
                matrixA_TV[element2-1][element2-1]-=value/time_step
                vectorB_TV[element2-1]-=L_pass[L_pointer][-1]*value/time_step
            elif(str2=="0"):
                matrixA_TV[element2-1][place1-1]+=1
                matrixA_TV[place1-1][element2-1]+=1 
                matrixA_TV[element2-1][element2-1]-=value/time_step
                vectorB_TV[element2-1]-=L_pass[L_pointer][-1]*value/time_step
            else:
                matrixA_TV[element2-1][place1-1]+=1
                matrixA_TV[element2-1][place2-1]-=1
                matrixA_TV[place1-1][element2-1]+=1
                matrixA_TV[place2-1][element2-1]-=1
                matrixA_TV[element2-1][element2-1]-=value/time_step
                vectorB_TV[element2-1]-=L_pass[L_pointer][-1]*value/time_step
            L_pointer+=1
        else:
            continue
    #填入
    #LU解法==========================================   
    #fill in mos
    if(len(nonlin_namelist)==0 and len(mos_device_list)==0):
        matrixA_mask=np.zeros((matrix_size,matrix_size),dtype=float)
        for s in range(element1):
            matrixA_mask[s][s]+=10**-12
        
        matrixA_mix=matrixA_TV+matrixA_Non+matrixA+matrixA_mask
        vectorB_mix=vectorB+vectorB_Non+vectorB_TV
        A=np.array(matrixA_mix,dtype=float)
        P, L, U = linalg.lu(A)
        B=np.array(vectorB_mix)
        Y=linalg.solve_triangular(L,P.T@B,lower=True)
        X=linalg.solve_triangular(U,Y,lower=False)
    else:
        cond=True
        iti=0
        #cap=0
        while(cond):
            #cond_false=0
            modified=0
            for j in range(len(nonlin_namelist)):
                ##print(t," ",j)
                ##print(j)
                ##print(x_0[diode_number])
                #用兩個str容器裝nodename
                str1=nonlin_nodelist[j][0]
                str2=nonlin_nodelist[j][1]
                #用place去找dic中nodename對應值
                place1=dic_node[str1]
                place2=dic_node[str2]
                ##print((nonlin_his[j][-1]/0.0258528413-32.2361913))
                #a=np.array(nonlin_his[j][-1]/0.0258528413-32.2361913, dtype=float)
                a=nonlin_his[j][-1]/0.0258528413-32.2361913
                tempa=np.exp(a)/0.0258528413
                tempb=np.exp(a)-10**(-14)-nonlin_his[j][-1]*np.exp(a)/0.0258528413
                #填入matrixA
                if(str1=="0"):
                    matrixA_Non[place2-1][place2-1]+=tempa
                elif(str2=="0"):
                    matrixA_Non[place1-1][place1-1]+=tempa
                else:
                    matrixA_Non[place1-1][place1-1]+=tempa
                    matrixA_Non[place1-1][place2-1]-=tempa
                    matrixA_Non[place2-1][place1-1]-=tempa
                    matrixA_Non[place2-1][place2-1]+=tempa
                #填入vectorB
                if(str1=="0"):
                    vectorB_Non[place2-1]+=tempb
                elif(str2=="0"):
                    vectorB_Non[place1-1]-=tempb
                else:
                    vectorB_Non[place2-1]+=tempb
                    vectorB_Non[place1-1]-=tempb
            #fill in mos
            for i in range(len(mos_device_list)):
                comp_number=dic_component[mos_device_list[i]]
                place1 = dic_node[mos_device_node[i][0]]
                place2 = dic_node[mos_device_node[i][1]]
                place3 = dic_node[mos_device_node[i][2]]
                place4 = dic_node[mos_device_node[i][3]]
                D = mos_device_node[i][0]
                G = mos_device_node[i][1]
                S = mos_device_node[i][2]
                B = mos_device_node[i][3]
                a11 = 0
                a12 = 0
                b1 = 0

                if mos_device_list[i][1] == "p":
        # =============================================================================
        #                 Vgs = vgs_his[u]
        #                 Vds = vds_his[u]
        # =============================================================================
                    Vgs = vgs_his[i][-1]
                    Vds = vds_his[i][-1]
                    ##print(vgs_his," ",vds_his)
                    #mos_ptr += 1
                    width = unit_symbol(components[comp_number].args["w"])
                    length = unit_symbol(components[comp_number].args["l"])
                    W_L = width / length
                    
                    # Cutoff
                    if Vgs < Vt:
                        
                        a12 = 0
                        a11 = sp.exp((Vgs - Vt) / 0.0258528413 - 32.2361913) / 0.0258528413
                        b1 = -(a11 * 0.0258528413 - a11 * (Vgs))
                    
                    # Triode
                    elif (Vgs - Vt) > Vds:

                        a11 = W_L * k_p * Vds
                        a12 = W_L * k_p * ((Vgs - Vt) - Vds)
                        b1 = -(
                            W_L * k_p * ((Vgs - Vt) * Vds - (Vds) ** 2 / 2)
                            - a11 * Vgs
                            - a12 * Vds
                        )

                    # Sat
                    elif (Vgs - Vt) <= Vds:
                        a11 = (
                            (1 / 2)
                            * W_L
                            * k_p
                            * (
                                2 * (Vgs - Vt) * (1 + 1 / va_p * (Vds - (Vgs - Vt)))
                                - 1 / va_p * (Vgs - Vt) ** 2
                            )
                        )
                        a12 = (1 / 2) * W_L * k_p * 1 / va_p * (Vgs - Vt) ** 2
                        b1 =-(
                            (1 / 2)
                            * W_L
                            * k_p
                            * (Vgs - Vt) ** 2
                            * (1 + 1 / va_p * (Vds - (Vgs - Vt)))
                            - a11 * Vgs
                            - a12 * Vds
                        )
                        
                elif mos_device_list[i][1] == "n":
                    
                    ##print(vgs_his)
                    Vgs = vgs_his[i][-1]
                    Vds = vds_his[i][-1]
                    ##print(vgs_his," ",vds_his)
                    #mos_ptr += 1
                    width = unit_symbol(components[comp_number].args["w"])
                    length = unit_symbol(components[comp_number].args["l"])
                    W_L = width / length
                    ##print(width," ",length," n",)
                    
                    # Cutoff
                    if Vgs < Vt:
                        
                        ##print("cutoff")
                        a12 = 0
                        a11 = sp.exp((Vgs - Vt) / 0.0258528413 - 32.2361913) / 0.0258528413
                        b1 = a11 * 0.0258528413 - a11 * (Vgs)
                    
                    # Triode
                    elif (Vgs - Vt) > Vds:
                        ##print("triode")
                        
                        a11 = W_L * k_n * Vds
                        a12 = W_L * k_n * ((Vgs - Vt) - Vds)
                        b1 = (
                            W_L * k_n * ((Vgs - Vt) * Vds - (Vds) ** 2 / 2)
                            - a11 * Vgs
                            - a12 * Vds
                        )
                    
                    # Sat
                    elif (Vgs - Vt) <= Vds:
                        ##print("sat")
                        
                        a11 = (
                            (1 / 2)
                            * W_L
                            * k_n
                            * (
                                2 * (Vgs - Vt) * (1 + 1 / va_n * (Vds - (Vgs - Vt)))
                                - 1 / va_n * (Vgs - Vt) ** 2
                            )
                        )
                        a12 = (1 / 2) * W_L * k_n * 1 / va_n * (Vgs - Vt) ** 2
                        b1 = (
                            (1 / 2)
                            * W_L
                            * k_n
                            * (Vgs - Vt) ** 2
                            * (1 + 1 / va_n * (Vds - (Vgs - Vt)))
                            - a11 * Vgs
                            - a12 * Vds
                        )
                ##print(a11," ",a12)
                if D == "0" and G == "0" and S == "0":
                    
                    pass
                elif D == "0":
                    if G == "0" and S != "0":
                        matrixA_mos[place3 - 1][place3 - 1] += a11 + a12
                        vectorB_mos[place3 - 1] += b1
                    elif S == "0" and G != "0":
                        pass
                    elif S != "0" and G != "0":
                        matrixA_mos[place3 - 1][place1 - 1] -= a11
                        matrixA_mos[place3 - 1][place3 - 1] += a11 + a12
                        vectorB_mos[place3 - 1] += b1
                elif G == "0":
                    if D == "0" and S != "0":
                        matrixA_mos[place3 - 1][place3 - 1] += a11 + a12
                        vectorB_mos[place3 - 1] += b1
                    elif S == "0" and D != "0":
                        matrixA_mos[place1 - 1][place1 - 1] += a12
                        vectorB_mos[place1 - 1] -= b1
                    elif S != "0" and D != "0":
                        matrixA_mos[place1 - 1][place1 - 1] += a12
                        matrixA_mos[place1 - 1][place3 - 1] -= a11 + a12
                        vectorB_mos[place1 - 1] -= b1
                        matrixA_mos[place3 - 1][place1 - 1] -= a12
                        matrixA_mos[place3 - 1][place3 - 1] += a11 + a12
                        vectorB_mos[place3 - 1] += b1
                elif S == "0":
                    if D == "0" and G != "0":
                        pass
                    elif G == "0" and D != "0":
                        matrixA_mos[place1 - 1][place1 - 1] += a12  # Change
                        vectorB_mos[place1 - 1] -= b1
                    elif G != "0" and D != "0":
                        matrixA_mos[place1 - 1][place2 - 1] += a11  # Change
                        matrixA_mos[place1 - 1][place1 - 1] += a12  # Change  
                        vectorB_mos[place1 - 1] -= b1
                else:
                    matrixA_mos[place1 - 1][place2 - 1] += a11
                    matrixA_mos[place1 - 1][place1 - 1] += a12
                    matrixA_mos[place1 - 1][place3 - 1] -= (a11 + a12)
                    vectorB_mos[place1 - 1] -= b1
                    matrixA_mos[place3 - 1][place2 - 1] -= a11
                    matrixA_mos[place3 - 1][place1 - 1] -= a12
                    matrixA_mos[place3 - 1][place3 - 1] += (a11 + a12)
                    vectorB_mos[place3 - 1] += b1  
                    #fill in mos end========================        
            matrixA_mask=np.zeros((matrix_size,matrix_size),dtype=float)
            for s in range(element1):
                matrixA_mask[s][s]+=10**-12
            matrixA_mix=matrixA_TV+matrixA_Non+matrixA+matrixA_mask+matrixA_mos
            vectorB_mix=vectorB+vectorB_Non+vectorB_TV+vectorB_mos
            A=np.array(matrixA_mix,dtype=float)
# =============================================================================
#             if(t==0):
#                 print("matrixA",A)
# =============================================================================
            
            condition_number=np.linalg.cond(A)
            ##print("cond: ",condition_number)
# =============================================================================
#             if((condition_number>=10**16 and iti==0) or time_step_new>increment ):
#                 cond_false=1
#                 if(time_step_new>10**-5):
#                     time_step_new/=10
#                 else:
#                     time_step_new=time_step_new
#                 #print("wrong",time_step_new)
#             elif(condition_number<10**8 and time_step_new<increment and iti==0):
#                 time_step_new*=10
#                 if(time_step_new>=increment):
#                     time_step_new=increment
#                 #print("speed_up",time_step_new)
# =============================================================================
            #elif(condition_number<10**8):
                #time_step*=10
            ##print(iti)
            P, L, U = linalg.lu(A)
            B=np.array(vectorB_mix)
            Y=linalg.solve_triangular(L,P.T@B,lower=True)
            X_list=linalg.solve_triangular(U,Y,lower=False)
            X=X_list
            ##print(X_list)
            #更新vgs vds=============================
            for u in range(len(mos_device_list)):
                comp_number=dic_component[mos_device_list[u]]
                place1 = dic_node[mos_device_node[u][0]]
                place2 = dic_node[mos_device_node[u][1]]
                place3 = dic_node[mos_device_node[u][2]]
                place4 = dic_node[mos_device_node[u][3]]
                D = mos_device_node[u][0]
                G = mos_device_node[u][1]
                S = mos_device_node[u][2]
                B = mos_device_node[u][3]
                #print(mos_device_list[u]," D: ",D," G:",G," S: ",S," B: ",B," ",t)
                if D == "0" and G == "0" and S == "0":
                        pass
                    
                elif D == "0":
                    if G == "0" and S != "0":
                        if mos_device_list[u][1] == "p":
                            vds_his[u][0]=vds_his[u][1]
                            vds_his[u][1] = (X[place3 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (X[place3 - 1])
                        else:
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = (-X[place3 - 1])
                            vgs_his[u][0]=vgs_his[u][1]
                            vgs_his[u][1] = (-X[place3 - 1])
                    elif S == "0" and G != "0":
                        if mos_device_list[u][1] == "p":
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1]=0
                            vgs_his[u][0] =vgs_his[u][1]
                            vgs_his[u][1] = (-X[place2 - 1])
                        else:
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = (0)
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (X[place2 - 1])
                    elif S != "0" and G != "0":
                        if mos_device_list[u][1] == "p":
                            vds_his[u][0] =vds_his[u][1]
                            vds_his[u][1] = (X[place3 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (
                                X[place3 - 1] - X[place2 - 1]
                            )
                        else:
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = (-X[place3 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (
                                X[place2 - 1] - X[place3 - 1]
                            )
                            
                elif G == "0":
                    if D == "0" and S != "0":
                        if mos_device_list[u][1] == "p":
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = (X[place3 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (X[place3 - 1])
                        else:
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = (-X[place3 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (-X[place3 - 1])
                    elif S == "0" and D != "0":
                        if mos_device_list[u][1] == "p":
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = (-X[place1 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (0)
                        else:
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = (X[place1 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (0)
                    elif S != "0" and D != "0":
                        if mos_device_list[u][1] == "p":
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = (
                                X[place3 - 1] - X[place1 - 1]
                            )
                            vgs_his[u][0] = vgs_his[u][1] 
                            vgs_his[u][1] = (X[place3 - 1])
                        else:
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = (
                                X[place1 - 1] - X[place3 - 1]
                            )
                            vgs_his[u][0] = vgs_his[u][1] 
                            vgs_his[u][1] = (-X[place3 - 1])
                            
                elif S == "0":
                    if D == "0" and G != "0":
                        if mos_device_list[u][1] == "p":
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = 0
                            vgs_his[u][0] = vgs_his[u][1] 
                            vgs_his[u][1] = (-X[place2 - 1])
                        else:
                            vds_his[u][0] =vds_his[u][1] 
                            vds_his[u][1] = (0)
                            vgs_his[u][0] =vgs_his[u][1] 
                            vgs_his[u][1] = (X[place2 - 1] )
                    elif G == "0" and D != "0":
                        if mos_device_list[u][1] == "p":
                            vds_his[u][0] =vds_his[u][1] 
                            vds_his[u][1] = (-X[place1 - 1])
                            vgs_his[u][0]=vgs_his[u][1]
                            vgs_his[u][1] = (0)
                        else:
                            vds_his[u][0] =vds_his[u][1] 
                            vds_his[u][1] = (X[place1 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (0)
                    elif G != "0" and D != "0":
                        if mos_device_list[u][1] == "p":
                            vds_his[u][0] = vds_his[u][1] 
                            vds_his[u][1] = (-X[place1 - 1])
                            vgs_his[u][0] = vgs_his[u][1] 
                            vgs_his[u][1] = (-X[place2 - 1])
                        else:
                            vds_his[u][0] = vds_his[u][1] 
                            vds_his[u][1] = (X[place1 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (X[place2 - 1])
                            
                else:
                    if mos_device_list[u][1] == "p":
                        vds_his[u][0] = vds_his[u][1] 
                        vds_his[u][1] = (
                            X[place3 - 1] - X[place1 - 1]
                        )
                        vgs_his[u][0] = vgs_his[u][1]
                        vgs_his[u][1] = (
                            X[place3 - 1] - X[place2 - 1]
                        )
                        ##print("check1: ",X[place3 - 1] - X[place1 - 1]," ",X[place3 - 1] - X[place2 - 1])
                    else:
                        # #print(X[place1-1]-X[place3-1])
                        vds_his[u][0] = vds_his[u][1]
                        vds_his[u][1] = (
                            X[place1 - 1] - X[place3 - 1]
                        )
                        vgs_his[u][0] = vgs_his[u][1]
                        vgs_his[u][1] = (
                            X[place2 - 1] - X[place3 - 1]
                        )
            #更新vgs vds END=====================
            #確認跨壓是否過大並修正
            for u in range(len(nonlin_nodelist)):
                ##print(x_0[diode_number])
                #用兩個str容器裝nodename
                str1=nonlin_nodelist[u][0]
                str2=nonlin_nodelist[u][1]
                #用place去找dic中nodename對應值
                place1=dic_node[str1]
                place2=dic_node[str2]
                
                #先填入解出值
                if(str1=="0"):
                    temp1=nonlin_his[u][-1]
                    temp2=-X_list[place2-1]
                    nonlin_his[u]=[temp1,temp2]
                elif(str2=="0"):
                    temp1=nonlin_his[u][-1]
                    temp2=X_list[place1-1]
                    nonlin_his[u]=[temp1,temp2]
                else:
                    temp1=nonlin_his[u][-1]
                    temp2=X_list[place1-1]-X_list[place2-1]
                    nonlin_his[u]=[temp1,temp2]
                ##print(nonlin_his)
                #檢查過大並修正
                if(nonlin_his[u][-1]>0.734):#更正掉過大的x_0
                    ##print(nonlin_his[u][-1])
                    #if(iti==0):
                        #time_step_new/=1000
                        ##print("wrong",time_step_new)
                        ##print(nonlin_his[u][-1])
                    #timefault.append([t,con])
                    x_klast=nonlin_his[u][-2]
                    ##print("fwve:",x_klast)
                    a=x_klast/0.0258528413-32.2361913
                    b=x_klast/0.0258528413
                    #a=np.array( a, dtype=float)
                    #b=np.array( b, dtype=float)
                    #a=x_klast/0.0258528413-32.2361913
                    #b=x_klast/0.0258528413
                    i_temp=nonlin_his[u][-1]*np.exp(a)/0.0258528413+np.exp(a)-10**(-14)-np.exp(a)*x_klast/0.0258528413
                    ##print(i_temp+10**(-14))
                    if(i_temp<0):
                        i_temp=10**-14
                    #tempk=1+i_temp/10**(-14)
                    queew=np.log(i_temp)-np.log(10**-14)
                    tempk=1+np.exp(queew)
                    ##print(tempk)
                    #tempk=np.array( tempk, dtype=float)
                    if(tempk<0):
                        ##print(condition_number)
                        tempk=-tempk
                        nonlin_his[u][-1]=0.0258528413*np.log(tempk)
                        ##print("a",0.0258528413*np.log(tempk))
                    elif(tempk==0.0):
                        ##print(condition_number)
                        #timefault.append(t)
                        tempk=2
                        nonlin_his[u][-1]=0.732
                        ##print("b",0.0258528413*np.log(tempk))
                    else:
                        nonlin_his[u][-1]=0.0258528413*np.log(tempk)
                       
                        ##print("c",0.0258528413*np.log(tempk))
                    modified=1
            
            ##print(nonlin_his)
            non_diff=[]
            for u in range(len(nonlin_nodelist)):
                non_diff.append(nonlin_his[u][-1]-nonlin_his[u][-2])
            if(len(non_diff)==0):
                non_diff=[0]
            ##print(non_diff)
            vgs_diff=[]
            vds_diff=[]
            for u in range (len(vgs_his)):
                vgs_diff.append(abs(vgs_his[u][1]-vgs_his[u][0]))
                vds_diff.append(abs(vds_his[u][1]-vds_his[u][0]))
            if(len(vgs_diff)==0):
                vgs_diff=[0]
                vds_diff=[0]
            ptr_max_vgs=0
            ptr_max_vds=0
            for u in range(len(vgs_diff)):
                if(max(vgs_diff)==vgs_diff[u]):
                    ptr_max_vgs=u
                if(max(vds_diff)==vds_diff[u]):
                    ptr_max_vds=u
            v_max_vgs=0
            v_max_vds=0
            if(len(vgs_his)!=0):
                v_max_vgs=max(abs(vgs_his[ptr_max_vgs][0]),abs(vgs_his[ptr_max_vgs][1]))
                v_max_vds=max(abs(vds_his[ptr_max_vds][0]),abs(vds_his[ptr_max_vds][1]))
                
            #print(max(vgs_diff)," x ",max(vds_diff)," ",iti)
            if((iti==0 or modified==1 or max(non_diff)>10**(-6) or 
                abs(min(non_diff))>10**(-6)
                or max(vgs_diff)> v_max_vgs*10**-3 or max(vds_diff)>v_max_vds*10**-3)
               and iti<=100):
                
                iti+=1
                iti_report=iti
                #生成Non線性元件矩陣A
                matrixA_Non=np.zeros((matrix_size,matrix_size),dtype=float)
                #生成Non線性元件矩陣B
                vectorB_Non=np.zeros(matrix_size,dtype=float)
                #生成Mos元件矩陣B
                vectorB_mos = np.zeros(matrix_size)
                #生成Mos元件矩陣A
                matrixA_mos = np.zeros((matrix_size, matrix_size))
                modified=0
                pass
            else:
                ##print("pass ",t/end_time*100,"%")
                for u in range(len(nonlin_nodelist)):
                    nonlin_his[u]=[nonlin_his[u][-1]]
                for u in range(len(mos_device_list)):
                    vgs_his[u]=vgs_his[u]
                    vds_his[u]=vds_his[u]
# =============================================================================
#                 if(iti<2):
#                     revive_step+=1
#                 if(iti>=2 and cond_false==0):
#                     revive_step=0
#                     if(time_step_new>10**-5):
#                         time_step_new/=1000
#                     else:
#                         time_step_new=time_step_new
#                     #print("wrong",time_step_new)
#                 elif(revive_step==1 and time_step_new<increment and cond_false==0):
#                     revive_step=0
#                     time_step_new*=1000
#                     if(time_step_new>=increment):
#                         time_step_new=increment
# =============================================================================
                ##print("iti over",iti)
                
                iti=0
                cond=False
                modified=0
                return (X_list),check_current,checknode,nonlin_his,vgs_his,vds_his,iti_report
def BDF2(L_pass,element1,nonlin_his,matrixA,vectorB,nodelist,C_pass,
         dic_node,matrix_size,typelist,dic_component,tv_namelist,tv_nodelist,
         element2,components,t,time_step,nonlin_namelist,nonlin_nodelist,
         mos_device_list,mos_device_node,vgs_his,vds_his):
    checknode=[]
    check_current=[]
    C_pointer=0
    L_pointer=0
    pi=np.pi
    radiant=2*np.pi*50
    #mos====================
    Vt = 0.4
    # W_L = 1/0.18
    Lambda = 0.02  # 為1/VA
    va_n = 1 / Lambda
    va_p = 1 / Lambda
    u_n = 600e-4  # ngspice預設(600cm^2/V sec)
    u_p=600e-4
    t_ox = 10**-7  # ngspice預設
    epsilon_0 = 8.854e-12  # 真空介電常數
    epsilon_SiO2 = 3.9  # SiO2相對介電常數
    C_ox = epsilon_0 * epsilon_SiO2 / t_ox
    k_n = u_n * C_ox
    k_p=u_p * C_ox
    #mos====================
    ##print(t)
    ##print(t/end_time,"%")
    #生成TV線性元件矩陣A
    matrixA_TV=np.zeros((matrix_size,matrix_size),dtype=float)
    #生成TV線性元件矩陣B
    vectorB_TV=np.zeros(matrix_size,dtype=float)#float maybe danger
    #生成Non線性元件矩陣A
    matrixA_Non=np.zeros((matrix_size,matrix_size),dtype=float)
    #生成Non線性元件矩陣B
    vectorB_Non=np.zeros(matrix_size,dtype=float)#float maybe danger
    #生成Mos元件矩陣B
    vectorB_mos = np.zeros(matrix_size)
    #生成Mos元件矩陣A
    matrixA_mos = np.zeros((matrix_size, matrix_size))
    temp=element2
    for i in range(len(tv_namelist)):
        comp_number=dic_component[tv_namelist[i]]
        ##print("KKK")
        #填入電阻的stamp
        #填入矩陣A的部分
        #填入電容的stamp
        #填入矩陣A的部分
        if(typelist[comp_number]=='C'):
           ##print("KKR")
            #用兩個str容器裝nodename
            str1=tv_nodelist[i][0]
            str2=tv_nodelist[i][1]
            #用place去找dic中nodename對應值
            place1=dic_node[str1]
            place2=dic_node[str2]
            comp_number=dic_component[tv_namelist[i]]
            value=components[comp_number].value
            checknode.append([str1,str2])
            if(str1=="0"):
                matrixA_TV[place2-1][place2-1]+=3*value/time_step/2
            elif(str2=="0"):
                matrixA_TV[place1-1][place1-1]+=3*value/time_step/2
            else:
                matrixA_TV[place1-1][place1-1]+=3*value/time_step/2
                matrixA_TV[place1-1][place2-1]-=3*value/time_step/2
                matrixA_TV[place2-1][place1-1]-=3*value/time_step/2
                matrixA_TV[place2-1][place2-1]+=3*value/time_step/2
            if(str1=="0"):
                vectorB_TV[place2-1]-=(2*C_pass[C_pointer][0]*value/time_step-C_pass[C_pointer][1]*value/time_step/2)
            elif(str2=="0"):
                vectorB_TV[place1-1]+=(2*C_pass[C_pointer][0]*value/time_step-C_pass[C_pointer][1]*value/time_step/2)
            else:
                vectorB_TV[place1-1]+=(2*C_pass[C_pointer][0]*value/time_step-C_pass[C_pointer][1]*value/time_step/2)
                vectorB_TV[place2-1]-=(2*C_pass[C_pointer][0]*value/time_step-C_pass[C_pointer][1]*value/time_step/2)
            C_pointer+=1
        #填入電流源的stamp
        #無須填入矩陣A的部分
        #僅填入vectorB
        elif(typelist[comp_number]=='I'):#暫時先不動，之後要改
            #用兩個str容器裝nodename
            str1=nodelist[i][0]
            str2=nodelist[i][1]
            #用place去找dic中nodename對應值
            place1=dic_node[str1]
            place2=dic_node[str2]
            if(components[i].sin!=0):
                if(str1=="0"):
                    vectorB[place2-1]+=components[i].sin*np.sin(t*radiant)
                elif(str2=="0"):
                    vectorB[place1-1]-=components[i].sin*np.sin(t*radiant)  
                else:
                    vectorB[place2-1]+=components[i].sin*np.sin(t*radiant)
                    vectorB[place1-1]-=components[i].sin*np.sin(t*radiant)
            else:
                if(str1=="0"):
                    vectorB[place2-1]+=components[i].dc
                elif(str2=="0"):
                    vectorB[place1-1]-=components[i].dc   
                else:
                    vectorB[place2-1]+=components[i].dc
                    vectorB[place1-1]-=components[i].dc
        #填入電壓源的stamp
        #須填入矩陣A的部分
        #同時填入vectorB
        elif(typelist[comp_number]=='V'):
            ##print("KKR")
            ##print(tv_namelist[i]," ",components[i].dc)
            element2=element2+1#作為一個指標用來定出元件所放置的位置
            ##print(element2)
            #用兩個str容器裝nodename
            str1=tv_nodelist[i][0]
            str2=tv_nodelist[i][1]
            #用place去找dic中nodename對應值
            place1=dic_node[str1]
            place2=dic_node[str2]
            comp_number=dic_component[tv_namelist[i]]
            #填入vectorB
            if(components[comp_number].sin[1]!=0):
                voff=components[comp_number].sin[0]
                amp=components[comp_number].sin[1]
                rads=components[comp_number].sin[2]*2*pi
                delay=components[comp_number].sin[3]
                damp=components[comp_number].sin[4]
                phase=components[comp_number].sin[5]
                vectorB_TV[element2-1]+=amp*np.sin(t*rads+phase)+voff
                #填入matrixA_TV
                if(str1=="0"):
                    matrixA_TV[element2-1][place2-1]-=1
                    matrixA_TV[place2-1][element2-1]-=1
                elif(str2=="0"):
                    matrixA_TV[element2-1][place1-1]+=1
                    matrixA_TV[place1-1][element2-1]+=1  
                else:
                    matrixA_TV[element2-1][place1-1]+=1
                    matrixA_TV[element2-1][place2-1]-=1
                    matrixA_TV[place1-1][element2-1]+=1
                    matrixA_TV[place2-1][element2-1]-=1
            elif(components[comp_number].pulse[6]!=0):
                #print("EUer ",components[comp_number].sin," ",components[comp_number])
                V1=components[comp_number].pulse[0]
                V2=components[comp_number].pulse[1]
                TD=components[comp_number].pulse[2]
                TR=components[comp_number].pulse[3]
                TF=components[comp_number].pulse[4]
                PW=components[comp_number].pulse[5]
                PER=components[comp_number].pulse[6]
                NP=components[comp_number].pulse[7]
                t_place=t-TD
                T_in_PER=t_place%PER
                if(t<=TD or T_in_PER==0):
                    vectorB_TV[element2-1]+=V1
                elif(t>TD and T_in_PER>0 and T_in_PER<(TR)):
                    vectorB_TV[element2-1]+=V1+(V2-V1)*(T_in_PER)/TR
                elif(t>TD and T_in_PER>=TR and T_in_PER<(TR+PW)):
                    vectorB_TV[element2-1]+=V2
                elif(t>TD and T_in_PER>=TR+PW and T_in_PER<(TR+PW+TF)):
                    vectorB_TV[element2-1]+=V2+(V1-V2)/TF*(T_in_PER-TR-PW)
                elif(t>TD and T_in_PER>=TR+PW+TF and T_in_PER<PER):
                    vectorB_TV[element2-1]+=V1
                
                #填入matrixA_TV
                if(str1=="0"):
                    matrixA_TV[element2-1][place2-1]-=1
                    matrixA_TV[place2-1][element2-1]-=1
                elif(str2=="0"):
                    matrixA_TV[element2-1][place1-1]+=1
                    matrixA_TV[place1-1][element2-1]+=1  
                else:
                    matrixA_TV[element2-1][place1-1]+=1
                    matrixA_TV[element2-1][place2-1]-=1
                    matrixA_TV[place1-1][element2-1]+=1
                    matrixA_TV[place2-1][element2-1]-=1    
            elif(len(components[comp_number].pwl)>0):
                time_table=[]
                V_table=[]
                for u in range(len(components[comp_number].pwl)):
                        V_table.append(components[comp_number].pwl[u][1])
                        
                        time_table.append(components[comp_number].pwl[u][0])
                #填入vector TV
                #print(V_table)
                #print(time_table)
                voltage_ptr=0
                for u in range(len(time_table)):
                    if(t<time_table[u]):
                        voltage_ptr=u-1#指向第n個電壓
                        break
                if(t>=time_table[-1]):
                    vectorB_TV[element2-1]+=V_table[-1]
                elif(V_table[voltage_ptr]==V_table[voltage_ptr+1]):
                    #print(V_table[voltage_ptr+1])
                    vectorB_TV[element2-1]+=V_table[voltage_ptr]
                
                else:
                    voltage_gap=(V_table[voltage_ptr+1]-V_table[voltage_ptr])
                    time_gap=time_table[voltage_ptr+1]-time_table[voltage_ptr]
                    voltage_in_variation=V_table[voltage_ptr]+(t-time_table[voltage_ptr])*voltage_gap/time_gap
                    print(voltage_in_variation)
                    vectorB_TV[element2-1]+=voltage_in_variation
                #填入matrixA_TV
                if(str1=="0"):
                    matrixA_TV[element2-1][place2-1]-=1
                    matrixA_TV[place2-1][element2-1]-=1
                elif(str2=="0"):
                    matrixA_TV[element2-1][place1-1]+=1
                    matrixA_TV[place1-1][element2-1]+=1  
                else:
                    matrixA_TV[element2-1][place1-1]+=1
                    matrixA_TV[element2-1][place2-1]-=1
                    matrixA_TV[place1-1][element2-1]+=1
                    matrixA_TV[place2-1][element2-1]-=1
            else:#需要修改，給DC電壓源在外面填
                ##print("xxx")
                vectorB[element2-1]+=components[comp_number].dc
            #填入matrixA
            
                if(str1=="0"):
                    matrixA[element2-1][place2-1]-=1
                    matrixA[place2-1][element2-1]-=1
                elif(str2=="0"):
                    matrixA[element2-1][place1-1]+=1
                    matrixA[place1-1][element2-1]+=1  
                else:
                    matrixA[element2-1][place1-1]+=1
                    matrixA[element2-1][place2-1]-=1
                    matrixA[place1-1][element2-1]+=1
                    matrixA[place2-1][element2-1]-=1
        #填入電感的stamp
        #須填入矩陣A的部分
        #無須填入vectorB
        elif(typelist[comp_number]=='L'):
            ##print("KKR")
            element2=element2+1#作為一個指標用來定出元件所放置的位置
            #用兩個str容器裝nodename
            str1=tv_nodelist[i][0]
            str2=tv_nodelist[i][1]
            #用place去找dic中nodename對應值
            place1=dic_node[str1]
            place2=dic_node[str2]
            comp_number=dic_component[tv_namelist[i]]
            value=components[comp_number].value
            check_current.append(element2-1)
            #填入matrixA
            if(str1=="0"):
                matrixA_TV[element2-1][place2-1]-=1
                matrixA_TV[place2-1][element2-1]-=1
                matrixA_TV[element2-1][element2-1]-=3*value/time_step/2
                vectorB_TV[element2-1]-=(2*L_pass[L_pointer][0]*value/time_step+L_pass[L_pointer][1]*value/time_step/2)
            elif(str2=="0"):
                matrixA_TV[element2-1][place1-1]+=1
                matrixA_TV[place1-1][element2-1]+=1 
                matrixA_TV[element2-1][element2-1]-=3*value/time_step/2
                vectorB_TV[element2-1]-=(2*L_pass[L_pointer][0]*value/time_step+L_pass[L_pointer][1]*value/time_step/2)

            else:
                matrixA_TV[element2-1][place1-1]+=1
                matrixA_TV[element2-1][place2-1]-=1
                matrixA_TV[place1-1][element2-1]+=1
                matrixA_TV[place2-1][element2-1]-=1
                matrixA_TV[element2-1][element2-1]-=3*value/time_step/2
                vectorB_TV[element2-1]-=(2*L_pass[L_pointer][0]*value/time_step+L_pass[L_pointer][1]*value/time_step/2)
            L_pointer+=1
        else:
            continue
    #填入
    #LU解法==========================================   
    ##print(matrixA_TV)
    #fill in mos
    if(len(nonlin_namelist)==0 and len(mos_device_list)==0):
        matrixA_mask=np.zeros((matrix_size,matrix_size),dtype=float)
        for s in range(element1):
            matrixA_mask[s][s]+=10**-12
        
        matrixA_mix=matrixA_TV+matrixA_Non+matrixA+matrixA_mask
        ##print(matrixA_mix)
        vectorB_mix=vectorB+vectorB_Non+vectorB_TV
        A=np.array(matrixA_mix,dtype=float)
        P, L, U = linalg.lu(A)
        B=np.array(vectorB_mix)
        Y=linalg.solve_triangular(L,P.T@B,lower=True)
        X=linalg.solve_triangular(U,Y,lower=False)
    else:
        cond=True
        iti=0
        cap=0
        while(cond):
            cond_false=0
            modified=0
            for j in range(len(nonlin_namelist)):
                ##print(t," ",j)
                ##print(j)
                ##print(x_0[diode_number])
                #用兩個str容器裝nodename
                str1=nonlin_nodelist[j][0]
                str2=nonlin_nodelist[j][1]
                #用place去找dic中nodename對應值
                place1=dic_node[str1]
                place2=dic_node[str2]
                ##print((nonlin_his[j][-1]/0.0258528413-32.2361913))
                #a=np.array(nonlin_his[j][-1]/0.0258528413-32.2361913, dtype=float)
                a=nonlin_his[j][-1]/0.0258528413-32.2361913
                tempa=np.exp(a)/0.0258528413
                tempb=np.exp(a)-10**(-14)-nonlin_his[j][-1]*np.exp(a)/0.0258528413
                #填入matrixA
                if(str1=="0"):
                    matrixA_Non[place2-1][place2-1]+=tempa
                elif(str2=="0"):
                    matrixA_Non[place1-1][place1-1]+=tempa
                else:
                    matrixA_Non[place1-1][place1-1]+=tempa
                    matrixA_Non[place1-1][place2-1]-=tempa
                    matrixA_Non[place2-1][place1-1]-=tempa
                    matrixA_Non[place2-1][place2-1]+=tempa
                #填入vectorB
                if(str1=="0"):
                    vectorB_Non[place2-1]+=tempb
                elif(str2=="0"):
                    vectorB_Non[place1-1]-=tempb
                else:
                    vectorB_Non[place2-1]+=tempb
                    vectorB_Non[place1-1]-=tempb
            #fill in mos
            for i in range(len(mos_device_list)):
                comp_number=dic_component[mos_device_list[i]]
                place1 = dic_node[mos_device_node[i][0]]
                place2 = dic_node[mos_device_node[i][1]]
                place3 = dic_node[mos_device_node[i][2]]
                place4 = dic_node[mos_device_node[i][3]]
                D = mos_device_node[i][0]
                G = mos_device_node[i][1]
                S = mos_device_node[i][2]
                B = mos_device_node[i][3]
                a11 = 0
                a12 = 0
                b1 = 0

                if mos_device_list[i][1] == "p":
        # =============================================================================
        #                 Vgs = vgs_his[u]
        #                 Vds = vds_his[u]
        # =============================================================================
                    Vgs = vgs_his[i][-1]
                    Vds = vds_his[i][-1]
                    ##print(vgs_his," ",vds_his)
                    #mos_ptr += 1
                    width = unit_symbol(components[comp_number].args["w"])
                    length = unit_symbol(components[comp_number].args["l"])
                    ##print(width," ",length," p",)
                    W_L = width / length
                    
                    # Cutoff
                    if Vgs < Vt:
                        ##print("cutoff_p",Vgs," ",Vds)
                        a12 = 0
                        a11 = sp.exp((Vgs - Vt) / 0.0258528413 - 32.2361913) / 0.0258528413
                        b1 = -(a11 * 0.0258528413 - a11 * (Vgs))
                    
                    # Triode
                    elif (Vgs - Vt) > Vds:
                        ##print("triode_p",Vgs," ",Vds)
                        a11 = W_L * k_p * Vds
                        a12 = W_L * k_p * ((Vgs - Vt) - Vds)
                        b1 = -(
                            W_L * k_p * ((Vgs - Vt) * Vds - (Vds) ** 2 / 2)
                            - a11 * Vgs
                            - a12 * Vds
                        )

                    # Sat
                    elif (Vgs - Vt) <= Vds:
                        ##print("sat_p",Vgs," ",Vds)
                        a11 = (
                            (1 / 2)
                            * W_L
                            * k_p
                            * (
                                2 * (Vgs - Vt) * (1 + 1 / va_p * (Vds - (Vgs - Vt)))
                                - 1 / va_p * (Vgs - Vt) ** 2
                            )
                        )
                        a12 = (1 / 2) * W_L * k_p * 1 / va_p * (Vgs - Vt) ** 2
                        b1 =-(
                            (1 / 2)
                            * W_L
                            * k_p
                            * (Vgs - Vt) ** 2
                            * (1 + 1 / va_p * (Vds - (Vgs - Vt)))
                            - a11 * Vgs
                            - a12 * Vds
                        )
                        
                elif mos_device_list[i][1] == "n":
                    ##print(vgs_his)
                    Vgs = vgs_his[i][-1]
                    Vds = vds_his[i][-1]
                    ##print(vgs_his," ",vds_his)
                    #mos_ptr += 1
                    width = unit_symbol(components[comp_number].args["w"])
                    length = unit_symbol(components[comp_number].args["l"])
                    W_L = width / length
                    ##print(width," ",length," n",)
                    
                    # Cutoff
                    if Vgs < Vt:
                        ##print("cutoff")
                        a12 = 0
                        a11 = sp.exp((Vgs - Vt) / 0.0258528413 - 32.2361913) / 0.0258528413
                        b1 = a11 * 0.0258528413 - a11 * (Vgs)
                    
                    # Triode
                    elif (Vgs - Vt) > Vds:
                        ##print("triode")
                        a11 = W_L * k_n * Vds
                        a12 = W_L * k_n * ((Vgs - Vt) - Vds)
                        b1 = (
                            W_L * k_n * ((Vgs - Vt) * Vds - (Vds) ** 2 / 2)
                            - a11 * Vgs
                            - a12 * Vds
                        )
                    
                    # Sat
                    elif (Vgs - Vt) <= Vds:
                        ##print("sat")
                        a11 = (
                            (1 / 2)
                            * W_L
                            * k_n
                            * (
                                2 * (Vgs - Vt) * (1 + 1 / va_n * (Vds - (Vgs - Vt)))
                                - 1 / va_n * (Vgs - Vt) ** 2
                            )
                        )
                        a12 = (1 / 2) * W_L * k_n * 1 / va_n * (Vgs - Vt) ** 2
                        b1 = (
                            (1 / 2)
                            * W_L
                            * k_n
                            * (Vgs - Vt) ** 2
                            * (1 + 1 / va_n * (Vds - (Vgs - Vt)))
                            - a11 * Vgs
                            - a12 * Vds
                        )
                ##print(a11," ",a12)
                if D == "0" and G == "0" and S == "0":
                    
                    pass
                elif D == "0":
                    if G == "0" and S != "0":
                        matrixA_mos[place3 - 1][place3 - 1] += a11 + a12
                        vectorB_mos[place3 - 1] += b1
                    elif S == "0" and G != "0":
                        pass
                    elif S != "0" and G != "0":
                        matrixA_mos[place3 - 1][place1 - 1] -= a11
                        matrixA_mos[place3 - 1][place3 - 1] += a11 + a12
                        vectorB_mos[place3 - 1] += b1
                elif G == "0":
                    if D == "0" and S != "0":
                        matrixA_mos[place3 - 1][place3 - 1] += a11 + a12
                        vectorB_mos[place3 - 1] += b1
                    elif S == "0" and D != "0":
                        matrixA_mos[place1 - 1][place1 - 1] += a12
                        vectorB_mos[place1 - 1] -= b1
                    elif S != "0" and D != "0":
                        matrixA_mos[place1 - 1][place1 - 1] += a12
                        matrixA_mos[place1 - 1][place3 - 1] -= a11 + a12
                        vectorB_mos[place1 - 1] -= b1
                        matrixA_mos[place3 - 1][place1 - 1] -= a12
                        matrixA_mos[place3 - 1][place3 - 1] += a11 + a12
                        vectorB_mos[place3 - 1] += b1
                elif S == "0":
                    if D == "0" and G != "0":
                        pass
                    elif G == "0" and D != "0":
                        matrixA_mos[place1 - 1][place1 - 1] += a12  # Change
                        vectorB_mos[place1 - 1] -= b1
                    elif G != "0" and D != "0":
                        matrixA_mos[place1 - 1][place2 - 1] += a11  # Change
                        matrixA_mos[place1 - 1][place1 - 1] += a12  # Change  
                        vectorB_mos[place1 - 1] -= b1
                else:
                    matrixA_mos[place1 - 1][place2 - 1] += a11
                    matrixA_mos[place1 - 1][place1 - 1] += a12
                    matrixA_mos[place1 - 1][place3 - 1] -= (a11 + a12)
                    vectorB_mos[place1 - 1] -= b1
                    matrixA_mos[place3 - 1][place2 - 1] -= a11
                    matrixA_mos[place3 - 1][place1 - 1] -= a12
                    matrixA_mos[place3 - 1][place3 - 1] += (a11 + a12)
                    vectorB_mos[place3 - 1] += b1  
                    #fill in mos end========================        
            matrixA_mask=np.zeros((matrix_size,matrix_size),dtype=float)
            for s in range(element1):
                matrixA_mask[s][s]+=10**-12
            matrixA_mix=matrixA_TV+matrixA_Non+matrixA+matrixA_mask+matrixA_mos
            vectorB_mix=vectorB+vectorB_Non+vectorB_TV+vectorB_mos
            A=np.array(matrixA_mix,dtype=float)
            ##print("matrixA",A)
            
            condition_number=np.linalg.cond(A)
            ##print("cond: ",condition_number)
# =============================================================================
#             if((condition_number>=10**16 and iti==0) or time_step_new>increment ):
#                 cond_false=1
#                 if(time_step_new>10**-5):
#                     time_step_new/=10
#                 else:
#                     time_step_new=time_step_new
#                 #print("wrong",time_step_new)
#             elif(condition_number<10**8 and time_step_new<increment and iti==0):
#                 time_step_new*=10
#                 if(time_step_new>=increment):
#                     time_step_new=increment
#                 #print("speed_up",time_step_new)
# =============================================================================
            #elif(condition_number<10**8):
                #time_step*=10
            ##print(iti)
            P, L, U = linalg.lu(A)
            B=np.array(vectorB_mix)
            Y=linalg.solve_triangular(L,P.T@B,lower=True)
            X_list=linalg.solve_triangular(U,Y,lower=False)
            X=X_list
            ##print(X_list)
            #更新vgs vds=============================
            for u in range(len(mos_device_list)):
                comp_number=dic_component[mos_device_list[u]]
                place1 = dic_node[mos_device_node[u][0]]
                place2 = dic_node[mos_device_node[u][1]]
                place3 = dic_node[mos_device_node[u][2]]
                place4 = dic_node[mos_device_node[u][3]]
                D = mos_device_node[u][0]
                G = mos_device_node[u][1]
                S = mos_device_node[u][2]
                B = mos_device_node[u][3]
                if D == "0" and G == "0" and S == "0":
                        pass
                    
                elif D == "0":
                    if G == "0" and S != "0":
                        if mos_device_list[u][1] == "p":
                            vds_his[u][0]=vds_his[u][1]
                            vds_his[u][1] = (X[place3 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (X[place3 - 1])
                        else:
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = (-X[place3 - 1])
                            vgs_his[u][0]=vgs_his[u][1]
                            vgs_his[u][1] = (-X[place3 - 1])
                    elif S == "0" and G != "0":
                        if mos_device_list[u][1] == "p":
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1]=0
                            vgs_his[u][0] =vgs_his[u][1]
                            vgs_his[u][1] = (-X[place2 - 1])
                        else:
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = (0)
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (X[place2 - 1])
                    elif S != "0" and G != "0":
                        if mos_device_list[u][1] == "p":
                            vds_his[u][0] =vds_his[u][1]
                            vds_his[u][1] = (X[place3 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (
                                X[place3 - 1] - X[place2 - 1]
                            )
                        else:
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = (-X[place3 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (
                                X[place2 - 1] - X[place3 - 1]
                            )
                            
                elif G == "0":
                    if D == "0" and S != "0":
                        if mos_device_list[u][1] == "p":
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = (X[place3 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (X[place3 - 1])
                        else:
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = (-X[place3 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (-X[place3 - 1])
                    elif S == "0" and D != "0":
                        if mos_device_list[u][1] == "p":
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = (-X[place1 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (0)
                        else:
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = (X[place1 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (0)
                    elif S != "0" and D != "0":
                        if mos_device_list[u][1] == "p":
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = (
                                X[place3 - 1] - X[place1 - 1]
                            )
                            vgs_his[u][0] = vgs_his[u][1] 
                            vgs_his[u][1] = (X[place3 - 1])
                        else:
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = (
                                X[place1 - 1] - X[place3 - 1]
                            )
                            vgs_his[u][0] = vgs_his[u][1] 
                            vgs_his[u][1] = (-X[place3 - 1])
                            
                elif S == "0":
                    if D == "0" and G != "0":
                        if mos_device_list[u][1] == "p":
                            vds_his[u][0] = vds_his[u][1]
                            vds_his[u][1] = 0
                            vgs_his[u][0] = vgs_his[u][1] 
                            vgs_his[u][1] = (-X[place2 - 1])
                        else:
                            vds_his[u][0] =vds_his[u][1] 
                            vds_his[u][1] = (0)
                            vgs_his[u][0] =vgs_his[u][1] 
                            vgs_his[u][1] = (X[place2 - 1] - 0)
                    elif G == "0" and D != "0":
                        if mos_device_list[u][1] == "p":
                            vds_his[u][0] =vds_his[u][1] 
                            vds_his[u][1] = (-X[place1 - 1] + 0)
                            vgs_his[u][0]=vgs_his[u][1]
                            vgs_his[u][1] = (0)
                        else:
                            vds_his[u][0] =vds_his[u][1] 
                            vds_his[u][1] = (X[place1 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (0)
                    elif G != "0" and D != "0":
                        if mos_device_list[u][1] == "p":
                            vds_his[u][0] = vds_his[u][1] 
                            vds_his[u][1] = (-X[place1 - 1])
                            vgs_his[u][0] = vgs_his[u][1] 
                            vgs_his[u][1] = (-X[place2 - 1])
                        else:
                            vds_his[u][0] = vds_his[u][1] 
                            vds_his[u][1] = (X[place1 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (X[place2 - 1])
                            
                else:
                    if mos_device_list[u][1] == "p":
                        vds_his[u][0] = vds_his[u][1] 
                        vds_his[u][1] = (
                            X[place3 - 1] - X[place1 - 1]
                        )
                        vgs_his[u][0] = vgs_his[u][1]
                        vgs_his[u][1] = (
                            X[place3 - 1] - X[place2 - 1]
                        )
                        ##print("check1: ",X[place3 - 1] - X[place1 - 1]," ",X[place3 - 1] - X[place2 - 1])
                    else:
                        # #print(X[place1-1]-X[place3-1])
                        vds_his[u][0] = vds_his[u][1]
                        vds_his[u][1] = (
                            X[place1 - 1] - X[place3 - 1]
                        )
                        vgs_his[u][0] = vgs_his[u][1]
                        vgs_his[u][1] = (
                            X[place2 - 1] - X[place3 - 1]
                        )
            #更新vgs vds END=====================
            #確認跨壓是否過大並修正
            for u in range(len(nonlin_nodelist)):
                ##print(x_0[diode_number])
                #用兩個str容器裝nodename
                str1=nonlin_nodelist[u][0]
                str2=nonlin_nodelist[u][1]
                #用place去找dic中nodename對應值
                place1=dic_node[str1]
                place2=dic_node[str2]
                
                #先填入解出值
                if(str1=="0"):
                    temp1=nonlin_his[u][-1]
                    temp2=-X_list[place2-1]
                    nonlin_his[u]=[temp1,temp2]
                elif(str2=="0"):
                    temp1=nonlin_his[u][-1]
                    temp2=X_list[place1-1]
                    nonlin_his[u]=[temp1,temp2]
                else:
                    temp1=nonlin_his[u][-1]
                    temp2=X_list[place1-1]-X_list[place2-1]
                    nonlin_his[u]=[temp1,temp2]
                ##print(nonlin_his)
                #檢查過大並修正
                if(nonlin_his[u][-1]>0.734):#更正掉過大的x_0
                    ##print(nonlin_his[u][-1])
                    #if(iti==0):
                        #time_step_new/=1000
                        ##print("wrong",time_step_new)
                        ##print(nonlin_his[u][-1])
                    #timefault.append([t,con])
                    x_klast=nonlin_his[u][-2]
                    ##print("fwve:",x_klast)
                    a=x_klast/0.0258528413-32.2361913
                    b=x_klast/0.0258528413
                    #a=np.array( a, dtype=float)
                    #b=np.array( b, dtype=float)
                    #a=x_klast/0.0258528413-32.2361913
                    #b=x_klast/0.0258528413
                    i_temp=nonlin_his[u][-1]*np.exp(a)/0.0258528413+np.exp(a)-10**(-14)-np.exp(a)*x_klast/0.0258528413
                    ##print(i_temp+10**(-14))
                    if(i_temp<0):
                        i_temp=10**-14
                    #tempk=1+i_temp/10**(-14)
                    queew=np.log(i_temp)-np.log(10**-14)
                    tempk=1+np.exp(queew)
                    ##print(tempk)
                    #tempk=np.array( tempk, dtype=float)
                    if(tempk<0):
                        ##print(condition_number)
                        tempk=-tempk
                        nonlin_his[u][-1]=0.0258528413*np.log(tempk)
                        ##print("a",0.0258528413*np.log(tempk))
                    elif(tempk==0.0):
                        ##print(condition_number)
                        #timefault.append(t)
                        tempk=2
                        nonlin_his[u][-1]=0.732
                        ##print("b",0.0258528413*np.log(tempk))
                    else:
                        nonlin_his[u][-1]=0.0258528413*np.log(tempk)
                       
                        ##print("c",0.0258528413*np.log(tempk))
                    modified=1
            
            ##print(nonlin_his)
            non_diff=[]
            for u in range(len(nonlin_nodelist)):
                non_diff.append(nonlin_his[u][-1]-nonlin_his[u][-2])
            if(len(non_diff)==0):
                non_diff=[0]
            ##print(non_diff)
            vgs_diff=[]
            vds_diff=[]
            for u in range (len(vgs_his)):
                vgs_diff.append(abs(vgs_his[u][1]-vgs_his[u][0]))
                vds_diff.append(abs(vds_his[u][1]-vds_his[u][0]))
            if(len(vgs_diff)==0):
                vgs_diff=[0]
                vds_diff=[0]
            ptr_max_vgs=0
            ptr_max_vds=0
            for u in range(len(vgs_diff)):
                if(max(vgs_diff)==vgs_diff[u]):
                    ptr_max_vgs=u
                if(max(vds_diff)==vds_diff[u]):
                    ptr_max_vds=u
            v_max_vgs=0
            v_max_vds=0
            if(len(vgs_his)!=0):
                v_max_vgs=max(abs(vgs_his[ptr_max_vgs][0]),abs(vgs_his[ptr_max_vgs][1]))
                v_max_vds=max(abs(vds_his[ptr_max_vds][0]),abs(vds_his[ptr_max_vds][1]))
            #print(max(vgs_diff)," x ",max(vds_diff)," ",iti)
            if((iti==0 or modified==1 or max(non_diff)>10**(-6) or 
                abs(min(non_diff))>10**(-6)
                or max(vgs_diff)> v_max_vgs*10**-3 or max(vds_diff)>v_max_vds*10**-3)
               and iti<=100):
                #
                
                iti+=1
                #生成Non線性元件矩陣A
                matrixA_Non=np.zeros((matrix_size,matrix_size),dtype=float)
                #生成Non線性元件矩陣B
                vectorB_Non=np.zeros(matrix_size,dtype=float)
                #生成Mos元件矩陣B
                vectorB_mos = np.zeros(matrix_size)
                #生成Mos元件矩陣A
                matrixA_mos = np.zeros((matrix_size, matrix_size))
                modified=0
                pass
            else:
                ##print("pass ",t/end_time*100,"%")
                for u in range(len(nonlin_nodelist)):
                    nonlin_his[u]=[nonlin_his[u][-1]]
                for u in range(len(mos_device_list)):
                    vgs_his[u]=vgs_his[u]
                    vds_his[u]=vds_his[u]
# =============================================================================
#                 if(iti<2):
#                     revive_step+=1
#                 if(iti>=2 and cond_false==0):
#                     revive_step=0
#                     if(time_step_new>10**-5):
#                         time_step_new/=1000
#                     else:
#                         time_step_new=time_step_new
#                     #print("wrong",time_step_new)
#                 elif(revive_step==1 and time_step_new<increment and cond_false==0):
#                     revive_step=0
#                     time_step_new*=1000
#                     if(time_step_new>=increment):
#                         time_step_new=increment
# =============================================================================
                ##print("iti over",iti)
                
                iti=0
                cond=False
                modified=0
                return (X_list)