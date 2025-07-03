from utils import unit_symbol, plot_picture, plot_picture_trans, plot_picture_dc
import sympy as sp
from mpmath import mp
import numpy as np
import matplotlib.pyplot as plt
import re
from scipy import linalg

def Dc_bias(
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
    x_0,
    start_V,
    C_pass,
    L_pass
):
    checknode=[]
    check_current=[]
    # PARAMETER
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
    print(k_n)
    result = []
    matrixA = np.zeros((matrix_size, matrix_size))
    matrixA_nonlin = np.zeros((matrix_size, matrix_size))
    
    # 生成向量B
    # Ax=B
    vectorB = np.zeros(matrix_size)
    vectorB_nonlin = np.zeros(matrix_size)
    
    x_axis = []
    x_0 = []
    x_0his = []
    
    iti_list = []
    x_0_list = []
    diode_number = 0
    # mos_ptr=0
    times = 0
    partial=-1
    firstTime = True
    first_print_maxtrix = True
    while  partial <100:
        partial+=1
        print("biasing: ",partial/100*100)
        vectorB_mos = np.zeros(matrix_size)
        matrixA_mos = np.zeros((matrix_size, matrix_size))
        iti_time = 0
        node_check = []
        # mos_check=[]
        mos_node = []
        mos_device = []
        mos_ptr = 0
        temp12 = element2
        
        
        # 填MNA
        for i in range(element):
            # 填入電阻的stamp
            # 填入矩陣A的部分
            if namelist[i][0] == "r":
                # 用兩個str容器裝nodename
                str1 = nodelist[i][0]
                str2 = nodelist[i][1]
                # 用place去找dic中nodename對應值
                place1 = dic_node[str1]
                place2 = dic_node[str2]
                value = components[i].value
                if str1 == "0":
                    matrixA[place2 - 1][place2 - 1] += 1 / value
                elif str2 == "0":
                    matrixA[place1 - 1][place1 - 1] += 1 / value
                else:
                    matrixA[place1 - 1][place1 - 1] += 1 / value
                    matrixA[place1 - 1][place2 - 1] -= 1 / value
                    matrixA[place2 - 1][place1 - 1] -= 1 / value
                    matrixA[place2 - 1][place2 - 1] += 1 / value
            # 填入電容的stamp
            # 填入矩陣A的部分
            elif namelist[i][0] == "c":
                str1 = nodelist[i][0]
                str2 = nodelist[i][1]
                if(partial==0):
                    checknode.append([str1,str2])
                pass

            # 填入電流源的stamp
            # 無須填入矩陣A的部分
            # 僅填入vectorB
            elif namelist[i][0] == "i":
                # 用兩個str容器裝nodename
                str1 = nodelist[i][0]
                str2 = nodelist[i][1]
                # 用place去找dic中nodename對應值
                place1 = dic_node[str1]
                place2 = dic_node[str2]
                value = components[i].dc
                if components[i].dc != 0:
                    if namelist[i] != object_source:
                        if str1 == "0":
                            vectorB[place2 - 1] += value*partial/100
                        elif str2 == "0":
                            vectorB[place1 - 1] -= value*partial/100
                        else:
                            vectorB[place2 - 1] += value*partial/100
                            vectorB[place1 - 1] -= value*partial/100
                    else:
                        if str1 == "0":
                            vectorB[place2 - 1] += 0
                        elif str2 == "0":
                            vectorB[place1 - 1] -= 0
                        else:
                            vectorB[place2 - 1] += 0
                            vectorB[place1 - 1] -= 0
                else:
                    if str1 == "0":
                        vectorB[place2 - 1] += 0
                    elif str2 == "0":
                        vectorB[place1 - 1] -= 0
                    else:
                        vectorB[place2 - 1] += 0
                        vectorB[place1 - 1] -= 0
            # 填入電壓源的stamp
            # 須填入矩陣A的部分
            # 同時填入vectorB
            elif namelist[i][0] == "v":
                element2 = element2 + 1  # 作為一個指標用來定出元件所放置的位置
                # rint(element2)
                # 用兩個str容器裝nodename
                str1 = nodelist[i][0]
                str2 = nodelist[i][1]
                # 用place去找dic中nodename對應值
                place1 = dic_node[str1]
                place2 = dic_node[str2]
                value = components[i].dc
                # 填入vectorB 
                if components[i].dc != 0:
                    if namelist[i] != object_source:
                        vectorB[element2 - 1] += value*(partial/100)
                    else:
                        vectorB[element2 - 1] += 0
                elif namelist[i] == object_source:
                    vectorB[element2 - 1] += start_V*(partial/100)
                elif (components[i].pulse[6] != 0 ):
                    #print("V2",components[i].pulse[1]*(partial/1000))
                    #vectorB[element2 - 1] += components[i].pulse[0]*(partial/1000)
                    if(components[i].pulse[2] == 0):    
                        vectorB[element2 - 1] += components[i].pulse[0]*(partial/100)
                        print("V2",components[i].pulse[1]*(partial/100))
                    else:
                        vectorB[element2 - 1] += components[i].pulse[0]*(partial/100)
                        print("V1",components[i].pulse[0]*(partial/100))
                elif(len(components[i].pwl)>0):
                    vectorB[element2 - 1]+=components[i].pwl[0][1]*partial/100
                else:
                    vectorB[element2 - 1] += 0
                # 填入matrixA
                if str1 == "0":
                    matrixA[element2 - 1][place2 - 1] -= 1
                    matrixA[place2 - 1][element2 - 1] -= 1
                elif str2 == "0":
                    matrixA[element2 - 1][place1 - 1] += 1
                    matrixA[place1 - 1][element2 - 1] += 1
                else:
                    matrixA[element2 - 1][place1 - 1] += 1
                    matrixA[element2 - 1][place2 - 1] -= 1
                    matrixA[place1 - 1][element2 - 1] += 1
                    matrixA[place2 - 1][element2 - 1] -= 1
            # 填入電感的stamp
            # 須填入矩陣A的部分
            # 無須填入vectorB
            elif namelist[i][0] == "l":
                element2 = element2 + 1  # 作為一個指標用來定出元件所放置的位置
                # 用兩個str容器裝nodename
                str1 = nodelist[i][0]
                str2 = nodelist[i][1]
                # 用place去找dic中nodename對應值
                place1 = dic_node[str1]
                place2 = dic_node[str2]
                if(partial==0):
                    check_current.append(element2-1)
                # 填入matrixA
                if str1 == "0":
                    matrixA[element2 - 1][place2 - 1] -= 1
                    matrixA[place2 - 1][element2 - 1] -= 1
                    matrixA[element2 - 1][element2 - 1] -= 0
                elif str2 == "0":
                    matrixA[element2 - 1][place1 - 1] += 1
                    matrixA[place1 - 1][element2 - 1] += 1
                    matrixA[element2 - 1][element2 - 1] -= 0
                else:
                    matrixA[element2 - 1][place1 - 1] += 1
                    matrixA[element2 - 1][place2 - 1] -= 1
                    matrixA[place1 - 1][element2 - 1] += 1
                    matrixA[place2 - 1][element2 - 1] -= 1
                    matrixA[element2 - 1][element2 - 1] -= 0
            elif namelist[i][0] == "d":
                # print(x_0[diode_number])
                # 用兩個str容器裝nodename
                str1 = nodelist[i][0]
                str2 = nodelist[i][1]
                # 用place去找dic中nodename對應值
                place1 = dic_node[str1]
                place2 = dic_node[str2]
                if partial == 0:
                    if (place1 - place2) <= 0:
                        x_0.append(0)
                    else:
                        x_0.append(0)
                else:
                    pass
                tempa = (
                    sp.exp(x_0[diode_number] / 0.0258528413 - 32.2361913) / 0.0258528413
                )
                tempb = (
                    sp.exp(x_0[diode_number] / 0.0258528413 - 32.2361913)
                    - 10 ** (-14)
                    - sp.exp(x_0[diode_number] / 0.0258528413 - 32.2361913)
                    / 0.0258528413
                    * x_0[diode_number]
                )
                # 填入matrixA
                if str1 == "0":
                    matrixA_nonlin[place2 - 1][place2 - 1] += tempa
                    node_check.append(["0", str2])
                elif str2 == "0":
                    matrixA_nonlin[place1 - 1][place1 - 1] += tempa
                    node_check.append([str1, "0"])
                else:
                    matrixA_nonlin[place1 - 1][place1 - 1] += tempa
                    matrixA_nonlin[place1 - 1][place2 - 1] -= tempa
                    matrixA_nonlin[place2 - 1][place1 - 1] -= tempa
                    matrixA_nonlin[place2 - 1][place2 - 1] += tempa
                    node_check.append([str1, str2])
                # 填入vectorB
                if str1 == "0":
                    vectorB_nonlin[place2 - 1] += tempb
                elif str2 == "0":
                    vectorB_nonlin[place1 - 1] -= tempb
                else:
                    vectorB_nonlin[place2 - 1] += tempb
                    vectorB_nonlin[place1 - 1] -= tempb
                diode_number += 1
            
            # 加mos
            elif namelist[i][0] == "m":
                # print()
                D = nodelist[i][0]
                G = nodelist[i][1]
                S = nodelist[i][2]
                B = nodelist[i][3]
                
                place1 = dic_node[D]
                place2 = dic_node[G]
                place3 = dic_node[S]
                place4 = dic_node[B]
                # 紀錄mos的各節點
                mos_node.append([D, G, S, B])
                mos_device.append(i)

            else:
                continue

        # =============================================================================
        #     A=np.array(matrixA,dtype=complex)
        #     P, L, U = linalg.lu(A)
        #     #print(L)
        #     #print(U)
        #     B=np.array(vectorB)
        #     Y=linalg.solve_triangular(L,P.T@B,lower=True)
        #     X=linalg.solve_triangular(U,Y,lower=False)
        # =============================================================================
        # print("maxtrixA: ", matrixA)
        # print("maxtrixA_nonlin: ", matrixA_nonlin)
        
        # print("mos_node: ", mos_node)
        for u in range(len(mos_node)):
            # 紀錄mos的各節點
            # mos_node.append([D,G,S,B])

            place1 = dic_node[mos_node[u][0]]
            place2 = dic_node[mos_node[u][1]]
            place3 = dic_node[mos_node[u][2]]
            place4 = dic_node[mos_node[u][3]]  #voltage==start_V
            #print(mos_node[u][0]," ",mos_node[u][1]," ",mos_node[u][2]," ",mos_node[u][3])
            # 計算初值firstTime == 
            if  partial == 0  :
                matrixA_mask=np.zeros((matrix_size, matrix_size))
                for k in range(element2):
                    matrixA_mask[k][k]+=10**-9
                vectorB_check = vectorB + vectorB_nonlin  # vectorB=[1,1,1,2,2]
                matrixA_check = matrixA + matrixA_nonlin + matrixA_mask
                #print(matrixA_check)
                # 先解一次矩陣 +matrixA_mask
                A = np.array(matrixA_check)
                P, L, U = linalg.lu(A)
                # print(L)
                # print(U)
                B = np.array(vectorB_check)
                Y = linalg.solve_triangular(L, P.T @ B, lower=True)
                D = mos_node[u][0]
                G = mos_node[u][1]
                S = mos_node[u][2]
                B = mos_node[u][3]
                print(u," G: ",G,"D: ",D,"S: ",S)
                vector_check = linalg.solve_triangular(U, Y, lower=False)
                #print("firstTime_A:", A)
                #print("firstTime_B:", B)
                #print(dic_node)
                #print("vector_check(first_X): ", vector_check)
                # mos_cnt+=1;
                if D == "0" and G == "0" and S == "0":
                    pass
                
                elif D == "0":
                    if G == "0" and S != "0":
                        if namelist[mos_device[u]][1] == "p":
                            vds_his.append([vector_check[place3 - 1],vector_check[place3 - 1]])
                            vgs_his.append([vector_check[place3 - 1],vector_check[place3 - 1]])
                        else:
                            vds_his.append([-vector_check[place3 - 1],-vector_check[place3 - 1]])
                            vgs_his.append([-vector_check[place3 - 1],-vector_check[place3 - 1]])
                    elif S == "0" and G != "0":
                        if namelist[mos_device[u]][1] == "p":
                            vds_his.append([0,0])
                            vgs_his.append([-vector_check[place2 - 1],-vector_check[place2 - 1]])
                        else:
                            vds_his.append([0,0])
                            vgs_his.append([vector_check[place2 - 1],vector_check[place2 - 1]])
                    elif S != "0" and G != "0":
                        if namelist[mos_device[u]][1] == "p":
                            vds_his.append([vector_check[place3 - 1],vector_check[place3 - 1]])
                            vgs_his.append(
                                [vector_check[place3 - 1] - vector_check[place2 - 1],
                                 vector_check[place3 - 1] - vector_check[place2 - 1]]
                            )
                        else:
                            vds_his.append([-vector_check[place3 - 1],-vector_check[place3 - 1]])
                            vgs_his.append(
                                [vector_check[place2 - 1] - vector_check[place3 - 1],
                                 vector_check[place2 - 1] - vector_check[place3 - 1]]
                            )
                            
                elif G == "0":
                    if D == "0" and S != "0":
                        if namelist[mos_device[u]][1] == "p":
                            vds_his.append([vector_check[place3 - 1],vector_check[place3 - 1]])
                            vgs_his.append([vector_check[place3 - 1],vector_check[place3 - 1]])
                        else:
                            vds_his.append([-vector_check[place3 - 1],-vector_check[place3 - 1]])
                            vgs_his.append([-vector_check[place3 - 1],-vector_check[place3 - 1]])
                    elif S == "0" and D != "0":
                        if namelist[mos_device[u]][1] == "p":
                            vds_his.append([-vector_check[place1 - 1],-vector_check[place1 - 1]])
                            vgs_his.append([0,0])
                        else:
                            vds_his.append([vector_check[place1 - 1],vector_check[place1 - 1]])
                            vgs_his.append([0,0])
                    elif S != "0" and D != "0":
                        if namelist[mos_device[u]][1] == "p":
                            vds_his.append(
                                [vector_check[place3 - 1] - vector_check[place1 - 1],
                                 vector_check[place3 - 1] - vector_check[place1 - 1]]
                            )
                            vgs_his.append([vector_check[place3 - 1],vector_check[place3 - 1]])
                        else:
                            vds_his.append(
                                [vector_check[place1 - 1] - vector_check[place3 - 1],
                                 vector_check[place1 - 1] - vector_check[place3 - 1]]
                            )
                            vgs_his.append([-vector_check[place3 - 1],-vector_check[place3 - 1]])
                
                elif S == "0":
                    if D == "0" and G != "0":
                        if namelist[mos_device[u]][1] == "p":
                            vds_his.append([0,0])
                            vgs_his.append([-vector_check[place2 - 1],-vector_check[place2 - 1]])
                        else:
                            vds_his.append([0,0])
                            vgs_his.append([vector_check[place2 - 1],vector_check[place2 - 1]])
                    elif G == "0" and D != "0":
                        if namelist[mos_device[u]][1] == "p":
                            vds_his.append([-vector_check[place1 - 1],-vector_check[place1 - 1]])
                            vgs_his.append([0,0])
                        else:
                            vds_his.append([vector_check[place1 - 1],vector_check[place1 - 1]])
                            vgs_his.append([0,0])
                    elif G != "0" and D != "0":
                        if namelist[mos_device[u]][1] == "p":
                            vds_his.append([-vector_check[place1 - 1],-vector_check[place1 - 1]])
                            vgs_his.append([-vector_check[place2 - 1],-vector_check[place2 - 1]])
                            #print(vgs_his[-1]," ",vds_his[-1])
                        else:
                            #print(vector_check[place1 - 1] - 0," ",vector_check[place2 - 1])
                            vds_his.append([vector_check[place1 - 1],vector_check[place1 - 1]])
                            vgs_his.append([vector_check[place2 - 1],vector_check[place2 - 1]])
                else:
                    if namelist[mos_device[u]][1] == "p":
                        vds_his.append(
                            [vector_check[place3 - 1] - vector_check[place1 - 1],
                             vector_check[place3 - 1] - vector_check[place1 - 1]]
                        )
                        vgs_his.append(
                            [vector_check[place3 - 1] - vector_check[place2 - 1],
                             vector_check[place3 - 1] - vector_check[place2 - 1]]
                        )
                    else:
                        # print(vector_check[place1-1]-vector_check[place3-1])
                        vds_his.append(
                            [vector_check[place1 - 1] - vector_check[place3 - 1],
                             vector_check[place1 - 1] - vector_check[place3 - 1]]
                        )
                        vgs_his.append(
                            [vector_check[place2 - 1] - vector_check[place3 - 1],
                             vector_check[place2 - 1] - vector_check[place3 - 1]]
                        )
                
                firstTime = False
                #continue
            
            else:
                #firstTime = True
                D = mos_node[u][0]
                G = mos_node[u][1]
                S = mos_node[u][2]
                B = mos_node[u][3]
                #print(u," G: ",G,"D: ",D,"S: ",S)
                pass
            
            
            a11 = 0
            a12 = 0
            b1 = 0

            #print(u," G: ",G,"D: ",D,"S: ",S)
            if namelist[mos_device[u]][1] == "p":
# =============================================================================
#                 Vgs = vgs_his[u]
#                 Vds = vds_his[u]
# =============================================================================
                #print(vgs_his)
                Vgs = vgs_his[u][-1]
                Vds = vds_his[u][-1]
                #print(vgs_his," ",vds_his)
                mos_ptr += 1
                width = unit_symbol(components[mos_device[u]].args["w"])
                length = unit_symbol(components[mos_device[u]].args["l"])
                #print(width," ",length," ",Vgs," ",Vds)
                W_L = width / length
                
                # Cutoff
                if Vgs < Vt:
                    #print("cutoff_p",Vgs," ",Vds)
                    a12 = 0
                    a11 = sp.exp((Vgs - Vt) / 0.0258528413 - 32.2361913) / 0.0258528413
                    b1 = -(a11 * 0.0258528413 - a11 * (Vgs))
                
                # Triode
                elif (Vgs - Vt) >= Vds:
                    #print("triode_p",Vgs," ",Vds)
                    a11 = W_L * k_p * Vds
                    a12 = W_L * k_p * ((Vgs - Vt) - Vds)
                    b1 = -(
                        W_L * k_p * ((Vgs - Vt) * Vds - (Vds) ** 2 / 2)
                        - a11 * Vgs
                        - a12 * Vds
                    )

                # Sat
                elif (Vgs - Vt) < Vds:
                    #print("sat_p",Vgs," ",Vds)
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
                    
            elif namelist[mos_device[u]][1] == "n":
                #print(vgs_his)
                Vgs = vgs_his[u][-1]
                Vds = vds_his[u][-1]
                #print("checking :",Vgs," ",Vds)
                mos_ptr += 1
                width = unit_symbol(components[mos_device[u]].args["w"])
                length = unit_symbol(components[mos_device[u]].args["l"])
                W_L = width / length
                
                # Cutoff
                if Vgs < Vt:
                    #print("cutoff")
                    a12 = 0
                    a11 = sp.exp((Vgs - Vt) / 0.0258528413 - 32.2361913) / 0.0258528413
                    b1 = a11 * 0.0258528413 - a11 * (Vgs)
                
                # Triode
                elif (Vgs - Vt) >= Vds:
                    #print("triode")
                    a11 = W_L * k_n * Vds
                    a12 = W_L * k_n * ((Vgs - Vt) - Vds)
                    b1 = (
                        W_L * k_n * ((Vgs - Vt) * Vds - (Vds) ** 2 / 2)
                        - a11 * Vgs
                        - a12 * Vds
                    )
                
                # Sat
                elif (Vgs - Vt) < Vds:
                    #print("sat")
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
            #print(a11," ",a12)
            if D == "0" and G == "0" and S == "0":
                
                pass
            elif D == "0":
                #print(place1," ",place2," ",place3," ")
                if G == "0" and S != "0":
                    matrixA_mos[place3 - 1][place3 - 1] += a11 + a12
                    vectorB_mos[place3 - 1] += b1
                elif S == "0" and G != "0":
                    matrixA_mos[place2 - 1][place2 - 1] += 0
                    vectorB_mos[place2 - 1] += 0
                elif S != "0" and G != "0":
                    matrixA_mos[place2 - 1][place2 - 1] += 0
                    matrixA_mos[place2 - 1][place3 - 1] += 0
                    vectorB_mos[place2 - 1] += 0
                    matrixA_mos[place3 - 1][place1 - 1] -= a12
                    matrixA_mos[place3 - 1][place3 - 1] += a11 + a12
                    vectorB_mos[place3 - 1] += b1
            elif G == "0":
                #print(place1," ",place2," ",place3," ")
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
                #print(place1," ",place2," ",place3," ")
                if D == "0" and G != "0":
                    matrixA_mos[place2 - 1][place2 - 1] += 0
                    vectorB_mos[place2 - 1] += 0
                elif G == "0" and D != "0":
                    matrixA_mos[place1 - 1][place1 - 1] += a12  # Change
                    vectorB_mos[place1 - 1] -= b1
                elif G != "0" and D != "0":
                    matrixA_mos[place2 - 1][place1 - 1] += 0
                    matrixA_mos[place2 - 1][place2 - 1] += 0
                    vectorB_mos[place2 - 1] += 0
                    matrixA_mos[place1 - 1][place2 - 1] += a11  # Change
                    matrixA_mos[place1 - 1][place1 - 1] += a12  # Change  
                    vectorB_mos[place1 - 1] -= b1
            else:
                #print(place1," ",place2," ",place3," ")
                matrixA_mos[place1 - 1][place2 - 1] += a11
                matrixA_mos[place1 - 1][place1 - 1] += a12
                matrixA_mos[place1 - 1][place3 - 1] -= (a11 + a12)
                vectorB_mos[place1 - 1] -= b1
                matrixA_mos[place3 - 1][place2 - 1] -= a11
                matrixA_mos[place3 - 1][place1 - 1] -= a12
                matrixA_mos[place3 - 1][place3 - 1] += (a11 + a12)
                vectorB_mos[place3 - 1] += b1
                
        
        X_list = []
        vectorB_mix = np.zeros(matrix_size)
        matrixA_mix = np.zeros((matrix_size, matrix_size))
        
        condition = True
        if len(node_check) != 0  or len(vgs_his)!=0:
            condition = True
        else:
            #print(abs(-voltage)/abs(start_V-end_V)*100)
            matrixA_mask=np.zeros((matrix_size, matrix_size))
            for k in range(element2):
                matrixA_mask[k][k]+=10**-9
            condition = False
            matrixA_mix = matrixA + matrixA_nonlin + matrixA_mos+matrixA_mask
            vectorB_mix = vectorB + vectorB_nonlin + vectorB_mos
            #print(matrixA_mix)
            # 歸零非線性元素矩陣
            matrixA_nonlin = np.zeros((matrix_size, matrix_size))
            vectorB_nonlin = np.zeros(matrix_size)
            # 解線性系統
            A = np.array(matrixA_mix)
            P, L, U = linalg.lu(A)
            B = np.array(vectorB_mix)
            
# =============================================================================
#             if first_print_maxtrix == True : 
#                 print("A: ", A)
#                 print("B: ", B)
#                 print(b1)
#                 first_print_maxtrix = False
# =============================================================================
            #print("A: ", A)
            #print("B: ", B)
            Y = linalg.solve_triangular(L, P.T @ B, lower=True)
            X = linalg.solve_triangular(U, Y, lower=False)
# =============================================================================
#             if(True):
#                 print("==============================================")
#                 print(f"Vd = {voltage},")
#                 # print(f"a11 = {a11},")
#                 # print(f"a12 = {a12},")
#                 print(f"A: {A}")
#                 print(f"B: {B}")
#                 print(f"X: {X}")
# =============================================================================
            
            

            # print("Vds_his: ", vds_his)
            # print("Vgs_his: ", vgs_his)
            
            # 檢驗收斂
            X_list.append(X)
            #print(X[2])
            
        
            
        while condition == True:
            #print(abs(start_V-voltage)/abs(start_V-end_V)*100)
            # 合併線性矩陣與非線性元件矩陣
            # =============================================================================
            #         for row in range(len(vectorB)):
            #             vectorB_mix[row]=vectorB[row]+vectorB_nonlin[row]
            #             for column in range(len(vectorB)):
            #                 matrixA_mix[row][column]=matrixA[row][column]+matrixA_nonlin[row][column]
            # =============================================================================
            matrixA_mask=np.zeros((matrix_size, matrix_size))
            for k in range(element2):
                matrixA_mask[k][k]+=10**-9
            matrixA_mix = matrixA + matrixA_nonlin+matrixA_mos+matrixA_mask
            vectorB_mix = vectorB + vectorB_nonlin+vectorB_mos
            #print("MA:",matrixA_mix)
            # print("Mb:",vectorB_mix)
            # 歸零非線性元素矩陣
            vectorB_mos = np.zeros(matrix_size)
            matrixA_mos = np.zeros((matrix_size, matrix_size))
            matrixA_nonlin = np.zeros((matrix_size, matrix_size))
            vectorB_nonlin = np.zeros(matrix_size)
            # print(matrixA_nonlin)
            # print(vectorB_nonlin)
            # 解線性系統
            A = np.array(matrixA_mix)
            #print("A: ",A)
            #print(max(A[0]))
            # print(A)
            P, L, U = linalg.lu(A)
            # print(L)
            # print(U)
            B = np.array(vectorB_mix)
            Y = linalg.solve_triangular(L, P.T @ B, lower=True)
            X = linalg.solve_triangular(U, Y, lower=False)
            # print("X:", X) #here
            # 檢驗收斂
            X_list.append(X)
            # x_0his為x(k-1)
            x_0his.append(x_0)
            # print("X_0:",x_0)
            # print("x_0his:",x_0his)
            # 根據node_check更新x_0
            x_0 = []  # 先清空x_0
            # x_0his為x(k-1)
            # 檢查x_0有沒有過大
            #print(dic_node)
            #print("x: ", X)
            for u in range(len(mos_node)):
                D = mos_node[u][0]
                G = mos_node[u][1]
                S = mos_node[u][2]
                B = mos_node[u][3]
                place1 = dic_node[D]
                place2 = dic_node[G]
                place3 = dic_node[S]
                place4 = dic_node[B] 
                #print("k"," G: ",G,"D: ",D,"S: ",S)
                if D == "0" and G == "0" and S == "0":
                        pass
                    
                elif D == "0":
                    if G == "0" and S != "0":
                        if namelist[mos_device[u]][1] == "p":
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
                        if namelist[mos_device[u]][1] == "p":
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
                        if namelist[mos_device[u]][1] == "p":
                            vds_his[u][0] =vds_his[u][1]
                            vds_his[u][1] = (X[place3 - 1])
                            vgs_his[u][0] = vgs_his[u][1]
                            vgs_his[u][1] = (
                                X[place3 - 1] - X[place3 - 1]
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
                        if namelist[mos_device[u]][1] == "p":
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
                        if namelist[mos_device[u]][1] == "p":
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
                        if namelist[mos_device[u]][1] == "p":
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
                        if namelist[mos_device[u]][1] == "p":
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
                        if namelist[mos_device[u]][1] == "p":
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
                        if namelist[mos_device[u]][1] == "p":
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
                    if namelist[mos_device[u]][1] == "p":
                        vds_his[u][0] = vds_his[u][1] 
                        vds_his[u][1] = (
                            X[place3 - 1] - X[place1 - 1]
                        )
                        vgs_his[u][0] = vgs_his[u][1]
                        vgs_his[u][1] = (
                            X[place3 - 1] - X[place2 - 1]
                        )
                        #print("check1: ",X[place3 - 1] - X[place1 - 1]," ",X[place3 - 1] - X[place2 - 1])
                    else:
                        # print(X[place1-1]-X[place3-1])
                        vds_his[u][0] = vds_his[u][1]
                        vds_his[u][1] = (
                            X[place1 - 1] - X[place3 - 1]
                        )
                        vgs_his[u][0] = vgs_his[u][1]
                        vgs_his[u][1] = (
                            X[place2 - 1] - X[place3 - 1]
                        )
                
                
            for u in range(len(node_check)):
                # 用兩個str容器裝nodename
                str1 = node_check[u][0]
                str2 = node_check[u][1]
                # 用place去找dic中nodename對應值
                place1 = dic_node[str1]
                place2 = dic_node[str2]
                if str1 == "0":
                    value = -X_list[-1][place2 - 1]
                    x_0.append(value)
                elif str2 == "0":
                    value = X_list[-1][place1 - 1]
                    x_0.append(value)
                else:
                    value = X_list[-1][place1 - 1] - X_list[-1][place2 - 1]
                    x_0.append(value)
                if x_0[-1] < 0.732:
                    pass
                else:  # 更正掉過大的x_0
                    x_klast = x_0his[-1][u]
                    i_temp = (
                        x_0[-1]
                        * sp.exp(x_klast / 0.0258528413 - 32.2361913)
                        / 0.0258528413
                        + 10 ** (-14) * sp.exp(x_klast / 0.0258528413)
                        - 10 ** (-14)
                        - 10 ** (-14)
                        * sp.exp(x_klast / 0.0258528413)
                        / 0.0258528413
                        * x_klast
                    )
                    temp = 1 + i_temp / 10 ** (-14)
                    x_0[-1] = 0.0258528413 * mp.log(temp)
            # print("X_0:",x_0)
            x_0fake = np.array(x_0)
            x_0fakehis = np.array(x_0his[-1])
            x_0diff = x_0fake - x_0fakehis
            if(len(x_0diff)==0):
                x_0diff=[0]
            vgs_diff=[]
            vds_diff=[]
            for u in range (len(vgs_his)):
                vgs_diff.append(abs(vgs_his[u][1]-vgs_his[u][0]))
                vds_diff.append(abs(vds_his[u][1]-vds_his[u][0]))
            if(len(vgs_diff)==0):
                vgs_diff=[0]
                vds_diff=[0]
            
            #print("maxi: ",max(vgs_diff)," ",max(vds_diff))
            #print(vgs_his," ",vds_his)
            ptr_max_vgs=0
            ptr_max_vds=0
            for u in range(len(vgs_diff)):
                if(max(vgs_diff)==vgs_diff[u]):
                    ptr_max_vgs=u
                if(max(vds_diff)==vds_diff[u]):
                    ptr_max_vds=u
            v_max_vgs=max(abs(vgs_his[ptr_max_vgs][1]),abs(vgs_his[ptr_max_vgs][1]))
            v_max_vds=max(abs(vds_his[ptr_max_vds][1]),abs(vds_his[ptr_max_vds][1]))
            #print(max(vgs_diff)," x ",max(vds_diff)," ",iti)
# =============================================================================
#             if((iti==0 or modified==1 or max(non_diff)>10**(-6) or 
#                 abs(min(non_diff))>10**(-6)
#                 or max(vgs_diff)> v_max_vgs*10**-3 or max(vds_diff)>v_max_vds*10**-3)
#                and iti<=100):
# =============================================================================
            if ((len(x_0his) == 1 or (
                max(x_0diff) >= 10 ** (-6) or abs(min(x_0diff)) >= 10 ** (-6)
                or max(vgs_diff)>=v_max_vgs*10**-6 or max(vds_diff)>=v_max_vds*10**-6) 
                    and iti_time<100
            )) :
                
                # print("length of x0",len(x_0))
                # print("length of x0his",len(x_0his))and iti_time<100
                iti_time += 1
                #print("iti_time: ",iti_time," ",voltage," ",max(vgs_diff)," ",max(vds_diff))
# =============================================================================
#                 [ 1.80176898e+00  0.00000000e+00  1.80175347e+00 -3.81658350e-30
#                   1.39262362e+00  1.80180000e+00 -6.79794607e-09  0.00000000e+00
#                   3.81658350e-21]
# =============================================================================
                if(iti_time==99):
                    print("biasing_false: ",partial/100*100)
                    print()
                for u in range(len(mos_node)):
                    # 紀錄mos的各節點
                    # mos_node.append([D,G,S,B])
                    D = mos_node[u][0]
                    G = mos_node[u][1]
                    S = mos_node[u][2]
                    B = mos_node[u][3]
                    place1 = dic_node[D]
                    place2 = dic_node[G]
                    place3 = dic_node[S]
                    place4 = dic_node[B]
                    
                    a11 = 0
                    a12 = 0
                    b1 = 0

                    if namelist[mos_device[u]][1] == "p":
        # =============================================================================
        #                 Vgs = vgs_his[u]
        #                 Vds = vds_his[u]
        # =============================================================================
                        Vgs = vgs_his[u][-1]
                        Vds = vds_his[u][-1]
                        #print(Vgs," ",Vds)
                        mos_ptr += 1
                        width = unit_symbol(components[mos_device[u]].args["w"])
                        length = unit_symbol(components[mos_device[u]].args["l"])
                        #print(width," ",length," ",Vgs," ",Vds)
                        W_L = width / length
                        #print(vgs_his[u]," ",vgs_his[u])
                        # Cutoff
                        if Vgs < Vt:
                            
                        
                            a12 = 0
                            a11 = sp.exp((Vgs - Vt) / 0.0258528413 - 32.2361913) / 0.0258528413
                            b1 = -(a11 * 0.0258528413 - a11 * (Vgs))
                        
                        # Triode
                        elif (Vgs - Vt) >= Vds:
                            
                            a11 = W_L * k_p * Vds
                            a12 = W_L * k_p * ((Vgs - Vt) - Vds)
                            b1 = -(
                                W_L * k_p * ((Vgs - Vt) * Vds - (Vds) ** 2 / 2)
                                - a11 * Vgs
                                - a12 * Vds
                            )

                        # Sat
                        elif (Vgs - Vt) < Vds:
                            
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
                            
                    elif namelist[mos_device[u]][1] == "n":
                        Vgs = vgs_his[u][-1]
                        Vds = vds_his[u][-1]
                        mos_ptr += 1
                        width = unit_symbol(components[mos_device[u]].args["w"])
                        length = unit_symbol(components[mos_device[u]].args["l"])
                        W_L = width / length
                        
                        # Cutoff
                        if Vgs < Vt:
                            #print("cutoff")
                            a12 = 0
                            a11 = sp.exp((Vgs - Vt) / 0.0258528413 - 32.2361913) / 0.0258528413
                            b1 = a11 * 0.0258528413 - a11 * (Vgs)
                        
                        # Triode
                        elif (Vgs - Vt) >= Vds:
                            #print("triode")
                            a11 = W_L * k_n * Vds
                            a12 = W_L * k_n * ((Vgs - Vt) - Vds)
                            b1 = (
                                W_L * k_n * ((Vgs - Vt) * Vds - (Vds) ** 2 / 2)
                                - a11 * Vgs
                                - a12 * Vds
                            )
                        
                        # Sat
                        elif (Vgs - Vt) < Vds:
                            #print("sat")
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
                    #print(a11," ",a12)
                    if D == "0" and G == "0" and S == "0":
                        
                        pass
                    elif D == "0":
                        if G == "0" and S != "0":
                            matrixA_mos[place3 - 1][place3 - 1] += a11 + a12
                            vectorB_mos[place3 - 1] += b1
                        elif S == "0" and G != "0":
                            matrixA_mos[place2 - 1][place2 - 1] += 0
                            vectorB_mos[place2 - 1] += 0
                        elif S != "0" and G != "0":
                            matrixA_mos[place2 - 1][place2 - 1] += 0
                            matrixA_mos[place2 - 1][place3 - 1] += 0
                            vectorB_mos[place2 - 1] += 0
                            matrixA_mos[place3 - 1][place1 - 1] -= a12
                            matrixA_mos[place3 - 1][place3 - 1] += a11 + a12
                            vectorB_mos[place3 - 1] += b1
                    elif G == "0":
                        if D == "0" and S != "0":
                            matrixA_mos[place3 - 1][place3 - 1] += a11 + a12
                            vectorB_mos[place3 - 1] += b1
                        elif S == "0" and D != "0":
                            matrixA_mos[place1 - 1][place2 - 1] += a11
                            vectorB_mos[place1 - 1] -= b1
                        elif S != "0" and D != "0":
                            matrixA_mos[place1 - 1][place2 - 1] += a11
                            matrixA_mos[place1 - 1][place3 - 1] -= a11 + a12
                            vectorB_mos[place1 - 1] -= b1
                            matrixA_mos[place3 - 1][place2 - 1] -= a11
                            matrixA_mos[place3 - 1][place3 - 1] += a11 + a12
                            vectorB_mos[place3 - 1] += b1
                    elif S == "0":
                        if D == "0" and G != "0":
                            matrixA_mos[place2 - 1][place2 - 1] += 0
                            vectorB_mos[place2 - 1] += 0
                        elif G == "0" and D != "0":
                            matrixA_mos[place1 - 1][place1 - 1] += a12 # Change
                            vectorB_mos[place1 - 1] -= b1
                        elif G != "0" and D != "0":
                            matrixA_mos[place2 - 1][place1 - 1] += 0
                            matrixA_mos[place2 - 1][place2 - 1] += 0
                            vectorB_mos[place2 - 1] += 0
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
                for node_pair in range(len(node_check)):
                    # 用兩個str容器裝nodename
                    str1 = node_check[node_pair][0]
                    str2 = node_check[node_pair][1]
                    # print("nodes:",str1," ",str2)
                    # print(str2)
                    # 用place去找dic中nodename對應值
                    place1 = dic_node[str1]
                    place2 = dic_node[str2]
                    # print("node1:" ,place1)
                    # print("node2:" ,place2)
                    tempa = (
                        sp.exp(x_0[node_pair] / 0.0258528413 - 32.2361913)
                        / 0.0258528413
                    )
                    tempb = (
                        sp.exp(x_0[node_pair] / 0.0258528413 - 32.2361913)
                        - 10 ** (-14)
                        - sp.exp(x_0[node_pair] / 0.0258528413 - 32.2361913)
                        / 0.0258528413
                        * x_0[node_pair]
                    )
                    # 填入matrixA
                    if str1 == "0":
                        matrixA_nonlin[place2 - 1][place2 - 1] += tempa
                    elif str2 == "0":
                        matrixA_nonlin[place1 - 1][place1 - 1] += tempa
                    else:
                        matrixA_nonlin[place1 - 1][place1 - 1] += tempa
                        matrixA_nonlin[place1 - 1][place2 - 1] -= tempa
                        matrixA_nonlin[place2 - 1][place1 - 1] -= tempa
                        matrixA_nonlin[place2 - 1][place2 - 1] += tempa
                    # 填入vectorB
                    if str1 == "0":
                        vectorB_nonlin[place2 - 1] += tempb
                    elif str2 == "0":
                        vectorB_nonlin[place1 - 1] -= tempb
                    else:
                        vectorB_nonlin[place2 - 1] += tempb
                        vectorB_nonlin[place1 - 1] -= tempb
            else:
                # 將解出的跨壓作為下一個voltage的初始值
                #print(x_0his)
                x_0 = []
                x_0his = []
                print(iti_time,end=" ")
                # result.append(X_list[-1])
                #print("end_iti 1: ",vgs_his)
                for u in range(len(mos_node)):
                    D = mos_node[u][0]
                    G = mos_node[u][1]
                    S = mos_node[u][2]
                    B = mos_node[u][3]
                    place1 = dic_node[D]
                    place2 = dic_node[G]
                    place3 = dic_node[S]
                    place4 = dic_node[B]
                    if D == "0" and G == "0" and S == "0":
                            pass
                        
                    elif D == "0":
                        if G == "0" and S != "0":
                            if namelist[mos_device[u]][1] == "p":
                                vds_his[u] = [X_list[-1][place3 - 1],X_list[-1][place3 - 1]]
                                vgs_his[u] = [X_list[-1][place3 - 1],X_list[-1][place3 - 1]]
                            else:
                                vds_his[u] = [-X_list[-1][place3 - 1],-X_list[-1][place3 - 1]]
                                vgs_his[u] = [-X_list[-1][place3 - 1],-X_list[-1][place3 - 1]]
                        elif S == "0" and G != "0":
                            if namelist[mos_device[u]][1] == "p":
                                vds_his[u]=[0,0]
                                vgs_his[u] = [-X_list[-1][place2 - 1],-X_list[-1][place2 - 1]]
                            else:
                                vds_his[u] = [0,0]
                                vgs_his[u] = [X_list[-1][place2 - 1],X_list[-1][place2 - 1]]
                        elif S != "0" and G != "0":
                            if namelist[mos_device[u]][1] == "p":
                                vds_his[u] = [X_list[-1][place3 - 1],X_list[-1][place3 - 1]]
                                vgs_his[u] = [
                                    X_list[-1][place3 - 1] - X_list[-1][place3 - 1],
                                    X_list[-1][place3 - 1] - X_list[-1][place3 - 1]
                                ]
                            else:
                                vds_his[u] = [-X_list[-1][place3 - 1],-X_list[-1][place3 - 1]]
                                vgs_his[u] = [
                                    X_list[-1][place2 - 1] - X_list[-1][place3 - 1],
                                    X_list[-1][place2 - 1] - X_list[-1][place3 - 1]
                                ]
                                
                    elif G == "0":
                        if D == "0" and S != "0":
                            if namelist[mos_device[u]][1] == "p":
                                vds_his[u]= [X_list[-1][place3 - 1],X_list[-1][place3 - 1]]
                                vgs_his[u] = [X_list[-1][place3 - 1],X_list[-1][place3 - 1]]
                            else:
                                vds_his[u]= [-X_list[-1][place3 - 1],-X_list[-1][place3 - 1]]
                                vgs_his[u] = [-X_list[-1][place3 - 1],-X_list[-1][place3 - 1]]
                        elif S == "0" and D != "0":
                            if namelist[mos_device[u]][1] == "p":
                                vds_his[u] = [-X_list[-1][place1 - 1],-X_list[-1][place1 - 1]]
                                vgs_his[u] = [0,0]
                            else:
                                vds_his[u] = [ X_list[-1][place1 - 1] , X_list[-1][place1 - 1] ]
                                vgs_his[u] = [0,0]
                        elif S != "0" and D != "0":
                            if namelist[mos_device[u]][1] == "p":
                                vds_his[u] = (
                                    [ X_list[-1][place3 - 1] - X_list[-1][place1 - 1],
                                     X_list[-1][place3 - 1] - X_list[-1][place1 - 1]]
                                     
                                )
                                vgs_his[u]= [X_list[-1][place3 - 1],X_list[-1][place3 - 1]]
                            else:
                                vds_his[u]= (
                                    [X_list[-1][place1 - 1] - X_list[-1][place3 - 1],
                                     X_list[-1][place1 - 1] - X_list[-1][place3 - 1] ]
                                )
                                vgs_his[u]= [-X_list[-1][place3 - 1],-X_list[-1][place3 - 1]]
                                
                    elif S == "0":
                        if D == "0" and G != "0":
                            if namelist[mos_device[u]][1] == "p":
                                vds_his[u] = [0,0]
                                vgs_his[u] = [-X_list[-1][place2 - 1],-X_list[-1][place2 - 1]]
                            else:
                                vds_his[u] = [0,0]
                                vgs_his[u] = [X_list[-1][place2 - 1],X_list[-1][place2 - 1]]
                        elif G == "0" and D != "0":
                            if namelist[mos_device[u]][1] == "p":
                                vds_his[u] = [-X_list[-1][place1 - 1],-X_list[-1][place1 - 1]]
                                vgs_his[u] = [0,0]
                            else:
                                vds_his[u] = [X_list[-1][place1 - 1],X_list[-1][place1 - 1]]
                                vgs_his[u] = [0,0]
                        elif G != "0" and D != "0":
                            if namelist[mos_device[u]][1] == "p":
                                vds_his[u] = [-X_list[-1][place1 - 1],-X_list[-1][place1 - 1]]
                                vgs_his[u] = [-X_list[-1][place2 - 1],-X_list[-1][place2 - 1]]
                            else:
                                
                                vds_his[u] = [X_list[-1][place1 - 1],X_list[-1][place1 - 1]]
                                vgs_his[u] = [X_list[-1][place2 - 1],X_list[-1][place2 - 1]]
                                #print()
                                #print("end_iti 3: ",vgs_his)
                                
                    else:
                        if namelist[mos_device[u]][1] == "p":
                            vds_his[u] = [
                                X_list[-1][place3 - 1] - X_list[-1][place1 - 1]
                                ,X_list[-1][place3 - 1] - X_list[-1][place1 - 1]
                                ]
                            vgs_his[u] = [
                                X_list[-1][place3 - 1] - X_list[-1][place2 - 1],
                                X_list[-1][place3 - 1] - X_list[-1][place2 - 1]
                                ]
                        else:
                            vds_his[u] = [
                                X_list[-1][place1 - 1] - X_list[-1][place3 - 1],
                                X_list[-1][place1 - 1] - X_list[-1][place3 - 1]
                                ]
                            vgs_his[u] = [
                                X_list[-1][place2 - 1] - X_list[-1][place3 - 1]
                                ,X_list[-1][place2 - 1] - X_list[-1][place3 - 1]
                                ]
                            
                #print("end_iti 2: ",vgs_his)
                for u in range(len(node_check)):
                    # 用兩個str容器裝nodename
                    str1 = node_check[u][0]
                    str2 = node_check[u][1]
                    # 用place去找dic中nodename對應值
                    place1 = dic_node[str1]
                    place2 = dic_node[str2]
                    if str1 == "0":
                        value = -X_list[-1][place2 - 1]
                        x_0.append(value)
                    elif str2 == "0":
                        value = X_list[-1][place1 - 1]
                        x_0.append(value)
                    else:
                        value = X_list[-1][place1 - 1] - X_list[-1][place2 - 1]
                        x_0.append(value)
                condition = False

        # 將矩陣及向量清空
        # print(X_list[-1])
        result.append(X_list[-1])
        # vgs_his = []
        # vds_his = []
        matrixA = np.zeros((matrix_size, matrix_size))
        vectorB = np.zeros(matrix_size)
        # print("voltage:",voltage)
        # print(result)
        # print(-result[-1][-1])
        iti_list.append(iti_time)
        

        
        X_list = []
        element2 = temp12
        diode_number = 0
        
        # 畫圖用的x矩陣加入當下頻率值
        x_0_list.append(x_0)
        times += 1
    counter=1
    dic_result={}
    for i in range(len(templist)):
        if(templist[i]=="0"):
            dic_result["gnd"]=0
            #counter=counter+1
            #continue
        else:
            dic_result["v("+templist[i]+")"]=counter
            counter=counter+1
        continue
    for i in range(len(namelist)):
        if(namelist[i][0]=="v"):
            dic_result["i("+namelist[i]+")"]=counter
            counter=counter+1
        elif(namelist[i][0]=="l"):
            dic_result["i("+namelist[i]+")"]=counter
            counter=counter+1
        else:
            continue
    #存C

    print(dic_result)   
    print("dc bias",result[-1])
    for i in range(len(checknode)):
        str1=checknode[i][0]
        str2=checknode[i][1]
        #用place去找dic中nodename對應值
        place1=dic_node[str1]
        place2=dic_node[str2]
        #print(i," : ",str1," ",str2," ",place1," ",place2)
        #check_current=[element2-1]
        #填入matrixA
        res=result[-1]
        #print(res)
        if(str1=="0"):
            print(i," : ",str1," ",str2," ",place1," ",place2)
            print("ccap voltage: ",-result[-1][place2-1])
            C_pass[i]=[-result[-1][place2-1],-result[-1][place2-1]]
        elif(str2=="0"):
            print(i," : ",str1," ",str2," ",place1," ",place2)
            print("ccap voltage: ",result[-1][place1-1])
            C_pass[i]=[result[-1][place1-1],result[-1][place1-1]]
        else:
            C_pass[i]=[result[-1][place1-1]-res[place2-1],result[-1][place1-1]-res[place2-1]]
    #存L
    #print(check_current)
    
    for i in range(len(check_current)):
        res=result[len(result)-1]
        L_pass[i]=[0,res[check_current[i]]]

    
    return x_0,vds_his,vgs_his,C_pass,L_pass,result[-1]


# def tran_analysis(templist,element2,element,dic_node,matrix_size,cmd,namelist,nodelist,circuit_name,valuelist,typelist,components,V_step,start_V,end_V,object_source):