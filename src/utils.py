import matplotlib.pyplot as plt
import numpy as np
def unit_symbol(symbol):
    if isinstance(symbol, (float, int)):
        return float(symbol)
    unit = {
        "n": 1e-9,
        "p": 1e-12,
        "f": 1e-15,
        "u": 1e-6,
        "m": 1e-3,
        "k": 1e3,
        "g": 1e9,
        "t": 1e12,
    }
    if symbol[-3:] == "meg":
        return float(symbol[:-3]) * 1e6
    elif symbol[-1] in unit:
        return float(symbol[:-1]) * unit[symbol[-1]]
    return float(symbol)
def plot_picture(templist,namelist,times,result,cmd,x,circuit_name):
    dic_result={}
    nodes=[]
    counter=1
    for i in range(len(templist)):
        if(templist[i]=="0"):
            dic_result["gnd"]=0
            nodes.append("gnd")
            #counter=counter+1
            #continue
        else:
            dic_result["v("+templist[i]+")"]=counter
            nodes.append("v("+templist[i]+")")
            counter=counter+1
        continue
    for i in range(len(namelist)):
        if(namelist[i][0]=="v"):
            dic_result["i("+namelist[i]+")"]=counter
            nodes.append("i("+namelist[i]+")")
            counter=counter+1
        elif(namelist[i][0]=="l"):
            dic_result["i("+namelist[i]+")"]=counter
            nodes.append("i("+namelist[i]+")")
            counter=counter+1
        else:
            continue
    while True:
        print(dic_result)
        #收集輸入字串
        object_value1 = ""
        object_value2 = ""
        while(True):
            input_value = input("what do you want \n")
            
            unit = ["+", "-", "*", "/"," "]
            if any(op in input_value for op in unit):
                    # 从运算符列表中找到实际使用的运算符
                operator = next(op for op in unit if op in input_value)
                test=1
            else:
                operator=""
                test=0
            # 首先检查是否以负号开头
            if input_value.startswith(operator) and test==1:
                object_value1 = "gnd"
                object_value2 = input_value[1:]
            else:
                # 检查输入值是否包含运算符
                if operator!="":
                    # 使用运算符进行分割
                    parts = input_value.split(operator)
                    object_value1 = parts[0]
                    object_value2 = parts[1] 
                else:
                    object_value1 = input_value
                    object_value2 = "gnd"
            if((object_value1 not in nodes) or (object_value2 not in nodes)):
                print("no such node exist! enter again!")
            else:
                break
        y=[]
        z=[]

        pointer1=dic_result[object_value1]
        pointer2=dic_result[object_value2]
        for i in range(times):
            if dic_result[object_value1]==0:
                y.append(0)
            else :
                y.append(result[i][pointer1-1])
            if dic_result[object_value2]==0:
                z.append(0)
            else :
                z.append(result[i][pointer2-1])
        for i in range(times):
            if(operator=="+"):
                y[i]=y[i]+z[i]
                y[i]=abs(y[i])
            elif(operator=="-"):
                y[i]=y[i]-z[i]
                y[i]=abs(y[i])
            elif(operator=="*"):
                y[i]=y[i]*z[i] 
                y[i]=abs(y[i])
            elif(operator=="/"):
                #y[i]=abs(y[i])/abs(z[i])
                #y[i]=abs(y[i]/z[i])
                y[i]=20*np.log10(abs(y[i])/abs(z[i]))
                
            else :
                y[i]=abs(y[i])
        if(cmd[0]=="ac" and (cmd[1]=="dec" or cmd[1]=="oct")):
            plt.xscale("log")
            plt.xlabel("Frequency (Hz in log scale)",fontsize=14)
        else:
            plt.xlabel("Frequency (Hz) ",fontsize=14)
        plt.title(circuit_name,fontsize=20)
# =============================================================================
# =============================================================================
        plt.ylabel("db(vo/vi)",fontsize=14)
        plt.plot(x,y,linewidth=1)
        # 實線
        plt.plot(x, y, linestyle='-', label=object_value1)
        
        # 虛線
        
        
        
        if(operator==" "):
            plt.plot(x, z, linestyle='--', label=object_value2)
            # 圖例
        file_name="RLC"+'_output.txt'
        file=open(file_name,'w')
        for i in range(len((x))):
            file.write(str(x[i])+" "+str(y[i])+"\n")
        file.close()
        plt.legend()     
        plt.grid()
        plt.rcParams['figure.dpi'] = 300
        plt.show()
        Continue=input("Continue or not, 1 for continue, 0 for end. ")
        if (Continue=="1"):
            continue
        else:
            break
def plot_picture_trans(templist,namelist,result,cmd,x,circuit_name):
    dic_result={}
    nodes=[]
    counter=1
    for i in range(len(templist)):
        if(templist[i]=="0"):
            dic_result["gnd"]=0
            nodes.append("gnd")
            #counter=counter+1
            #continue
        else:
            dic_result["v("+templist[i]+")"]=counter
            nodes.append("v("+templist[i]+")")
            counter=counter+1
        continue
    for i in range(len(namelist)):
        if(namelist[i][0]=="v"):
            dic_result["i("+namelist[i]+")"]=counter
            nodes.append("i("+namelist[i]+")")
            counter=counter+1
        elif(namelist[i][0]=="l"):
            dic_result["i("+namelist[i]+")"]=counter
            nodes.append("i("+namelist[i]+")")
            counter=counter+1
        else:
            continue
    while True:
        print(dic_result)
        #收集輸入字串
        object_value1 = ""
        object_value2 = ""
        while(True):
            input_value = input("what do you want \n")
            
            unit = ["+", "-", "*", "/"," "]
            if any(op in input_value for op in unit):
                    # 从运算符列表中找到实际使用的运算符
                operator = next(op for op in unit if op in input_value)
                test=1
            else:
                operator=""
                test=0
            # 首先检查是否以负号开头
            if input_value.startswith(operator) and test==1:
                object_value1 = "gnd"
                object_value2 = input_value[1:]
            else:
                # 检查输入值是否包含运算符
                if operator!="":
                    # 使用运算符进行分割
                    parts = input_value.split(operator)
                    object_value1 = parts[0]
                    object_value2 = parts[1] 
                else:
                    object_value1 = input_value
                    object_value2 = "gnd"
            if((object_value1 not in nodes) or (object_value2 not in nodes)):
                print("no such node exist! enter again!")
            else:
                break
        y=[]
        z=[]

        pointer1=dic_result[object_value1]
        pointer2=dic_result[object_value2]
        for i in range(len(x)):
            if dic_result[object_value1]==0:
                y.append(0)
            else :
                y.append(result[i][pointer1-1])
            if dic_result[object_value2]==0:
                z.append(0)
            else :
                z.append(result[i][pointer2-1])
        #print(y[-1]," ",z[-1])
        for i in range(len(x)):
            if(operator=="+"):
                y[i]=y[i]+z[i]
            elif(operator=="-"):
                y[i]=y[i]-z[i] 
            elif(operator=="*"):
                y[i]=y[i]*z[i] 
            elif(operator=="/"):
                y[i]=y[i]/z[i] 
            else :
                y[i]=y[i]
        plt.title(circuit_name,fontsize=20)
        plt.xlabel("time (sec)",fontsize=14)
# =============================================================================
# =============================================================================
        plt.ylabel(object_value1[0].upper(),fontsize=14)
        plt.plot(x,y,linewidth=1)
        # 實線
        plt.plot(x, y, linestyle='-', label=object_value1)
        
        # 虛線
        
        
        
        if(operator==" "):
            plt.plot(x, z, linestyle='--', label=object_value2)
            file_name=object_value1+'_output.txt'
            file=open(file_name,'w')
            for i in range(len((x))):
                file.write(str(x[i])+" "+str(z[i]))
        else:
            file_name=object_value1+'_output.txt'
            file=open(file_name,'w')
            for i in range(len((x))):
                file.write(str(x[i])+" "+str(y[i])+"\n")
            file.close()
            # 圖例
        
        plt.legend()     
        plt.grid()
        plt.rcParams['figure.dpi'] = 300
        plt.show()
        Continue=input("Continue or not, 1 for continue, 0 for end. ")
        while(Continue!="1" and Continue!="0"):
            print("false cmd")
            Continue=input("Continue or not, 1 for continue, 0 for end. ")
        if (Continue=="1"):
            continue
        elif(Continue=="0"):
            break
def plot_picture_dc(templist,namelist,times,result,cmd,x,circuit_name,object_source):
    dic_result={}
    nodes=[]
    counter=1
    for i in range(len(templist)):
        if(templist[i]=="0"):
            dic_result["gnd"]=0
            nodes.append("gnd")
            #counter=counter+1
            #continue
        else:
            dic_result["v("+templist[i]+")"]=counter
            nodes.append("v("+templist[i]+")")
            counter=counter+1
        continue
    for i in range(len(namelist)):
        if(namelist[i][0]=="v"):
            dic_result["i("+namelist[i]+")"]=counter
            nodes.append("i("+namelist[i]+")")
            counter=counter+1
        elif(namelist[i][0]=="l"):
            dic_result["i("+namelist[i]+")"]=counter
            nodes.append("i("+namelist[i]+")")
            counter=counter+1
        else:
            continue
    while True:
        print(dic_result)
        #收集輸入字串
        object_value1 = ""
        object_value2 = ""
        while(True):
            input_value = input("what do you want \n")
            
            unit = ["+", "-", "*", "/"," "]
            if any(op in input_value for op in unit):
                    # 从运算符列表中找到实际使用的运算符
                operator = next(op for op in unit if op in input_value)
                test=1
            else:
                operator=""
                test=0
            # 首先检查是否以负号开头
            if input_value.startswith(operator) and test==1:
                object_value1 = "gnd"
                object_value2 = input_value[1:]
            else:
                # 检查输入值是否包含运算符
                if operator!="":
                    # 使用运算符进行分割
                    parts = input_value.split(operator)
                    object_value1 = parts[0]
                    object_value2 = parts[1] 
                else:
                    object_value1 = input_value
                    object_value2 = "gnd"
            if((object_value1 not in nodes) or (object_value2 not in nodes)):
                print("no such node exist! enter again!")
            else:
                break
        y=[]
        z=[]

        pointer1=dic_result[object_value1]
        pointer2=dic_result[object_value2]
        for i in range(times):
            if dic_result[object_value1]==0:
                y.append(0)
            else :
                y.append(result[i][pointer1-1])
            if dic_result[object_value2]==0:
                z.append(0)
            else :
                z.append(result[i][pointer2-1])
        for i in range(times):
            if(operator=="+"):
                y[i]=y[i]+z[i]
            elif(operator=="-"):
                y[i]=y[i]-z[i] 
            elif(operator=="*"):
                y[i]=y[i]*z[i] 
            elif(operator=="/"):
                y[i]=y[i]/z[i] 
            else :
                y[i]=y[i]
        plt.title(circuit_name,fontsize=20)
        plt.xlabel("sweep source: "+object_source+"(V)",fontsize=14)
# =============================================================================
# =============================================================================
        plt.ylabel(object_value1[0].upper(),fontsize=14)
        plt.plot(x,y,linewidth=1)
        # 實線
        plt.plot(x, y, linestyle='-', label=object_value1)
        
        # 虛線
        
        
        
        if(operator==" "):
            plt.plot(x, z, linestyle='--', label=object_value2)
            # 圖例
        
        plt.legend()     
        plt.grid()
        plt.rcParams['figure.dpi'] = 300
        plt.show()
        Continue=input("Continue or not, 1 for continue, 0 for end. ")
        while(Continue!="1" and Continue!="0"):
            print("false cmd")
            Continue=input("Continue or not, 1 for continue, 0 for end. ")
        if (Continue=="1"):
            continue
        elif(Continue=="0"):
            break
        
            