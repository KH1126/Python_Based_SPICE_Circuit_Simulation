from schema import Capacitor, ComponentBase, Current, Diode, Inductor, Resistor, Voltage, Mos  # fmt: skip
from linker import link_program
# from utils import unit_symbol


def parse_args(args: list[str]) -> dict:
    kargs = {}
    for arg in args:
        arg = arg.split("=")
        if len(arg) == 1:
            #print(arg)
            raise ValueError("args must be key=value")
        kargs[arg[0]] = arg[1]
    return kargs


def parse_VandI(line: list[str]) -> Voltage | Current:
    args = {}
    args["name"] = line[0]
    args["nodes"] = [line[1], line[2]]
    p = 3
    if line[3] not in {"dc", "ac", "pulse", "sin"}:
        args["dc"] = line[3]
        p += 1
    else:
        while p < len(line) and line[p] in {"dc", "ac", "pulse", "sin"}:
            if line[p] == "sin":
               i=1
               while p+i<len(line) and line[p+i] not in {"dc", "ac", "pulse", "sin"}:
                   #print(i)
                   i+=1
               if(i-1<3):
                   raise ValueError("sin must be 3 args")
                   break
               else:
                   args[line[p]] = line[p + 1 : p+i]
               p+=i
            elif line[p] == "pulse":
               i=1
               while p+i<len(line) and line[p+i] not in {"dc", "ac", "pulse", "sin"}:
                   i+=1
               if(i-1<8):
                   raise ValueError("pulse must be 8 args")
                   break
               else:
                   args[line[p]] = line[p + 1 : p+i]
               p+=i

                
            else:
                args[line[p]] = line[p + 1]
                p += 2
    args["args"] = parse_args(line[p:])

    if line[0][0] == "v":
        return Voltage(**args)
    elif line[0][0] == "i":
        return Current(**args)
    elif line[0][0] == "r":
        return Resistor(**args)


def parse_diode(line: list[str]):
    args = {}
    args["name"] = line[0]
    args["nodes"] = [line[1], line[2]]
    args["mname"] = line[3]
    args["args"] = parse_args(line[4:])

    return Diode(**args)


def parse_mos(line: list[str]) -> ComponentBase:
    args = {}
    args["name"] = line[0]
    args["nodes"] = [line[i] for i in range(1, 5)]
    args["mname"] = line[5]
    args["args"] = parse_args(line[6:])

    return Mos(**args)


def parse_other(line: list[str]) -> ComponentBase:
    args = {}
    args["name"] = line[0]
    args["nodes"] = [line[1], line[2]]
    args["value"] = line[3]
    args["args"] = parse_args(line[4:])

    if line[0][0] == "l":
        return Inductor(**args)
    elif line[0][0] == "c":
        return Capacitor(**args)
    elif line[0][0] == "r":
        return Resistor(**args)


def parse(file_name):
    # 相對位置加檔名
    with open(file_name, "r") as f:
        #raw_lines = f.readlines()
        raw_lines=link_program(f.read()).split("\n")
        #print(raw_lines)
    # 删除第一行
    circuit_name = raw_lines.pop(0).strip()

    lines, model_array = [], []
    # 到".end"为止
    for line in raw_lines[: raw_lines.index(".end")]:
        # 删除注释
        if line[0] == "*":
            continue

        line = line.lower().split()
        if line[0] == ".model":
            model_array.append(line[1])
            model_array.append(line[2])
        else:
            lines.append(line)

    # 創立namelist確定元件及元件數rc
    namelist = list(zip(*lines))[0]
    print(lines)
    print(list(zip(*lines)))
    components = []
    for line in lines:
        if line[0][0] in {"v", "i"}:
            components.append(parse_VandI(line))
        elif line[0][0] == "d":
            components.append(parse_diode(line))
        elif line[0][0] == "m":
            components.append(parse_mos(line))
        else:
            components.append(parse_other(line))

    return circuit_name, components, namelist


if __name__ == "__main__":
    # file_name = "diode.txt"
    file_name = input("file name: ")
    circuit_name, components, namelist = parse(file_name)
    nodelist = [c.nodes for c in components]
    valuelist = [c.value for c in components]
    print("circuit_name: ",circuit_name)
    print("namelist: ",namelist)
    print("nodelist: ",nodelist)
    print("valuelist: ",valuelist)
    for component in components:
        print(component)
