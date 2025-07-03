from schema import Capacitor, ComponentBase, Current, Diode, Inductor, Resistor, Voltage, Mos  # fmt: skip
from linker import link_program
import os
# from utils import unit_symbol


def parse_args(args: list[str]) -> dict:
    kargs = {}
    for arg in args:
        arg = arg.split("=")
        if len(arg) == 1:
            raise ValueError("args must be key=value")
        kargs[arg[0]] = arg[1]
    return kargs


def parse_VandI(line: list[str]) -> Voltage | Current:
    args = {}
    args["name"] = line[0]
    args["nodes"] = [line[1], line[2]]
    _type = {"dc", "ac", "pulse", "sin", "pwl"}
    p = 3
    if line[3] not in _type:
        args["dc"] = line[3]
        p += 1
    else:
        while p < len(line) and line[p] in _type:
            if line[p] == "sin":
                i = 1
                while p + i < len(line) and line[p + i] not in _type:
                    # print(i)
                    i += 1
                if i - 1 < 3:
                    raise ValueError("sin must be 3 args")
                else:
                    args[line[p]] = line[p + 1 : p + i]
                p += i
            elif line[p] == "pulse":
                i = 1
                while p + i < len(line) and line[p + i] not in _type:
                    i += 1
                if i - 1 < 8:
                    raise ValueError("pulse must be 8 args")
                else:
                    args[line[p]] = line[p + 1 : p + i]
                p += i
            elif line[p] == "pwl":
                i = 1
                while p + i < len(line) and line[p + i] not in _type:
                    i += 1
                v_len = i - 1
                if line[p + v_len].startswith("td="):
                    args["td"] = line[p + v_len][3:]
                    v_len -= 1
                if line[p + v_len].startswith("r="):
                    args["r"] = line[p + v_len][2:]
                    v_len -= 1

                if v_len % 2 != 0:
                    raise ValueError("pwl must have an even number of arguments")
                else:
                    args["pwl"] = [
                        (line[p + j], line[p + j + 1]) for j in range(1, v_len, 2)
                    ]
                p += i
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


def parse_subckt(subckt_name, raw_lines, line_name):
    subckt_lines = []  # Initialize an empty list to store lines inside the subcircuit
    while True:
        line = (
            raw_lines.pop(0).strip().lower().split()
        )  # Read the next line and split it into a list
        if not line:
            continue  # Skip empty lines
        if line[0] == ".ends":
            break  # Exit the loop when ".ends" is encountered
        subckt_lines.append(line)  # Add the line to the list of subcircuit lines

    # Write the subcircuit content to a file named after the subcircuit name
    output_file_name = f"{subckt_name}.txt"
    with open(output_file_name, "w") as output_file:
        output_file.write(line_name)
        for subckt_line in subckt_lines:
            output_file.write(
                " ".join(subckt_line) + "\n"
            )  # Write each line of the subcircuit to the file
        output_file.write(
            ".ends\n"
        )  # Write ".ends" to indicate the end of the subcircuit

    # Return the remaining lines after encountering .ends
    return raw_lines


def subckt_read(components, line_main):
    subckt_name = line_main[-1]
    output_file_name = f"{subckt_name}.txt"
    with open(output_file_name, "r") as f:
        raw_lines = f.readlines()
    for line in raw_lines:
        if line[0][0] in {"v", "i"}:
            components.append(parse_VandI(line))
        elif line[0][0] == "d":
            components.append(parse_diode(line))
        elif line[0][0] == "m":
            components.append(parse_mos(line))
        elif line[0][0] == "x":
            subckt_read(components, line)
        else:
            components.append(parse_other(line))
    return


def parse(file_name):
    path = os.path.join("..", "testcase", file_name)
    # 相對位置加檔名
    # with open(file_name, "r") as f:
    # with open(r"C:\Users\USER\Desktop\project_spice\netlist"+ "\\" + file_name, "r") as f:
    with open(path, "r") as f:
        #raw_lines = f.readlines()
        raw_lines=link_program(f.read()).split("\n")
    # 删除第一行
    circuit_name = raw_lines.pop(0).strip()

    lines, model_array = [], []
    # 到".end"为止
    while raw_lines:
        line_long = raw_lines.pop(0).strip()  # all line
        # line = line_subcktname.strip()
        if not line_long:
            continue  # Skip empty lines
        line = line_long.lower().split()
        # 删除注释
        if line_long and line[0].startswith("*"):
            continue
        if line_long and line[0] == ".end":
            break
        if line_long and line[0] == ".model":
            model_array.append(line[1])
            model_array.append(line[2])
        elif line_long and line[0] == ".subckt":
            subckt_name = line[1]
            raw_lines = parse_subckt(subckt_name, raw_lines, line_long)
        else:
            lines.append(line)

    # 創立namelist確定元件及元件數rc
    namelist = list(zip(*lines))[0]
    # print(lines)
    # print(list(zip(*lines)))
    components = []
    for line in lines:
        if line[0][0] in {"v", "i"}:
            components.append(parse_VandI(line))
        elif line[0][0] == "d":
            components.append(parse_diode(line))
        elif line[0][0] == "m":
            components.append(parse_mos(line))
        elif line[0][0] == "x":
            subckt_read(components, line)
        else:
            components.append(parse_other(line))

    return circuit_name, components, namelist


if __name__ == "__main__":
    #file_name = "test.txt"
    file_name = input("file name: ")
    circuit_name, components, namelist = parse(file_name)
    nodelist = [c.nodes for c in components]
    valuelist = [c.value for c in components]
    print(circuit_name)
    print(namelist)
    print(nodelist)
    print(valuelist)
    for component in components:
        print(component)
