from typing import TypeAlias

# 第一個列表是子電路的參數，第二個列表是子電路的內容
Func: TypeAlias = dict[str, tuple[list[str], list[str]]]
# 定義一個字符串，用於後續對內部元件名稱進行後綴的格式化
FORMAT_suffix = "{X_name}"


def add_prefix(line: str, params):
    """將内部元件加上子元件名前綴"""
    # 將傳入的 line 字符串分割成一個列表，取第一個字符
    element = line.split()
    # 多upper
    comp = element[0][0].upper()

    # 如果調用其他子電路，檢查所有參數
    if comp == "X":
        count = len(element) - 2
    # 根據不同的元件類型，確定參數的個數
    # Mos 元件檢查 4 個參數，其他檢查 2 個參數
    elif comp == "M":
        count = 4
    else:
        count = 2

    # 內部元件名稱加上後綴
    # 如果元件不是子電路（X），則對元件名稱進行後綴處理
    if comp != "X":
        element[0] = element[0] + FORMAT_suffix

    # 檢查參數
    # 對參數進行處理，如果不存在於 params 列表中，則進行後綴處理
    for i in range(1, count + 1):
        if element[i] in params or element[i].isdigit():
            continue
        element[i] = element[i] + FORMAT_suffix
    return " ".join(element)


def parse_function(lines: list[str]):
    """分離子電路名稱、參數、內容"""
    # 從第一行獲取子電路名稱和參數
    name, *params = lines[0].split()[1:]
    body = lines[1:-1]

    # 修改內部節點名稱
    for i, b in enumerate(body):
        body[i] = add_prefix(b, params)

    return name, params, body


def parse_program(program: str) -> tuple[Func, list[str]]:
    """分離子電路與主程式"""
    lines = program.strip().split("\n")
    functions = {}
    main_program = []
    current_function: list | None = None  # 保存當前處理的子電路

    for line in lines:
        if line == "":
            continue
        if line.startswith(".subckt"):
            current_function = [line]

        # 如果在子電路中，且遇到 .ends 表示子電路結束
        elif line.startswith(".ends"):
            # 子電路結束，解析子電路並保存到 functions 字典中
            if not isinstance(current_function, list):
                raise SyntaxError("Unexpected .ends")

            current_function.append(line)
            name, params, body = parse_function(current_function)
            functions[name] = (params, body)
            current_function = None

        elif line.startswith(".end"):
            main_program.append(".end")
            break

        # 如果還在子電路中，加入當前行
        elif current_function is not None:
            current_function.append(line)
        # current_function 不存在，表示在主程式中
        else:
            main_program.append(line)
    # 返回子電路字典和主程式列表的元組
    return functions, main_program


def expand_line(line: str, functions: Func, func_name: str = "") -> str:
    # 如果不是子電路，直接返回
    # 原本    if not line.startswith("X"):
    if not line.upper().startswith("X"):
        return line

    line_expand = []
    name, *input_args, func = line[1:].split()
    total_name = f"{func_name}_{name}"
    params, body = functions[func]  # 取得調用子電路參數與內容

    # 將參數替換成輸入參數
    arg_map = dict(zip(params, input_args))
    for body_line in body:
        body_line = body_line.format(X_name=total_name).split()
        for i, arg in enumerate(body_line):
            if arg in arg_map:
                body_line[i] = arg_map[arg]
        body_line = " ".join(body_line)

        # 如果可能，繼續展開子電路
        line_expand.append(expand_line(body_line, functions, total_name))
    # 返回展開後的子電路內容
    return "\n".join(line_expand)


def expand_functions(functions, main_program):
    program_name = main_program.pop(0)  # 第一會是程式名稱
    # 對主程式中的每一行進行展開
    expanded_program = [program_name]
    for line in main_program:
        expanded_program.append(expand_line(line, functions))
    # 返回展開後的主程式內容
    return "\n".join(expanded_program)


def link_program(program):
    # 解析整個程序，包括子電路和主程式
    functions, main_program = parse_program(program)
    return expand_functions(functions, main_program)


if __name__ == "__main__":
    # file_name = input("file name: ")
    file_name = "test.txt"
    with open(file_name, "r") as f:
        program = f.read()
    print(link_program(program))
