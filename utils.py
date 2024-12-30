import sympy as sp
import math


def isInTriangle(A, B, C, P):
    def calArea(A0, B0, C0):
        c = math.dist(A0, B0)
        b = math.dist(A0, C0)
        a = math.dist(B0, C0)
        s = (a + b + c) / 2
        return math.sqrt(s * (s - a) * (s - b) * (s - c))

    area0 = calArea(A, B, C)
    area1 = calArea(A, B, P)
    area2 = calArea(A, P, C)
    area3 = calArea(P, B, C)
    return math.fabs(area1 + area2 + area3 - area0) < 0.1


def segmentInTriangle(A, B, C, m, n):  # 判断线段是否在三角形内
    if isInTriangle(A, B, C, m) and isInTriangle(A, B, C, n):
        return True
    else:
        return False


def triangleInTriangle(A, B, C, m, n, l):
    if isInTriangle(A, B, C, m) and isInTriangle(A, B, C, n) and isInTriangle(A, B, C, l):
        return True
    else:
        return False


def prTriangle(A, B, C):  # 判断三角形是否有两顶点在同一条竖线上
    if A[0] == B[0]:
        return C
    elif A[0] == C[0]:
        return B
    elif C[0] == B[0]:
        return A
    else:
        numbers = sorted([A[0], B[0], C[0]])
        mid_0 = numbers[1]
        add_p = [mid_0]
        if mid_0 == A[0]:
            add_p.append(B[1] + (C[1] - B[1]) * (mid_0 - B[0]) / (C[0] - B[0]))
        elif mid_0 == B[0]:
            add_p.append(A[1] + (C[1] - A[1]) * (mid_0 - A[0]) / (C[0] - A[0]))
        elif mid_0 == C[0]:
            add_p.append(B[1] + (A[1] - B[1]) * (mid_0 - B[0]) / (A[0] - B[0]))
        add_p_t = tuple(add_p)
        return add_p_t


def integrateArea(A, B, C):
    x, y = sp.symbols('x y')
    third_point = prTriangle(A, B, C)
    if third_point == A:
        x_range = [min(A[0], C[0]), max(A[0], C[0])]
        min_y, max_y = min(B[1], C[1]), max(B[1], C[1])
        y_range = [A[1] + (max_y - A[1]) * (x - A[0]) / (B[0] - A[0]),
                   A[1] + (min_y - A[1]) * (x - A[0]) / (B[0] - A[0])]
        return x_range, y_range, 1
    elif third_point == B:
        x_range = [min(B[0], C[0]), max(B[0], C[0])]
        min_y, max_y = min(A[1], C[1]), max(A[1], C[1])
        y_range = [B[1] + (max_y - B[1]) * (x - B[0]) / (C[0] - B[0]),
                   B[1] + (min_y - B[1]) * (x - B[0]) / (C[0] - B[0])]
        return x_range, y_range, 1
    elif third_point == C:
        x_range = [min(A[0], C[0]), max(A[0], C[0])]
        min_y, max_y = min(B[1], A[1]), max(B[1], A[1])
        y_range = [C[1] + (max_y - C[1]) * (x - C[0]) / (A[0] - C[0]),
                   C[1] + (min_y - C[1]) * (x - C[0]) / (A[0] - C[0])]
        return x_range, y_range, 1
    else:
        if third_point[0] == A[0]:
            x_range_1 = [min(third_point[0], B[0]), max(third_point[0], B[0])]
            x_range_2 = [min(third_point[0], C[0]), max(third_point[0], C[0])]
            min_y, max_y = min(third_point[1], A[1]), max(third_point[1], A[1])
            y_range_1 = [B[1] + (max_y - B[1]) * (x - B[0]) / (A[0] - B[0]),
                         B[1] + (min_y - B[1]) * (x - B[0]) / (A[0] - B[0])]
            y_range_2 = [C[1] + (max_y - C[1]) * (x - C[0]) / (A[0] - C[0]),
                         C[1] + (min_y - C[1]) * (x - C[0]) / (A[0] - C[0])]
            return [x_range_1, y_range_1], [x_range_2, y_range_2], 0
        elif third_point[0] == B[0]:
            x_range_1 = [min(third_point[0], A[0]), max(third_point[0], A[0])]
            x_range_2 = [min(third_point[0], C[0]), max(third_point[0], C[0])]
            min_y, max_y = min(third_point[1], B[1]), max(third_point[1], B[1])
            y_range_1 = [A[1] + (max_y - A[1]) * (x - A[0]) / (B[0] - A[0]),
                         A[1] + (min_y - A[1]) * (x - A[0]) / (B[0] - A[0])]
            y_range_2 = [C[1] + (max_y - C[1]) * (x - C[0]) / (B[0] - C[0]),
                         C[1] + (min_y - C[1]) * (x - C[0]) / (B[0] - C[0])]
            return [x_range_1, y_range_1], [x_range_2, y_range_2], 0
        elif third_point[0] == C[0]:
            x_range_1 = [min(third_point[0], B[0]), max(third_point[0], B[0])]
            x_range_2 = [min(third_point[0], A[0]), max(third_point[0], A[0])]
            min_y, max_y = min(third_point[1], C[1]), max(third_point[1], C[1])
            y_range_1 = [B[1] + (max_y - B[1]) * (x - B[0]) / (C[0] - B[0]),
                         B[1] + (min_y - B[1]) * (x - B[0]) / (C[0] - B[0])]
            y_range_2 = [A[1] + (max_y - A[1]) * (x - A[0]) / (C[0] - A[0]),
                         A[1] + (min_y - A[1]) * (x - A[0]) / (C[0] - A[0])]
            return [x_range_1, y_range_1], [x_range_2, y_range_2], 0


def calX_num(ass):
    sum_u, sum_v = 0, 0
    u = []
    v = []
    for val in list(ass.values()):
        if val == 1:
            sum_v += 1
        elif val == 2:
            sum_u += 1
        elif val == 0:
            sum_u += 1
            sum_v += 1
    for i in range(0, sum_u):
        u.append(sp.symbols(f'u{i}'))
    for i in range(0, sum_v):
        v.append(sp.symbols(f'v{i}'))

    cu, cv = 0, 0
    key_list = list(ass.keys())
    val_list = list(ass.values())
    for i in range(0, len(val_list)):
        if val_list[i] == 1:
            val_list[i] = [0, v[cv]]
            cv += 1
        elif val_list[i] == 2:
            val_list[i] = [u[cu], 0]
            cu += 1
        elif val_list[i] == 3:
            val_list[i] = [0, 0]
        elif val_list[i] == 0:
            val_list[i] = [u[cu], v[cv]]
            cu += 1
            cv += 1
    ass_0 = dict(zip(key_list, val_list))
    return ass_0


def main():
    A = tuple((2, 1))
    B = tuple((0, 1))
    C = tuple((1, 0))
    P = tuple((0.333333, 1))
    # print(integrateArea(A, B, C))
    ass = {(360, 300): 2, (450, 300): 2, (390, 240): 0, (480, 240): 0, (420, 180): 3}
    ass = calX_num(ass)
    # print(ass)


if __name__ == '__main__':
    main()
