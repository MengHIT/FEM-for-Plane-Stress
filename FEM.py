import numpy as np
import sympy as sp
import utils


class UnitFEM:
    def __init__(self, co_0, Elast, v, thick, F_assemble):
        self.co_0 = co_0
        self.E = Elast
        self.v = v
        self.thick = thick

        self.x, self.y = sp.symbols('x y')
        self.F_assemble = F_assemble

    def findAnticlockwise(self):  # 逆时针标号 co为三个坐标
        co = self.co_0.copy()
        co_x = [co[0][0], co[1][0], co[2][0]]
        minx = co_x.index(min(co_x))
        i = co[minx]  # x最小的
        co.remove(i)

        co_x2 = [co[0][0], co[1][0]]
        co_y = [co[0][1], co[1][1]]
        if co_y[1] == co_y[0]:
            j = co[co_x2.index(min(co_x2))]
        else:
            miny = co_y.index(min(co_y))
            j = co[miny]
        co.remove(j)
        m = co[0]

        return [i, j, m]

    def shapeFunction(self, co, x, y):  # 形函数N
        i, j, m = co[0], co[1], co[2]
        a, b, c = [], [], []
        N = []
        for count in range(0, 3):
            a.append(co[(count + 1) % 3][0] * co[(count + 2) % 3][1] - co[(count + 2) % 3][0] * co[(count + 1) % 3][1])
            b.append(co[(count + 1) % 3][1] - co[(count + 2) % 3][1])
            c.append(co[(count + 2) % 3][0] - co[(count + 1) % 3][0])
        # 面积
        Area = 0.5 * (a[0] + i[0] * b[0] + i[1] * c[0])
        for count in range(0, 3):
            N.append((a[count] + b[count] * x + c[count] * y) / (2 * Area))

        N_array = np.array([
            [N[0], 0, N[1], 0, N[2], 0],
            [0, N[0], 0, N[1], 0, N[2]],
        ])

        return a, b, c, Area, N, N_array

    def gradFunction(self, b, c, Area):  # 梯度矩阵B
        B = []
        for count in range(0, 3):
            B.append(np.array(
                [[b[count], 0],
                 [0, c[count]],
                 [c[count], b[count]]]) / (2 * Area))

        return B

    def stressConvertFunction(self, b, c, Area):  # 应力转换矩阵S
        S = []
        for count in range(0, 3):
            S.append(np.array(
                [[b[count], self.v * c[count]],
                 [self.v * b[count], c[count]],
                 [(0.5 - self.v / 2) * c[count], (0.5 - self.v / 2) * b[count]]])
                     * self.E / (2 * (1 - self.v * self.v) * Area))
        S_array = np.array([
            [b[0], self.v * c[0], b[1], self.v * c[1], b[2], self.v * c[2]],
            [self.v * b[0], c[0], self.v * b[1], c[1], self.v * b[2], c[2]],
            [(0.5 - self.v / 2) * c[0], (0.5 - self.v / 2) * b[0], (0.5 - self.v / 2) * c[1], (0.5 - self.v / 2) * b[1], (0.5 - self.v / 2) * c[2], (0.5 - self.v / 2) * b[2]]
        ]) * self.E / (2 * (1 - self.v * self.v) * Area)

        return S, S_array

    def stiffnessFunction(self, B, S, Area):  # 单元劲度系数矩阵k
        k_list = []
        for i in range(0, 3):
            for j in range(0, 3):
                k_list.append(np.dot(B[i].T, S[j]) * self.thick * Area)

        return k_list

    def calUnit_k(self):
        co_anti = self.findAnticlockwise()
        a, b, c, Area, N, N_array = self.shapeFunction(co_anti, self.x, self.y)
        B = self.gradFunction(b, c, Area)
        S, S_array = self.stressConvertFunction(b, c, Area)
        k = self.stiffnessFunction(B, S, Area)

        # 力的转置
        f_dic = {tuple(co_anti[0]): [0, 0], tuple(co_anti[1]): [0, 0], tuple(co_anti[2]): [0, 0]}
        for i in range(0, len(self.F_assemble[0])):  # 判断集中力
            f_place = list(self.F_assemble[0].keys())
            if utils.isInTriangle(tuple(co_anti[0]), tuple(co_anti[1]), tuple(co_anti[2]), f_place[i]):
                f_mea = np.array(list(self.F_assemble[0].values())[i])
                f = np.dot(N_array.T, f_mea.T) * self.thick
                for j in range(0, len(f)):
                    f[j] = f[j].subs([(self.x, f_place[i][0]), (self.y, f_place[i][1])])
                f_dic[tuple(co_anti[0])][0] += f[0]
                f_dic[tuple(co_anti[0])][1] += f[1]
                f_dic[tuple(co_anti[1])][0] += f[2]
                f_dic[tuple(co_anti[1])][1] += f[3]
                f_dic[tuple(co_anti[2])][0] += f[4]
                f_dic[tuple(co_anti[2])][1] += f[5]

        for i in range(0, len(self.F_assemble[1])):    # 判断面力
            f_place = list(self.F_assemble[1].keys())
            if utils.segmentInTriangle(tuple(co_anti[0]), tuple(co_anti[1]), tuple(co_anti[2]), f_place[i][0], f_place[i][1]):
                f_mea = np.array(list(self.F_assemble[1].values())[i])
                if f_place[i][1][0]-f_place[i][0][0] != 0:   # 斜率存在
                    dif = (f_place[i][1][1]-f_place[i][0][1]) / (f_place[i][1][0]-f_place[i][0][0])
                    expr = sp.Eq(self.y, dif * self.x + f_place[i][1][1] - dif * f_place[i][1][0])
                    f_sol = []
                    f = np.dot(N_array.T, f_mea.T) * self.thick * np.sqrt(dif * dif + 1)
                    if isinstance(sp.solve(expr)[0], dict) is False:
                        sp_f = sp.Matrix(f).subs({self.y: sp.solve(expr)[0]})
                        for ele in sp_f:
                            temp_f = sp.integrate(ele, (self.x, f_place[i][0][0], f_place[i][1][0]))
                            f_sol.append(temp_f)
                    else:
                        sp_f = sp.Matrix(f).subs({self.x: list(sp.solve(expr)[0].values())[0]})
                        for ele in sp_f:
                            temp_f = sp.integrate(ele, (self.y, f_place[i][0][1], f_place[i][1][1]))
                            f_sol.append(temp_f)

                else:     # 斜率不存在
                    expr = sp.Eq(self.x, f_place[i][1][0])
                    f = np.dot(N_array.T, f_mea.T) * self.thick
                    sp_f = sp.Matrix(f).subs({self.x: sp.solve(expr)[0]})
                    f_sol = []
                    for ele in sp_f:
                        temp_f = sp.integrate(ele, (self.y, f_place[i][0][1], f_place[i][1][1]))
                        f_sol.append(temp_f)
                f_dic[tuple(co_anti[0])][0] += f_sol[0]
                f_dic[tuple(co_anti[0])][1] += f_sol[1]
                f_dic[tuple(co_anti[1])][0] += f_sol[2]
                f_dic[tuple(co_anti[1])][1] += f_sol[3]
                f_dic[tuple(co_anti[2])][0] += f_sol[4]
                f_dic[tuple(co_anti[2])][1] += f_sol[5]

        for i in range(0, len(self.F_assemble[2])):
            f_place = list(self.F_assemble[2].keys())
            if utils.triangleInTriangle(tuple(co_anti[0]), tuple(co_anti[1]), tuple(co_anti[2]), f_place[i][0], f_place[i][1], f_place[i][2]):
                f_mea = np.array(list(self.F_assemble[2].values())[i])
                f = np.dot(N_array.T, f_mea.T) * self.thick
                range_1, range_2, pr = utils.integrateArea(tuple(co_anti[0]), tuple(co_anti[1]), tuple(co_anti[2]))
                f_integrate = []
                if pr == 1:
                    for ele in f:
                        temp_f = sp.integrate(ele, (self.y, range_2[0], range_2[1]), (self.x, range_1[0], range_1[1]))
                        f_integrate.append(temp_f)
                else:
                    for ele in f:
                        temp_f = sp.integrate(ele, (self.y, range_1[1][0], range_1[1][1]), (self.x, range_1[0][0], range_1[0][1])) \
                                 + sp.integrate(ele, (self.y, range_2[1][0], range_2[1][1]), (self.x, range_2[0][0], range_2[0][1]))
                        f_integrate.append(temp_f)
                f_dic[tuple(co_anti[0])][0] += f_integrate[0]
                f_dic[tuple(co_anti[0])][1] += f_integrate[1]
                f_dic[tuple(co_anti[1])][0] += f_integrate[2]
                f_dic[tuple(co_anti[1])][1] += f_integrate[3]
                f_dic[tuple(co_anti[2])][0] += f_integrate[4]
                f_dic[tuple(co_anti[2])][1] += f_integrate[5]

        return k, f_dic, S_array


class FEM:
    def __init__(self, co_assemble, Elast, v, thick, F_assemble, delta):
        self.co_assemble = co_assemble
        self.Elast = Elast
        self.v = v
        self.thick = thick
        self.F_assemble = F_assemble
        self.delta = delta

    def FEM_K(self):  # 合并所有三角形的角点，以索引作为标号，并把各个三角形的角点标号
        co_xj, co_xjj, co_anti_list = [], [], []
        co_copy = self.co_assemble.copy()

        for i in range(0, len(self.co_assemble)):
            co_xj += self.co_assemble[i]
        [co_xjj.append(i) for i in co_xj if i not in co_xjj]  # co_xjj所有角点，以其索引为标号

        indexes, index_corner = [], []
        for i in range(0, len(self.co_assemble)):
            Unit_FEM_1 = UnitFEM(self.co_assemble[i], self.Elast, self.v, self.thick, self.F_assemble)
            co_anti = Unit_FEM_1.findAnticlockwise()
            co_anti_list.append(co_anti)
            for j in range(0, 3):
                indexes += [index for index, value in enumerate(co_xjj) if value == co_anti[j]]
                index_corner = [indexes[i:i + 3] for i in range(0, len(indexes), 3)]  # [i, j, m]在整体中的标号

        whole_K = np.zeros((len(co_xjj) * 2, len(co_xjj) * 2))  # 整体劲度矩阵赋0
        whole_F = np.zeros(len(co_xjj) * 2)
        delta_list = [0] * len(co_xjj) * 2

        for iii in range(0, len(co_xjj)):
            for (xxx, yyy), value_2 in self.delta.items():
                if (xxx, yyy) == tuple(co_xjj[iii]):
                    delta_list[2 * iii] += value_2[0]
                    delta_list[2 * iii + 1] += value_2[1]
        delta_array = np.array(delta_list)

        S_total = []
        for i in range(0, len(index_corner)):
            key_in = [
                (index_corner[i][0], index_corner[i][0]), (index_corner[i][0], index_corner[i][1]),
                (index_corner[i][0], index_corner[i][2]),
                (index_corner[i][1], index_corner[i][0]), (index_corner[i][1], index_corner[i][1]),
                (index_corner[i][1], index_corner[i][2]),
                (index_corner[i][2], index_corner[i][0]), (index_corner[i][2], index_corner[i][1]),
                (index_corner[i][2], index_corner[i][2])
            ]

            Unit_FEM_2 = UnitFEM(co_copy[i], self.Elast, self.v, self.thick, self.F_assemble)
            k_aa, f_dic, S_a = Unit_FEM_2.calUnit_k()
            S_total.append(S_a)

            for ii in range(0, len(co_xjj)):
                for (xx, yy), value_1 in f_dic.items():
                    if (xx, yy) == tuple(co_xjj[ii]):
                        whole_F[2 * ii] += value_1[0]
                        whole_F[2 * ii + 1] += value_1[1]

            index2k = dict(zip(key_in, k_aa))
            for (row, col), k_ijm in index2k.items():
                whole_K[row * 2, col * 2] += k_ijm[0][0]
                whole_K[row * 2, col * 2 + 1] += k_ijm[0][1]
                whole_K[row * 2 + 1, col * 2] += k_ijm[1][0]
                whole_K[row * 2 + 1, col * 2 + 1] += k_ijm[1][1]

        return whole_K, whole_F, delta_array, S_total, co_anti_list   # whole_K 整体劲度矩阵

    def solveFEM(self):
        whole_K, whole_F, delta_array, S_total, co_anti_list = self.FEM_K()
        cal_F = np.dot(whole_K, delta_array.T)
        retain_delta_list = []
        retain_list = []
        for i in range(0, len(delta_array)):
            if delta_array[i] != 0:
                retain_list.append(i)
                retain_delta_list.append(delta_array[i])

        retain_delta = np.array(retain_delta_list)

        retain_K = np.zeros((len(retain_list), len(retain_list)))
        retain_F = np.zeros(len(retain_list))

        for i in range(0, len(retain_list)):
            retain_F[i] = whole_F[retain_list[i]]

        for i in range(0, len(retain_list)):
            for j in range(0, len(retain_list)):
                retain_K[i][j] = whole_K[retain_list[i]][retain_list[j]]

        solve_0 = np.dot(retain_K, retain_delta.T)
        sp_solve_0 = sp.Matrix(solve_0)
        sp_retain_F = sp.Matrix(retain_F)
        sp_v = sp.Matrix(retain_delta)
        equation = sp.Eq(sp_solve_0, sp_retain_F)
        solution_0 = sp.solve(equation, sp_v)

        sp_cal_F = sp.Matrix(cal_F)
        solution_F = np.array(sp_cal_F.subs(solution_0)).T[0] - whole_F

        solution_delta = self.delta.copy()
        temp_s = []
        temp_ke = []
        for i in range(0, len(solution_delta)):
            temp_s.append(list(sp.Matrix(list(solution_delta.values())[i]).subs(solution_0)))
            temp_ke.append(list(solution_delta.keys())[i])
        s_delta = dict(zip(temp_ke, temp_s))

        # 求各单元应力
        stress_assemble = []
        ttt = 0
        for co_ev in co_anti_list:
            s_tt = []
            for i in range(0, 3):
                for j in range(0, len(s_delta)):
                    if tuple(co_ev[i]) == list(s_delta.keys())[j]:
                        s_tt.append(list(s_delta.values())[j])
            stress_assemble.append(np.dot(S_total[ttt], np.array(s_tt).reshape(1, -1).T))
            ttt += 1

        return solution_F, s_delta, stress_assemble


def main():
    x, y = sp.symbols('x y')
    co_assemble = [
        [[1, 0], [2, 1], [0, 1]],
        [[1, 0], [0, 0], [0, 1]],
        [[0, 0], [1, 0], [1, -1]]
    ]

    F_assemble = [
        {(1/3, 1): [0, -0.5]},           # 集中力
        {((0, 1), (1, 1)): [0, x-y]},           # 面力 起止点的坐标，q的函数,必须作用在单元边界上
        {((0, 1), (2, 1), (1, 0)): [0, 1]}            # 体力 作用的单元标号
    ]

    v1, v2 = sp.symbols('v1 v2')
    delta = {
        (0, 0): [0, v1], (1, 0): [0, 0], (2, 1): [0, 0], (0, 1): [0, v2], (1, -1): [0, 0]
    }
    Elast = 35 / 36
    v = 1 / 6
    thick = 2

    FEMxx = FEM(co_assemble,  Elast, v, thick, F_assemble, delta)
    solution_F, s_delta, stress_assemble = FEMxx.solveFEM()   # 结点力， 位移， 单元应力
    print(solution_F, s_delta, stress_assemble)


if __name__ == '__main__':
    main()
