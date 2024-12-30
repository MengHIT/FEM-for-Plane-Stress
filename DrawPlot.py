import tkinter as tk
from tkinter import messagebox, simpledialog

from numpy import array
from sympy import symbols
import utils
import FEM
import Results


class TriangleDrawer:
    def __init__(self, root, grid_size=30, snap_radius=10):
        self.x, self.y = symbols('x y')

        self.root = root
        self.grid_size = grid_size
        self.snap_radius = snap_radius
        self.origin_x = 360  # 坐标原点 x
        self.origin_y = 300  # 坐标原点 y
        self.points = []
        self.all_triangles = []  # 保存所有三角形的顶点
        self.all_triangles_xx = []
        self.all_triangles_yy = []
        self.point_types = {}  # 记录约束类型，格式为 {(x, y): type}
        self.point_types_xx = {}
        self.inserted_points = {}  # 记录集中力及其数组，格式为 {(x, y): [x_value, y_value]}
        self.edge_arrays = {}  # 记录面力及其数组，格式为 {((x1, y1), (x2, y2)): [x, y]}
        self.triangle_arrays = {}  # 记录体力及其数组，格式为 {((x1, y1), (x2, y2), (x3, y3)): [x, y]}

        self.selected_type = None  # 当前选择的类型
        self.active_feature = None  # 当前激活的功能 (1, 2, 3, or None)

        # 主框架布局
        self.frame = tk.Frame(root)
        self.frame.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # 画布
        self.canvas = tk.Canvas(self.frame, width=800, height=600, bg="white")
        self.canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # 按钮
        self.buttons_frame = tk.Frame(root)
        self.buttons_frame.pack(side=tk.RIGHT, fill=tk.Y)

        self.create_buttons()

        self.option_frame = None  # 存储选项框架
        self.draw_grid()
        self.canvas.bind("<Button-1>", self.on_click)

        self.param1 = None
        self.param2 = None
        self.param3 = None

    def create_buttons(self):
        # 按钮：结束绘制
        btn_end = tk.Button(
            self.buttons_frame, text="结束绘制", command=self.finish_drawing
        )
        btn_end.pack(fill=tk.X)

        btn_feature1 = tk.Button(self.buttons_frame, text="添加约束", command=lambda: self.toggle_feature(1))
        btn_feature1.pack(fill=tk.X)

        btn_feature2 = tk.Button(self.buttons_frame, text="集中力", command=lambda: self.toggle_feature(2))
        btn_feature2.pack(fill=tk.X)

        btn_feature3 = tk.Button(self.buttons_frame, text="面力", command=lambda: self.toggle_feature(3))
        btn_feature3.pack(fill=tk.X)

        btn_feature4 = tk.Button(self.buttons_frame, text="体力", command=lambda: self.toggle_feature(4))
        btn_feature4.pack(fill=tk.X)

        # 按钮：清空画布
        btn_clear = tk.Button(
            self.buttons_frame, text="清空画布", command=self.clear_canvas
        )
        btn_clear.pack(fill=tk.X)

        btn_input_params = tk.Button(
            self.buttons_frame, text="输入参数", command=self.input_parameters
        )
        btn_input_params.pack(fill=tk.X)

    def input_parameters(self):
        try:
            # 获取三个参数
            self.param1 = simpledialog.askstring("输入弹性模量", "请输入弹性模量：")
            self.param2 = simpledialog.askstring("输入泊松比", "请输入泊松比：")
            self.param3 = simpledialog.askstring("输入厚度", "请输入厚度：")

            # 检查输入是否为空
            if self.param1 is None or self.param2 is None or self.param3 is None:
                raise ValueError("所有参数都必须输入！")

            # 显示输入结果
            messagebox.showinfo("参数输入结果", f"弹性模量: {self.param1}\n泊松比: {self.param2}\n厚度: {self.param3}")
            print(f"用户输入的参数：\n弹性模量: {self.param1}\n泊松比: {self.param2}\n厚度: {self.param3}")
        except ValueError as e:
            messagebox.showerror("错误", str(e))

    def draw_grid(self):
        width = int(self.canvas["width"])
        height = int(self.canvas["height"])

        for x in range(0, width, self.grid_size):
            self.canvas.create_line(x, 0, x, height, fill="lightgray")
        for y in range(0, height, self.grid_size):
            self.canvas.create_line(0, y, width, y, fill="lightgray")

            # 绘制坐标原点
        origin_screen_x = self.origin_x
        origin_screen_y = self.origin_y
        self.canvas.create_oval(
            origin_screen_x - 5, origin_screen_y - 5,
            origin_screen_x + 5, origin_screen_y + 5,
            fill="red"
        )
        self.canvas.create_text(
            origin_screen_x + 10, origin_screen_y - 10,
             text="原点", fill="black"
        )

    def snap_to_grid(self, x, y):
        grid_x = round((x - self.origin_x) / self.grid_size) * self.grid_size + self.origin_x
        grid_y = round((y - self.origin_y) / self.grid_size) * self.grid_size + self.origin_y
        if abs(grid_x - x) <= self.snap_radius:
            x = grid_x
        if abs(grid_y - y) <= self.snap_radius:
            y = grid_y
        return x, y

    def toggle_feature(self, feature):
        if self.active_feature == feature:
            self.active_feature = None
            messagebox.showinfo("提示", f"已退出 {feature} 模式。")
        else:
            self.active_feature = feature
            if feature == 1:
                self.activate_feature1()
            elif feature == 2:
                self.activate_feature2()
            elif feature == 3:
                self.activate_feature3()
            elif feature == 4:
                self.activate_feature4()

    def activate_feature1(self):
        if self.option_frame:
            self.option_frame.destroy()
        self.option_frame = tk.Frame(self.frame, bg="lightgray")
        self.option_frame.pack(side=tk.TOP, fill=tk.X)
        tk.Button(
            self.option_frame, text="横向链杆", bg="green",
            command=lambda: self.set_type(1)
        ).pack(side=tk.LEFT, padx=5, pady=5)
        tk.Button(
            self.option_frame, text="竖向链杆", bg="yellow",
            command=lambda: self.set_type(2)
        ).pack(side=tk.LEFT, padx=5, pady=5)
        tk.Button(
            self.option_frame, text="铰支", bg="purple",
            command=lambda: self.set_type(3)
        ).pack(side=tk.LEFT, padx=5, pady=5)
        messagebox.showinfo("提示", "已进入添加约束模式，请选择类型并点击画布。")

    def activate_feature2(self):
        messagebox.showinfo("提示", "已进入结点力模式，点击画布以插入点并输入数据。")

    def activate_feature3(self):
        messagebox.showinfo("提示", "已进入面力模式，请选择三角形的一条边。")

    def activate_feature4(self):
        messagebox.showinfo("提示", "已进入体力模式，请点击一个三角形")

    def on_click(self, event):
        x, y = self.snap_to_grid(event.x, event.y)
        if self.active_feature == 1:
            if self.selected_type is not None:
                self.point_types[(x, y)] = self.selected_type
                self.point_types_xx[(x - self.origin_x, self.origin_y - y)] = self.selected_type
                color = {1: "green", 2: "yellow", 3: "purple"}[self.selected_type]
                self.canvas.create_oval(x - 5, y - 5, x + 5, y + 5, outline=color, width=2)
        elif self.active_feature == 2:
            self.insert_point(x, y)
        elif self.active_feature == 3:
            self.select_edge(x, y, self.x, self.y)
        elif self.active_feature == 4:
            self.select_triangle(x, y)
        else:
            self.points.append((x, y))
            self.canvas.create_oval(x - 3, y - 3, x + 3, y + 3, fill="red")
            if len(self.points) == 3:
                self.draw_triangle()

    def select_triangle(self, x, y):
        nearest_triangle = None
        min_dist = float("inf")
        for triangle in self.all_triangles:
            centroid = tuple(sum(coord) / 3 for coord in zip(*triangle))
            dist = ((centroid[0] - x) ** 2 + (centroid[1] - y) ** 2) ** 0.5
            if dist < min_dist:
                min_dist = dist
                nearest_triangle = triangle
        if nearest_triangle and min_dist <= self.snap_radius * 2:
            self.highlight_triangle(nearest_triangle)
            self.add_triangle_array(nearest_triangle, self.x, self.y)
        else:
            messagebox.showerror("提示", "未找到选中的三角形，请点击靠近三角形中心的位置！")

    def highlight_triangle(self, triangle):
        # 添加斜线阴影效果
        x1, y1 = triangle[0]
        x2, y2 = triangle[1]
        x3, y3 = triangle[2]
        self.canvas.create_polygon(triangle, outline="blue", fill="", width=2)
        for i in range(0, 21, 2):  # 绘制斜线阴影
            self.canvas.create_line(x1 + i, y1 + i, x2 + i, y2 + i, fill="gray")

    def add_triangle_array(self, triangle, x, y):
        x1, y1 = triangle[0]
        x2, y2 = triangle[1]
        x3, y3 = triangle[2]
        x1x, x2x, x3x = x1 - self.origin_x, x2 - self.origin_x, x3 - self.origin_x
        y1x, y2x, y3x = - y1 + self.origin_y, - y2 + self.origin_y, - y3 + self.origin_y
        # 为选中的三角形输入 1×2 数组
        array_input = simpledialog.askstring(
            "设置力的大小", f"为三角形 {triangle} 输入 1×2 数组（格式：[x, y]）:"
        )
        try:
            array = eval(array_input)
            if isinstance(array, list) and len(array) == 2:
                triangle_key = tuple(sorted(triangle))
                triangle_key_xx = ((x1x, y1x), (x2x, y2x), (x3x, y3x))
                self.triangle_arrays[triangle_key_xx] = array
                print(f"添加三角形: {triangle_key_xx}, 数组: {array}")
            else:
                raise ValueError
        except Exception:
            messagebox.showerror("错误", "输入格式错误，请输入类似 [0, -0.5] 的数组！")

    def insert_point(self, x, y):
        array_input = simpledialog.askstring("设置力的大小", "请输入 1×2 数组（格式：[x, y]）:")
        try:
            array = eval(array_input)
            if isinstance(array, list) and len(array) == 2:
                fractional_x = int(f"{x}")
                fractional_y = int(f"{y}")
                self.inserted_points[(fractional_x - self.origin_x, - fractional_y + self.origin_y)] = array
                self.canvas.create_oval(x - 6, y - 6, x + 6, y + 6, outline="black", width=2)
                print(f"插入点: ({fractional_x - self.origin_x}, {- fractional_y + self.origin_y}) -> {array}")
            else:
                raise ValueError
        except Exception:
            messagebox.showerror("错误", "输入格式错误，请输入类似 [0, -0.5] 的数组！")

    def select_edge(self, x1, y1, x, y):
        nearest_edge = None
        nearest_edge_xx = None
        min_dist = float("inf")
        for triangle in self.all_triangles:
            for i in range(3):
                p1, p2 = triangle[i], triangle[(i + 1) % 3]
                mid_x, mid_y = (p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2
                dist = ((mid_x - x1) ** 2 + (mid_y - y1) ** 2) ** 0.5
                if dist < min_dist:
                    min_dist = dist
                    nearest_edge = (p1, p2)
                    nearest_edge_xx = ((p1[0] - self.origin_x, - p1[1] + self.origin_y), (p2[0] - self.origin_x, - p2[1] + self.origin_y))
        if nearest_edge and min_dist <= self.snap_radius:
            self.canvas.create_line(*nearest_edge[0], *nearest_edge[1], fill="red", width=3)
            array_input = simpledialog.askstring("设置力的大小", f"为边 {nearest_edge} 输入 1×2 数组:")
            try:
                array = eval(array_input)
                if isinstance(array, list) and len(array) == 2:
                    self.edge_arrays[nearest_edge_xx] = array
                    print(f"边: {nearest_edge_xx} -> {array}")
                else:
                    raise ValueError
            except Exception:
                messagebox.showerror("错误", "输入格式错误，请输入类似 [0, -0.5] 的数组！")

    def set_type(self, type_id):
        self.selected_type = type_id
        print(f"Selected type: {type_id}")

    def draw_triangle(self):
        self.canvas.create_polygon(self.points, outline="blue", fill="", width=2)
        self.all_triangles.append(self.points)
        self.all_triangles_xx.append([[self.points[0][0] - self.origin_x, self.origin_y - self.points[0][1]],
                                      [self.points[1][0] - self.origin_x, self.origin_y - self.points[1][1]],
                                      [self.points[2][0] - self.origin_x, self.origin_y - self.points[2][1]]])
        self.all_triangles_yy.append([(self.points[0][0] - self.origin_x, self.origin_y - self.points[0][1]),
                                      (self.points[1][0] - self.origin_x, self.origin_y - self.points[1][1]),
                                      (self.points[2][0] - self.origin_x, self.origin_y - self.points[2][1])])
        for point in self.points:
            self.point_types[point] = 0
            self.point_types_xx[point[0] - self.origin_x, self.origin_y - point[1]] = 0
        self.points = []

    def finish_drawing(self):
        print("全部三角形：", self.all_triangles)
        print("约束及类型：", self.point_types)
        print("结点力：", self.inserted_points)
        print("面力：", self.edge_arrays)
        print("体力：", self.triangle_arrays)
        messagebox.showinfo(
            "绘制完成",
            f"约束及类型:\n{self.point_types}\n\n结点力:\n{self.inserted_points}\n\n面力:\n{self.edge_arrays}\n\n体力:\n{self.triangle_arrays} "
        )

    def clear_canvas(self):
        self.canvas.delete("all")
        self.points = []
        self.all_triangles = []
        self.point_types = {}
        self.inserted_points = {}
        self.edge_arrays = {}
        self.active_feature = None
        self.draw_grid()

        messagebox.showinfo("提示", "画布已清空！")


if __name__ == "__main__":
    root = tk.Tk()
    root.title("绘制工具")
    app = TriangleDrawer(root)
    root.mainloop()
    co_assemble = app.all_triangles_xx
    co_assemble_tuple = app.all_triangles_yy
    delta = utils.calX_num(app.point_types_xx)
    F_assemble = [app.inserted_points, app.edge_arrays, app.triangle_arrays]

    Elast = float(app.param1)
    v = float(app.param2)
    thick = float(app.param3)
    FEMxx = FEM.FEM(co_assemble,  Elast, v, thick, F_assemble, delta)
    solution_F, s_delta, stress_assemble = FEMxx.solveFEM()   # 结点力， 位移， 单元应力
    for i in range(0, len(stress_assemble)):
        stress_assemble[i] = array(stress_assemble[i], dtype=float)
    Results.drawStressAndDisplacement(co_assemble_tuple, stress_assemble, s_delta)
