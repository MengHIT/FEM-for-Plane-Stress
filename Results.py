import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation


def drawStressAndDisplacement(triangles, stresses, displacements):
    plt.rcParams['font.family'] = 'SimHei'
    plt.rcParams['axes.unicode_minus'] = False
    # 将应力提取为 x、y 和切应力的列表
    x_stress = np.array([s.flatten() for s in stresses])[:, 0]
    y_stress = np.array([s.flatten() for s in stresses])[:, 1]
    shear_stress = np.array([s.flatten() for s in stresses])[:, 2]

    # 准备绘制数据
    def prepare_plot_data(triangles, stress_values):
        points = []
        values = []
        triangle_indices = []

        for i, tri in enumerate(triangles):
            tri_points = [p for p in tri]
            points.extend(tri_points)
            values.extend([stress_values[i]] * len(tri_points))
            triangle_indices.append([len(points) - 3, len(points) - 2, len(points) - 1])

        points = np.array(points)
        values = np.array(values).flatten()
        triangle_indices = np.array(triangle_indices)

        return points, values, triangle_indices

    # 准备绘制每个应力的热力图
    def plot_stress_heatmap(ax, points, values, triangle_indices, title):
        triang = Triangulation(points[:, 0], points[:, 1], triangle_indices)
        heatmap = ax.tricontourf(triang, values, cmap='hot')
        plt.colorbar(heatmap, ax=ax, label='Stress')
        ax.set_title(title)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_aspect('equal')

    # 获取位移后的三角形
    def apply_displacement(triangles, displacements):
        displaced_triangles = []
        for tri in triangles:
            displaced_tri = []
            for point in tri:
                dx, dy = displacements.get(point, [0, 0])
                displaced_tri.append((point[0] + dx, point[1] + dy))
            displaced_triangles.append(displaced_tri)
        return displaced_triangles

    displaced_triangles = apply_displacement(triangles, displacements)

    # 绘制位移图
    def plot_displacement(ax, triangles, displaced_triangles, displacements):
        for tri in triangles:
            tri = np.array(tri + [tri[0]])  # 闭合三角形
            ax.plot(tri[:, 0], tri[:, 1], 'b-', label='Original' if 'Original' not in ax.get_legend_handles_labels()[1] else "")

        for tri in displaced_triangles:
            tri = np.array(tri + [tri[0]])  # 闭合三角形
            ax.plot(tri[:, 0], tri[:, 1], 'r-', label='Displaced' if 'Displaced' not in ax.get_legend_handles_labels()[1] else "")

        for point, disp in displacements.items():
            ax.plot(point[0], point[1], 'bo')  # 原始点
            displaced_point = (point[0] + disp[0], point[1] + disp[1])
            ax.plot(displaced_point[0], displaced_point[1], 'ro')  # 位移点

        ax.legend()
        ax.set_title("位移图")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_aspect('equal')
        ax.grid(True)

    # 准备绘图
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # 绘制 X 应力热力图
    points, values, triangle_indices = prepare_plot_data(triangles, x_stress)
    plot_stress_heatmap(axes[0, 0], points, values, triangle_indices, 'X方向正应力')

    # 绘制 Y 应力热力图
    points, values, triangle_indices = prepare_plot_data(triangles, y_stress)
    plot_stress_heatmap(axes[0, 1], points, values, triangle_indices, 'Y方向正应力')

    # 绘制切应力热力图
    points, values, triangle_indices = prepare_plot_data(triangles, shear_stress)
    plot_stress_heatmap(axes[1, 0], points, values, triangle_indices, '切应力')

    # 绘制位移图
    plot_displacement(axes[1, 1], triangles, displaced_triangles, displacements)

    plt.tight_layout()
    plt.show()


def main():
    triangles = [[(0, 0), (30, 30), (60, 0)], [(60, 0), (90, 30), (30, 30)]]

    # 应力值列表
    stresses =  [np.array([[-0.0400000000000000],
       [0.0266666666666667],
       [0.0400000000000000]], dtype=object),
                 np.array([[-0.0266666666666667],
       [-0.0266666666666667],
       [0.0266666666666667]], dtype=object)]

    displacements = {(0, 0): [2.72000000000000, 0], (30, 30): [4.24000000000000, 1.04000000000000],
                     (60, 0): [0, 0], (90, 30): [2.96000000000000, -2.32000000000000]}

    drawStressAndDisplacement(triangles, stresses, displacements)


if __name__ == '__main__':
    main()