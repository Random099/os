import numpy as np
import matplotlib.pyplot as plt
import math
import copy


class DE:
    def __init__(self, initial_cond, boundary_cond_1, boundary_cond_2):
        self.initial_cond = initial_cond
        self.boundary_cond_1 = boundary_cond_1
        self.boundary_cond_2 = boundary_cond_2

    def euler_method(self, t: int, x: int, t_step: float, x_step: float, t_start: int, x_start: int):
        t_steps = int((t/t_step)+1)
        x_steps = int((x/x_step)+1)

        v = np.zeros(x_steps)
        if isinstance(self.initial_cond, (float, int)):
            u = np.full(x_steps, float(self.initial_cond))
        else:
            u = [float(self.initial_cond(i * x_step)) for i in range(x_steps)]
        for j in range(t_start, t_steps):
            if isinstance(self.boundary_cond_1, (float, int)):
                v[0] = float(self.boundary_cond_1)
            else:
                v[0] = self.boundary_cond_1(j * t_step)
            if isinstance(self.boundary_cond_2, (float, int)):
                v[-1] = float(self.boundary_cond_2)
            else:
                v[-1] = self.boundary_cond_2(j * t_step)
            for i in range(x_start, x_steps-1):
                v[i] = u[i] - ((t_step / (2 * x_step)) * (u[i+1] - u[i-1]))
            u = v[:]
        return u

    def mod_euler_method(self, t: int, x: int, t_step: float, x_step: float, t_start: int, x_start: int):
        t_steps = int((t/t_step)+1)
        x_steps = int((x/x_step)+1)

        intermediate, v = np.zeros(x_steps), np.zeros(x_steps)
        if isinstance(self.initial_cond, (float, int)):
            u = np.full(x_steps, float(self.initial_cond))
        else:
            u = [float(self.initial_cond(i * x_step)) for i in range(x_steps)]
        for j in range(t_start, t_steps):
            if isinstance(self.boundary_cond_1, (float, int)):
                intermediate[0], v[0] = float(self.boundary_cond_1), float(self.boundary_cond_1)
            else:
                intermediate[0] = self.boundary_cond_1((j-0.5) * t_step)
                v[0] = self.boundary_cond_1(j * t_step)
            if isinstance(self.boundary_cond_2, (float, int)):
                intermediate[-1], v[-1] = float(self.boundary_cond_2), float(self.boundary_cond_2)
            else:
                intermediate[-1] = self.boundary_cond_2((j-0.5) * t_step)
                v[-1] = self.boundary_cond_2(j * t_step)
            for i in range(x_start, x_steps - 1):
                intermediate[i] = u[i] - (t_step / (4 * x_step)) * (u[i + 1] - u[i - 1])
            for i in range(x_start, x_steps - 1):
                v[i] = u[i] - (t_step / (2 * x_step)) * (intermediate[i + 1] - intermediate[i - 1])
            u = v[:]
        return u


def euler_2_bound_cond(x: float) -> float:
    return math.exp(-(pow((x - 10), 2)))


if __name__ == '__main__':
    time = 5 # int
    displacement = 20  # int
    equation_1 = DE(0, math.sin, 0)
    equation_2 = DE(euler_2_bound_cond, 0, 0)
    x_axis = np.arange(201)

    euler_u_1 = equation_1.euler_method(time, displacement, 0.01, 0.1, 1, 1)
    euler_u_2 = equation_2.euler_method(time, displacement, 0.01, 0.1, 1, 1)

    mod_euler_u_1 = equation_1.mod_euler_method(time, displacement, 0.01, 0.1, 1, 1)
    mod_euler_u_2 = equation_2.mod_euler_method(time, displacement, 0.01, 0.1, 1, 1)

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, figsize=(20, 7))
    line1, = ax1.plot(x_axis, euler_u_1, color='green')
    line2, = ax2.plot(x_axis, euler_u_2, color='red')
    line3, = ax3.plot(x_axis, mod_euler_u_1, color='blue')
    line4, = ax4.plot(x_axis, mod_euler_u_2, color='yellow')
    ax1.set_title("euler_u_1")
    ax2.set_title("euler_u_2")
    ax3.set_title("mod_euler_u_1")
    ax4.set_title("mod_euler_u_2")

    plt.show()

