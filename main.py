import numpy as np
import matplotlib.pyplot as plt
import math
import filecmp

class PartialDE():
    def __init__(self, initial_cond, boundary_cond_1, boundary_cond_2):
        self.initial_cond = initial_cond
        self.boundary_cond_1 = boundary_cond_1
        self.boundary_cond_2 = boundary_cond_2

    def euler_method(self, t, x, t_step, x_step, t_start, x_start):
        t_steps = int((t/t_step)+1)
        x_steps = int((x/x_step)+1)

        u = np.zeros(x_steps)
        if isinstance(self.initial_cond, (float, int)):
            v = np.full(x_steps, float(self.initial_cond))
        else:
            v = [self.initial_cond(i * x_step) for i in range(x_steps+1)]
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

    def modified_euler_method(self, t, x, t_step, x_step, t_start, x_start):  # wip
        t_steps = int((t/t_step)+1)
        x_steps = int((x/x_step)+1)

        v = u = np.zeros(x_steps)
        if isinstance(self.initial_cond, (float, int)):
            u = np.full(x_steps, self.initial_cond)
        else:
            u = [self.initial_cond(i * x_step) for i in range(x_steps+1)]
            print(x)
        for j in range(t_start, t_steps):
            if isinstance(self.boundary_cond_1, (float, int)):
                u[0] = self.boundary_cond_1
            else:
                u[0] = self.boundary_cond_1(j * t_step)
            if isinstance(self.boundary_cond_2, (float, int)):
                u[-1] = self.boundary_cond_2
            else:
                u[-1] = self.boundary_cond_2(j * t_step)
            for i in range(x_start, x_steps-1):
                v[i] = u[i] - ((t_step/(2*x_step)) * (u[i+1] - u[i-1]))
            u = v

def mod_euler_bound_cond(x):
    return math.exp(-(pow((x - 10), 2)))

if __name__ == '__main__':
    time = 5
    displacement = 20
    equation_1 = PartialDE(0, math.sin, 0)
    #equation_2 = PartialDE(mod_euler_bound_cond, 0, 0)
    x_axis = np.arange(201)
    u = equation_1.eulers_method(time, displacement, 0.01, 0.1, 1, 1)
    #= equation_2.eulers_method(time, displacement, 0.01, 0.1, 1, 1)

    fig, ax = plt.subplots(figsize=(10, 7))
    line, = ax.plot(x_axis, u, color='green')
    ax.set_title("u")

    plt.show()

