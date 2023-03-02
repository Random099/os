import numpy as np
import matplotlib.pyplot as plt
import math

class PartialDE():
    def __init__(self, initial_cond, boundary_cond_1, boundary_cond_2):
        self.initial_cond = initial_cond
        self.boundary_cond_1 = boundary_cond_1
        self.boundary_cond_2 = boundary_cond_2

    def eulers_method(self, t, x, t_step, x_step, t_start, x_start):
        t_steps = int((t/t_step)+1)
        x_steps = int((x/x_step)+1)
        v = np.zeros(shape=(t_steps, x_steps))
        print(v.shape)

        u_1 = np.zeros(x_steps)
        if isinstance(self.initial_cond, (float, int)):
            v_1 = np.full(x_steps, float(self.initial_cond))
        else:
            v_1 = [self.initial_cond(i * x_step) for i in range(x_steps+1)]
        for j in range(t_start, t_steps):
            if isinstance(self.boundary_cond_1, (float, int)):
                v_1[0] = float(self.boundary_cond_1)
            else:
                v_1[0] = self.boundary_cond_1(j * t_step)
            if isinstance(self.boundary_cond_2, (float, int)):
                v_1[-1] = float(self.boundary_cond_2)
            else:
                v_1[-1] = self.boundary_cond_2(j * t_step)
            for i in range(x_start, x_steps-1):
                v_1[i] = u_1[i] - ((t_step / (2 * x_step)) * (u_1[i+1] - u_1[i-1]))
            for i in range(0, x_steps):
                u_1[i] = v_1[i]

        v_2 = u_2 = np.zeros(x_steps)
        for j in range(t_start, t_steps):
            v_2[0] = math.sin(j * t_step)
            v_2[-1] = 0
            for i in range(x_start, x_steps-1):
                v_2[i] = u_2[i] - (t_step / (2 * x_step) * (u_2[i + 1] - u_2[i - 1]))
            for i in range(0, x_steps):
                u_2[i] = v_2[i]
        np.savetxt('u_1.txt', u_1)
        np.savetxt('u_2.txt', u_2)
        return u_1, u_2

    def modified_eulers_method(self, t, x, t_step, x_step, t_start, x_start):
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
    equation_2 = PartialDE(mod_euler_bound_cond, 0, 0)
    x_axis = np.arange(201)
    u_1, u_2 = equation_1.eulers_method(time, displacement, 0.01, 0.1, 1, 1)
    #= equation_2.eulers_method(time, displacement, 0.01, 0.1, 1, 1)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7))
    line1, = ax1.plot(x_axis, u_1, color='green')
    ax1.set_title("u_1")
    line2, = ax2.plot(x_axis, u_2, color='red')
    ax2.set_title("u_2")
    plt.show()

