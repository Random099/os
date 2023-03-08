import numpy as np
import matplotlib.pyplot as plt
import math


class DE:
    def __init__(self, initial_cond, boundary_cond_1, boundary_cond_2):
        self.condition_list = {'initial_cond': initial_cond, 'boundary_cond_1': boundary_cond_1,
                               'boundary_cond_2': boundary_cond_2}

    @staticmethod
    def initialize_conditions(condition):
        if isinstance(condition, (float, int)):
            return lambda *args: condition
        else:
            return lambda x: condition(x)

    def euler_method(self, t: int, x: int, t_step: float, x_step: float) -> list[float]:
        t_steps = int((t/t_step)+1)
        x_steps = int((x/x_step)+1)

        l_conditions = {key: DE.initialize_conditions(condition)
                        for key, condition in self.condition_list.items()}
        v = np.zeros(x_steps)
        u = [float(l_conditions['initial_cond'](i * x_step)) for i in range(x_steps)]
        for j in range(1, t_steps):
            v[0] = l_conditions['boundary_cond_1'](j * t_step)
            v[-1] = l_conditions['boundary_cond_2'](j * t_step)
            for i in range(1, x_steps-1):
                v[i] = u[i] - ((t_step / (2 * x_step)) * (u[i+1] - u[i-1]))
            u = v[:]
        return u

    def mod_euler_method(self, t: int, x: int, t_step: float, x_step: float) -> list[float]:
        t_steps = int((t/t_step)+1)
        x_steps = int((x/x_step)+1)

        l_conditions = {key: DE.initialize_conditions(condition)
                        for key, condition in self.condition_list.items()}
        intermediate, v = np.zeros(x_steps), np.zeros(x_steps)
        u = [float(l_conditions['initial_cond'](i * x_step)) for i in range(x_steps)]
        for j in range(1, t_steps):
            intermediate[0] = l_conditions['boundary_cond_1']((j-0.5) * t_step)
            v[0] = l_conditions['boundary_cond_1'](j * t_step)
            intermediate[-1] = l_conditions['boundary_cond_2']((j-0.5) * t_step)
            v[-1] = l_conditions['boundary_cond_2'](j * t_step)
            for i in range(1, x_steps - 1):
                intermediate[i] = u[i] - (t_step / (4 * x_step)) * (u[i + 1] - u[i - 1])
            for i in range(1, x_steps - 1):
                v[i] = u[i] - (t_step / (2 * x_step)) * (intermediate[i + 1] - intermediate[i - 1])
            u = v[:]
        return u


def euler_2_bound_cond(x: float) -> float:
    return math.exp(-(pow((x - 10), 2)))


def plot(data: dict, x_axis=None) -> None:
    lines = [None] * len(data)
    fig, (axs) = plt.subplots(1, len(data), figsize=(15, 7))
    for n, value in enumerate(data.items()):
        lines[n], = axs[n].plot(x_axis, value[1])
        axs[n].locator_params(nbins=5)
        axs[n].set_title(value[0])
    del x_axis
    plt.show()


if __name__ == '__main__':
    TIME, DISPLACEMENT = 6, 20  # int, int
    T_STEP, X_STEP = 0.01, 0.1  # float, float
    T_START, X_START = 1, 1  # int, int
    X_AXIS = np.arange(0, ((DISPLACEMENT/X_STEP)+1)/10, 0.1)
    equation_1 = DE(0, math.sin, 0)
    equation_2 = DE(euler_2_bound_cond, 0, 0)
    TO_PLOT = {'euler_u_1': equation_1.euler_method(TIME, DISPLACEMENT, T_STEP, X_STEP),
    'euler_u_2': equation_2.euler_method(TIME, DISPLACEMENT, T_STEP, X_STEP),
    'mod_euler_u_1': equation_1.mod_euler_method(TIME, DISPLACEMENT, T_STEP, X_STEP),
    'mod_euler_u_2': equation_2.mod_euler_method(TIME, DISPLACEMENT, T_STEP, X_STEP)}
    plot(TO_PLOT, X_AXIS)

