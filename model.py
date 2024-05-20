import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


class PsyllidModel:
    def __init__(
        self,
        a=0.57,
        c=0.5,
        lambdu=0.3,
        p=1.22,
        mua=0.76,
        f=7000,
        v=0.75,
        omega=np.pi / 6,
        w=479.72,
        q=0.82,
        mui=0.44,
        g=56,
        gamma=1e-7,
        beta=0.21,
        initial_conditions=None,
        t_span=[0, 350],
    ):

        # Guardamos los params en una lista
        self.params = [a, c, lambdu, p, mua, f, v, omega, w, q, mui, g, gamma, beta]

        # Condiciones iniciales default
        self.default_initial_conditions = [1e7, 1e2, 0, 0, 13000, 0]
        self.initial_conditions = (
            initial_conditions
            if initial_conditions is not None
            else self.default_initial_conditions
        )

        # Time span
        self.t_span = t_span
        self.t_eval = np.linspace(self.t_span[0], self.t_span[1], 1000)

        # Placeholder por si no corren la solución
        self.solution = None

    def model(self, t, y):
        Pu, Pi, Iu, Ii, Tu, Ti = y
        a, c, lambdu, p, mua, f, v, omega, w, q, mui, g, gamma, beta = self.params

        dPu_dt = a * Iu - c * lambdu * Pu * Ti + (1 - c) * p * Pi * Tu - mua * Pu
        dPi_dt = a * Ii + c * lambdu * Pu * Ti - (1 - c) * p * Pi * Tu - mua * Pi
        dIu_dt = (
            w
            * (Pu + Pi)
            * (1 - c)
            * Tu
            * (1 - ((Iu + Ii) / (f * (1 + v * np.sin(omega * t)) * (Ti + Tu))))
            + w
            * (Pu + Pi)
            * c
            * (1 - q)
            * Ti
            * (1 - ((Iu + Ii) / (f * (1 + v * np.sin(omega * t)) * (Ti + Tu))))
            - a * Iu
            - mui * Iu
        )
        dIi_dt = (
            w
            * (Pu + Pi)
            * c
            * q
            * Ti
            * (1 - ((Iu + Ii) / (f * (1 + v * np.sin(omega * t)) * (Ti + Tu))))
            - a * Ii
            - mui * Ii
        )
        dTu_dt = g - gamma * Tu * Pi
        dTi_dt = gamma * Tu * Pi - beta * Ti

        return [dPu_dt, dPi_dt, dIu_dt, dIi_dt, dTu_dt, dTi_dt]

    def solve(self):
        self.solution = solve_ivp(
            self.model,
            self.t_span,
            self.initial_conditions,
            t_eval=self.t_eval,
            method="LSODA",
        )
        return self.solution

    def plot(self):
        if self.solution is None:
            raise ValueError(
                "El modelo no ha sido resuelto aún. Llamá 'model.solve()' primero."
            )

        fig, axs = plt.subplots(1, 2, figsize=(14, 6))

        # Plot P_uninfected and P_infected
        axs[0].plot(
            self.solution.t, self.solution.y[0], label="P_uninfected(t)", linewidth=2.5
        )
        axs[0].plot(
            self.solution.t, self.solution.y[1], label="P_infected(t)", linewidth=2.5
        )
        axs[0].set_xlabel("Meses")
        axs[0].set_ylabel("Población (P)")
        axs[0].set_title("Dinámica poblacional de Psyllid (P)")
        axs[0].legend()

        # Plot T_uninfected and T_infected
        axs[1].plot(
            self.solution.t, self.solution.y[4], label="T_uninfected(t)", linewidth=2.5
        )
        axs[1].plot(
            self.solution.t, self.solution.y[5], label="T_infected(t)", linewidth=2.5
        )
        axs[1].set_xlabel("Meses")
        axs[1].set_ylabel("Población (T)")
        axs[1].set_title("Dinámica poblacional de Árboles (T)")
        axs[1].legend()

        plt.tight_layout()
        plt.show()
