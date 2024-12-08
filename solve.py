import sympy as sym
import numpy as np
import matplotlib.pyplot as plt

class Solver():
    def __init__(self):
        self.g, self.m1, self.m2, self.m3, self.m4, self.m_box, self.t, self.r = sym.symbols('g m1 m2 m3 m4 mbox t r')

        self.x_box = sym.Function('x_box')(self.t)
        self.y_box = sym.Function('y_box')(self.t)
        self.x_jack = sym.Function('x_jack')(self.t)
        self.y_jack = sym.Function('y_jack')(self.t)
        self.theta_box = sym.Function('theta_box')(self.t)
        self.theta_jack = sym.Function('theta_jack')(self.t)
        
        self.y_matrix = sym.Matrix([0, 1, 0, 0])
        self.z_matrix = sym.Matrix([0, 0, 0, 1])
        
        self.inertia_matrix_box = sym.Matrix([[self.m_box, 0, 0, 0, 0, 0],
                                         [0, self.m_box, 0, 0, 0, 0],
                                         [0, 0, self.m_box, 0, 0, 0],
                                         [0, 0, 0, 0, 0, 0],
                                         [0, 0, 0, 0, 0, 0],
                                         [0, 0, 0, 0, 0, 1]])
        
        self.q = sym.Matrix([self.x_box, self.y_box, self.theta_box, self.x_jack, self.y_jack, self.theta_jack])
        self.qdot = self.q.diff(self.t)
        self.qddot = self.qdot.diff(self.t)
        
        self.g_w_box_rotation = sym.Matrix([[sym.cos(self.theta_box), -sym.sin(self.theta_box), 0, 0],
                                       [sym.sin(self.theta_box), sym.cos(self.theta_box), 0, 0],
                                       [0, 0, 1, 0],
                                       [0, 0, 0, 1]])
        
        self.g_w_box_translation = sym.Matrix([[1, 0, 0, self.x_box],
                                          [0, 1, 0, self.y_box],
                                          [0, 0, 1, 0],
                                          [0, 0, 0, 1]])
        
        self.g_box_jack_rotation = sym.Matrix([[sym.cos(self.theta_jack), -sym.sin(self.theta_jack), 0, 0],
                                          [sym.sin(self.theta_jack), sym.cos(self.theta_jack), 0, 0],
                                          [0, 0, 1, 0],
                                          [0, 0, 0, 1]])
        
        self.g_box_jack_translation = sym.Matrix([[1, 0, 0, self.x_jack],
                                             [0, 1, 0, self.y_jack],
                                             [0, 0, 1, 0],
                                             [0, 0, 0, 1]])
        
        self.g_w_box = self.g_w_box_rotation @ self.g_box_jack_translation
        self.g_w_jack = self.g_w_box @ (self.g_box_jack_rotation @ self.g_box_jack_translation)
        
        self.subs_dict = {self.g: 9.8, self.m1: 1, self.m2: 1, self.m3: 1, self.m4: 1, self.m_box: 2, self.r: 1}
        self.initial_conditions = [3, 3, np.pi/6, 2.5, 3.2, np.pi/2]
        
    def compute_Lagrangian(self):
        ke_r1 = sym.simplify(1/2 * self.m1 * ((self.g_w_jack @ sym.Matrix([self.r, 0, 0, 1])).diff(self.t)).dot((self.g_w_jack @ sym.Matrix([self.r, 0, 0, 1])).diff(self.t)))
        ke_r2 = sym.simplify(1/2 * self.m2 * ((self.g_w_jack @ sym.Matrix([0, -self.r, 0, 1])).diff(self.t)).dot((self.g_w_jack @ sym.Matrix([0, -self.r, 0, 1])).diff(self.t)))
        ke_r3 = sym.simplify(1/2 * self.m3 * ((self.g_w_jack @ sym.Matrix([-self.r, 0, 0, 1])).diff(self.t)).dot((self.g_w_jack @ sym.Matrix([-self.r, 0, 0, 1])).diff(self.t)))
        ke_r4 = sym.simplify(1/2 * self.m4 * ((self.g_w_jack @ sym.Matrix([0, self.r, 0, 1])).diff(self.t)).dot((self.g_w_jack @ sym.Matrix([0, self.r, 0, 1])).diff(self.t)))
        
        print("ke_r1:", type(ke_r1))
        print("ke_r2:", type(ke_r2))
        print("ke_r3:", type(ke_r3))
        print("ke_r4:", type(ke_r4))
        
        pe_r1 = sym.simplify(self.m1 * self.g * (self.g_w_jack @ sym.Matrix([self.r, 0, 0, 1])).dot(self.y_matrix))
        pe_r2 = sym.simplify(self.m2 * self.g * (self.g_w_jack @ sym.Matrix([0, -self.r, 0, 1])).dot(self.y_matrix))
        pe_r3 = sym.simplify(self.m3 * self.g * (self.g_w_jack @ sym.Matrix([-self.r, 0, 0, 1])).dot(self.y_matrix))
        pe_r4 = sym.simplify(self.m4 * self.g * (self.g_w_jack @ sym.Matrix([0, self.r, 0, 1])).dot(self.y_matrix))
        
        print("pe_r1:", type(pe_r1))
        print("pe_r2:", type(pe_r2))
        print("pe_r3:", type(pe_r3))
        print("pe_r4:", type(pe_r4))
        
        ke_jack = ke_r1 + ke_r2 + ke_r3 + ke_r4
        pe_jack = pe_r1 + pe_r2 + pe_r3 + pe_r4
        
        l_jack = ke_jack - pe_jack
        
        g_w_box_inv = self.inv(self.g_w_box)
        g_w_box_dot = self.g_w_box.diff(self.t)
        
        V_w_box = self.unhat(g_w_box_inv @ g_w_box_dot)
        
        print("V_w_box: ", type(V_w_box))
        
        ke_box = sym.simplify(1/2 * (V_w_box.T @ self.inertia_matrix_box).dot(V_w_box))
        pe_box = sym.simplify(self.m_box * self.g * (self.g_w_box @ self.z_matrix).dot(self.y_matrix))
        
        print("ke box: ", type(ke_box))
        print("pe box: ", type(pe_box))
                
        l_box = ke_box - pe_box
        
        l = l_jack + l_box
        
        return l
        
    def unhat(self, g):
        unhatted_g = sym.Matrix([g[0, -1], g[1, -1], g[2, -1], g[2, 1], g[0, 2], g[1, 0]])
        return unhatted_g
    
    def inv(self, g):
        R_t = sym.Matrix([[g[0, 0], g[1, 0], g[2, 0]],
                    [g[0, 1], g[1, 1], g[2, 1]],
                    [g[0, 2], g[1, 2], g[2, 2]]])

        p = sym.Matrix([g[0, 3], g[1, 3], g[2, 3]])

        R_t_p = -(R_t * p)

        inverse_g = sym.Matrix([[R_t[0, 0], R_t[0, 1], R_t[0, 2], R_t_p[0]],
                                [R_t[1, 0], R_t[1, 1], R_t[1, 2], R_t_p[1]],
                                [R_t[2, 0], R_t[2, 1], R_t[2, 2], R_t_p[2]],
                                [0, 0, 0, 1]])

        return inverse_g
    
    def solve_EL(self):
        l = self.compute_Lagrangian()
        lhs = sym.simplify(sym.Matrix([(l.diff(self.qdot) - l.diff(self.q))]))
        rhs = sym.simplify(sym.Matrix([0, 0, 0, 0, 0, 0]))
        
        el = sym.Eq(lhs, rhs)
        
        solved = sym.solve(el, [self.qddot], dict=False)
        print(solved)
        
        solved_func = sym.Matrix([solved[self.qddot[0]], solved[self.qddot[1]], solved[self.qddot[2]], solved[self.qddot[3]], solved[self.qddot[4]], solved[self.qddot[5]]])
        solved_subs = solved_func.subs(self.subs_dict)
        
        self.lamb_func = sym.lambdify([self.q[0], self.q[1], self.q[2], self.q[3], self.q[4], self.q[5], self.qdot[0], self.qdot[1]. self.qdot[2], self.qdot[3], self.qdot[4], self.qdot[5]], solved_subs)
        
        return solved
        
    def simulate(self, f, x0, tspan, dt, integrate):
        N = int((max(tspan)-min(tspan))/dt)
        x = np.copy(x0)
        tvec = np.linspace(min(tspan),max(tspan),N)
        xtraj = np.zeros((len(x0),N))
        for i in range(1,N):
            xtraj[:,i]=integrate(f,x,dt)
            xtraj[-1, i] = tvec[i]
            x = np.copy(xtraj[:,i])
        return xtraj, tvec
    
    def integrate(self, f, xt, dt):
        k1 = dt * f(xt)
        k2 = dt * f(xt+k1/2.)
        k3 = dt * f(xt+k2/2.)
        k4 = dt * f(xt+k3)
        new_xt = xt + (1/6.) * (k1+2.0*k2+2.0*k3+k4)
        return new_xt
    
    def xddot(self, s):
        return self.lamb_func(s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7], s[8], s[9], s[10], s[11])
    
    def dyn(self, s):
        return np.array([s[6], s[7], s[8], s[9], s[10], s[11], self.xddot(s)[0][0], self.xddot(s)[1][0], self.xddot(s)[2][0], self.xddot(s)[3][0], self.xddot(s)[4][0], self.xddot(s)[5][0]])
    
    def plot(self):
        traj, tvec = self.simulate(self.dyn, self.initial_conditions, [0, 10], 0.01, self.integrate)
        plt.plot(tvec, traj[0])
        plt.plot(tvec, traj[1])
        plt.plot(tvec, traj[2])
        plt.plot(tvec, traj[3])
        plt.plot(tvec, traj[4])
        plt.plot(tvec, traj[5])
        plt.plot(tvec, traj[6])
        plt.plot(tvec, traj[7])
        plt.plot(tvec, traj[8])
        plt.plot(tvec, traj[9])
        plt.plot(tvec, traj[10])
        plt.plot(tvec, traj[11])
        plt.title(r'Theta vs time')
        plt.xlabel("time")
        plt.ylabel(r'Trajectory')
        plt.show()
        
solver = Solver()

solver.solve_EL()