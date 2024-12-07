import sympy as sym

class Solver():
    def __init__(self):
        self.g, self.m1, self.m2, self.m3, self.m4, self.m_box, self.t, self.r = sym.symbols('g m1 m2 m3 m4 mbox t x y r')

        self.x_box = sym.Function('x_box')(self.t)
        self.y_box = sym.Function('y_box')(self.t)
        self.x_jack = sym.Function('x_jack')(self.t)
        self.y_jack = sym.Function('y_jack')(self.t)
        self.theta_box = sym.Function('theta_box')(self.t)
        self.theta_jack = sym.Function('theta_jack')(self.t)
        
        self.y_matrix = sym.Matrix([0, 1, 0, 0])
        
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
        
    def compute_Lagrangian(self):
        ke_r1 = 1/2 * self.m1 * self.g_w_jack @ (self.g_w_jack @ sym.Matrix(self.r, 0, 0, 1))
        ke_r2 = 1/2 * self.m2 * self.g_w_jack @ (self.g_w_jack @ sym.Matrix(0, -self.r, 0, 1))
        ke_r3 = 1/2 * self.m3 * self.g_w_jack @ (self.g_w_jack @ sym.Matrix(-self.x, 0, 0, 1))
        ke_r4 = 1/2 * self.m4 * self.g_w_jack @ (self.g_w_jack @ sym.Matrix(0, self.y, 0, 1))
        
        pe_r1 = self.m1 * self.g * self.y_matrix @ (self.g_w_jack @ sym.Matrix(self.r, 0, 0, 1))
        pe_r2 = self.m2 * self.g * self.y_matrix @ (self.g_w_jack @ sym.Matrix(0, -self.r, 0, 1))
        pe_r3 = self.m3 * self.g * self.y_matrix @ (self.g_w_jack @ sym.Matrix(-self.x, 0, 0, 1))
        pe_r4 = self.m4 * self.g * self.y_matrix @ (self.g_w_jack @ sym.Matrix(0, self.y, 0, 1))
        
        ke_jack = ke_r1 + ke_r2 + ke_r3 + ke_r4
        pe_jack = pe_r1 + pe_r2 + pe_r3 + pe_r4
        
        l_jack = ke_jack - pe_jack
        
        
        
    def unhat(self):
        pass
    def inv(self):
        pass
    def find_transforms(self):
        pass
    def EL_equations(self):
        pass
    def simulate(self):
        pass
    def integrate(self):
        pass
    def xddot(self):
        pass
    def dyn(self, s):
        pass
    
        
    