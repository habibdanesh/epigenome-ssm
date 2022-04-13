import numpy as np
import copy
import sys
import time

class ssm(object):
    def __init__(self, seed=0, E=5, G=100, K=3, lambda_1_l2=0.1, lambda_2_l1=0.1, lambda_2_l2=0.1, lambda_3_l2=0.1, positive_state=False, sumone_state=False, positive_em=False, message_passing=True, verbose=False):
        """
        :param E: # assay 
        :param G: # bp
        :param K: # feature
        :param postive_flag: 0:no constraint, 1:abs(Y), 2:max(Y,0)
        :param lambda_1_l2: l2 reg weight of Y
        :param lambda_2_l1: l1 reg weight of Theta
        :param lambda_2_l2: l2 reg weight of Theta
        :param lambda_3_l2: whether to update lambda
        """
        np.random.seed(seed)
        self.E = E
        self.G = G
        self.K = K
        self.lambda_1_l2 = lambda_1_l2
        self.lambda_2_l1 = lambda_2_l1
        self.lambda_2_l2 = lambda_2_l2
        self.lambda_3_l2 = lambda_3_l2
        self.message_passing = message_passing
        self.positive_state = positive_state
        self.sumone_state = sumone_state
        self.positive_em = positive_em
        self.precision = sys.float_info.epsilon
        
        self.x_m = np.asmatrix(np.random.rand(self.E, self.G))
        self.y_m = np.asmatrix(np.random.dirichlet(np.ones(self.K), size=self.G)).T
        self.theta_m = np.asmatrix(np.random.rand(self.K, self.E))
        self.lambda_m = np.asmatrix(np.eye(self.K))
        self.error_m = []
        self.opt_time_m = []
        self.message_dic = {"a_m_f": [], "b_m_f": [], "c_m_f": [], "c_m_b": [], "a_m_b": [], "b_m_b": []}
        
    def initialize(self):
        self.y_m = np.asmatrix(np.random.dirichlet(np.ones(self.K), size=self.G)).T
        g = 0
        while g < self.G:
            self.x_m[:, g] = (
                np.matrix(np.random.multivariate_normal((self.theta_m.T * self.y_m[:, g]).flatten().tolist()[0], np.eye(self.E)))).T
            g += 1
        self.message_dic = {"a_m_f": [], "b_m_f": [], "c_m_f": [], "c_m_b": [], "a_m_b": [], "b_m_b": []}

    def set_x(self, x):
        self.x_m = copy.deepcopy(x)

    def set_y(self, y):
        self.y_m = copy.deepcopy(y)

    def set_theta(self, theta_m):
        self.theta_m = copy.deepcopy(theta_m)

    def set_lambda(self, lambda_m):
        self.lambda_m = copy.deepcopy(lambda_m)

    def total_error(self):
        t_error, e_error, r_error = self.get_error()
        return t_error + e_error + r_error

    def get_error(self):
        g = 0
        t_error = 0 # error comes from transition
        e_error = 0 # error comes from emission
        r_error = 0 # error comes from regularization
        while g < self.G:
            e_error += np.sum((self.x_m[:, g] - self.theta_m.T * self.y_m[:, g]).T \
                     * (self.x_m[:, g] - self.theta_m.T * self.y_m[:, g]))
            r_error += np.sum(self.lambda_1_l2 * self.y_m[:, g].T * self.y_m[:, g]) # penalization
          
            if g < self.G - 1:
                t_error += np.sum((self.y_m[:, g + 1] - self.lambda_m * self.y_m[:, g]).T * (
                        self.y_m[:, g + 1] - self.lambda_m * self.y_m[:, g]))
               
            g += 1
        r_error += np.sum(self.lambda_2_l2*np.multiply(self.theta_m, self.theta_m))
        r_error += np.sum(self.lambda_3_l2*np.multiply(self.lambda_m, self.lambda_m))
        return t_error, e_error, r_error

    def print_error(self):
        t_error, e_error, r_error = self.get_error()
        print("\n======Error=====")
        print("Error total: {}".format(t_error+e_error+r_error))
        print("Em: {}".format(e_error))
        print("Trans: {}".format(t_error))
        print("Reg: {}".format(r_error))
        print("=================\n")

    def get_theta_msg(self):
        lhs_m = np.asmatrix(np.zeros((self.K, self.K)))
        rhs_m = np.asmatrix(np.zeros((self.K, self.E)))
        g = 0
        while g < self.G:
            lhs_m += self.y_m[:, g] * self.y_m[:, g].T
            rhs_m += self.y_m[:, g] * self.x_m[:, g].T
            g += 1
        lhs_m += self.lambda_2_l2 * np.eye(self.K)
        lhs_m = np.linalg.pinv(lhs_m)
        # rhs_m + 1/2 lambda_2_l1 * J'
        J = np.sign(self.theta_m)
        J[np.where(J==0)] = 1
        rhs_m -= J * self.lambda_2_l1/2
        return lhs_m, rhs_m
    
    def get_lambda_msg(self):
        lhs_m = np.asmatrix(np.zeros((self.K, self.K)))
        rhs_m = np.asmatrix(np.zeros((self.K, self.K)))
        g = 0
        while g < self.G - 1:
            lhs_m += self.y_m[:,g] * self.y_m[:,g].T
            rhs_m += self.y_m[:,g] * self.y_m[:,g+1].T
            g += 1
        return lhs_m, rhs_m

    def forward(self):
        g = 0
        # initialize
        a_m = self.theta_m * self.theta_m.T + self.lambda_1_l2 * np.eye(self.K)
        b_m = (-2 * self.x_m[:, g].T * self.theta_m.T).T
        # c_m = self.x_m[:, g].T * self.x_m[:, g]
        self.message_dic["a_m_f"].append(copy.deepcopy(a_m))
        self.message_dic["b_m_f"].append(copy.deepcopy(b_m))
        # self.message_dic["c_m_f"].append(copy.deepcopy(c_m))
        while g <= self.G - 2:
            # update
            g += 1
            g_m = (self.lambda_m.T * self.lambda_m + a_m.T).I
            a_m_next = self.lambda_1_l2 * np.eye(self.K) \
                       + self.theta_m * self.theta_m.T \
                       + (np.eye(self.K) - self.lambda_m * g_m * self.lambda_m.T).T * (
                               np.eye(self.K) - self.lambda_m * g_m * self.lambda_m.T) \
                       + self.lambda_m * g_m.T * a_m * g_m * self.lambda_m.T
            b_m_next = (b_m.T * g_m.T * self.lambda_m.T * (np.eye(self.K) - self.lambda_m * g_m * self.lambda_m.T) \
                        - b_m.T * g_m.T * a_m * g_m * self.lambda_m.T \
                        + b_m.T * g_m * self.lambda_m.T \
                        - 2 * self.x_m[:, g].T * self.theta_m.T).T
            # c_m_next = self.x_m[:, g].T * self.x_m[:, g] \
            #            + 1 / 4 * b_m.T * g_m.T * self.lambda_m.T * self.lambda_m * g_m * b_m \
            #            + 1 / 4 * b_m.T * g_m.T * a_m * g_m * b_m \
            #            - 1 / 2 * b_m.T * g_m * b_m \
            #            + c_m
            a_m = a_m_next
            b_m = b_m_next
            # c_m = c_m_next
            self.message_dic["a_m_f"].append(copy.deepcopy(a_m))
            self.message_dic["b_m_f"].append(copy.deepcopy(b_m))
            # self.message_dic["c_m_f"].append(copy.deepcopy(c_m))

    def backward(self):
        g = self.G - 1
        # initialize
        a_m = self.theta_m * self.theta_m.T + self.lambda_1_l2 * np.eye(self.K)
        b_m = (-2 * self.x_m[:, g].T * self.theta_m.T).T
        # c_m = self.x_m[:, g].T * self.x_m[:, g]
        self.message_dic["a_m_b"].insert(0, copy.deepcopy(a_m))
        self.message_dic["b_m_b"].insert(0, copy.deepcopy(b_m))
        # self.message_dic["c_m_b"].insert(0, copy.deepcopy(c_m))
        while g > 0:
            g -= 1
            g_m = (np.eye(self.K) + a_m.T).I
            a_m_next = self.lambda_1_l2 * np.eye(self.K) \
                       + self.theta_m * self.theta_m.T \
                       + (g_m * self.lambda_m - self.lambda_m).T * (g_m * self.lambda_m - self.lambda_m) \
                       + self.lambda_m.T * g_m.T * a_m * g_m * self.lambda_m
            b_m_next = (- 2 * self.x_m[:, g].T * self.theta_m.T \
                        + b_m.T * g_m * self.lambda_m \
                        - b_m.T * g_m.T * a_m * g_m * self.lambda_m \
                        - b_m.T * g_m.T * (g_m * self.lambda_m - self.lambda_m)).T
            # c_m_next = self.x_m[:, g].T * self.x_m[:, g] \
            #            + 1 / 4 * b_m.T * g_m.T * g_m * b_m \
            #            + 1 / 4 * b_m.T * g_m.T * a_m * g_m * b_m \
            #            - 1 / 2 * b_m.T * g_m * b_m \
            #            + c_m
            a_m = a_m_next
            b_m = b_m_next
            # c_m = c_m_next
            self.message_dic["a_m_b"].insert(0, copy.deepcopy(a_m))
            self.message_dic["b_m_b"].insert(0, copy.deepcopy(b_m))
            # self.message_dic["c_m_b"].insert(0, copy.deepcopy(c_m))

    def update_y(self):
        g = 0
        while g < self.G:
            a_m_f = self.message_dic["a_m_f"][g]
            b_m_f = self.message_dic["b_m_f"][g]
            a_m_b = self.message_dic["a_m_b"][g]
            b_m_b = self.message_dic["b_m_b"][g]
            # c_m_f = self.message_dic["c_m_f"][g]
            # c_m_b = self.message_dic["c_m_b"][g]
            if self.sumone_state:
                D = np.asmatrix(np.zeros((1, self.K)))
                d = np.asmatrix(np.zeros((1, 1)))
                # sum one constraint
                for k in range(self.K):
                    D[0, k] = 1
                d[0, 0] = 1
                # constraint D.T D Y_g = D.T d
                lhs_m = (-a_m_f - a_m_b + self.theta_m * self.theta_m.T + self.lambda_1_l2*np.eye(self.K))
                rhs_m = (1 / 2 * b_m_f + 1 / 2 * b_m_b + self.theta_m * self.x_m[:, g])
#                 print("lhs_m.I:\n", lhs_m.I)
#                 print("rhs_m:\n", rhs_m)
#                 print("D:\n", D)
#                 print("D*lhs_m.I*D.T):\b", (D*lhs_m.I*D.T))
#                 print("D*lhs_m.I):\b", (D*lhs_m.I))
                lagrange_term = (D*lhs_m.I*D.T).I*(-d + D*lhs_m.I*rhs_m)
                new_y_m = lhs_m.I * (-D.T*lagrange_term + rhs_m)
#                 print("[state]D*Y_g:", D*self.y_m[:, g])
#                 print("[state]before lag:", (-a_m_f - a_m_b + self.theta_m * self.theta_m.T + self.lambda_1_l2*np.eye(self.K)).I * (
#                     1 / 2 * b_m_f + 1 / 2 * b_m_b + self.theta_m * self.x_m[:, g]))
#                 print("[state]after lag:", self.y_m[:, g])
            else:
                new_y_m = (-a_m_f - a_m_b + self.theta_m * self.theta_m.T + self.lambda_1_l2*np.eye(self.K)).I * (
                    1 / 2 * b_m_f + 1 / 2 * b_m_b + self.theta_m * self.x_m[:, g])
            if self.positive_state:
                """
                r = new_y_m - self.y_m[:, g]
                min_scale = 1
                min_k = 0
                for k in range(self.K):
                    if self.y_m[k,g] < 0:
                        if r[k,0] < 0:
                            min_scale = 0
                            break
                    else:
                        if new_y_m[k,0] < 0:
                            scale = (0.0 - self.y_m[k,g]) / r[k,0]
                            if scale < min_scale and scale >= 0:
                                min_scale = scale
                                min_k = k
                self.y_m[:, g] = self.y_m[:, g] + min_scale * r
                # some value can be left as negative, enforce them to be 0.0
                for k in range(self.K):
                    if self.y_m[k, g] < 0:
                        self.y_m[k, g] = 0
                """
                active_set = set()
                for _ in range(self.K):    
                    if self.check_positive(new_y_m):
                        break
                    D = np.asmatrix(np.zeros((self.K, self.K)))
                    d = np.asmatrix(np.zeros((self.K, 1)))
                    for k in range(self.K):
                        if new_y_m[k,0] <= self.precision:
                            active_set.add(k)
                    for k in active_set:        
                        D[k, k] = 1
                        d[k, 0] = 0
                    # constraint D.T D Y_g = D.T d
                    lhs_m = (-a_m_f - a_m_b + self.theta_m * self.theta_m.T + self.lambda_1_l2*np.eye(self.K))
                    rhs_m = (1 / 2 * b_m_f + 1 / 2 * b_m_b + self.theta_m * self.x_m[:, g])
    #                 print("lhs_m.I:\n", lhs_m.I)
    #                 print("rhs_m:\n", rhs_m)
    #                 print("D:\n", D)
    #                 print("D*lhs_m.I*D.T):\b", (D*lhs_m.I*D.T))
    #                 print("D*lhs_m.I):\b", (D*lhs_m.I))
                    lagrange_term = np.linalg.pinv((D*lhs_m.I*D.T))*(-d + D*lhs_m.I*rhs_m)
                    new_y_m = lhs_m.I * (-D.T*lagrange_term + rhs_m)
    #                 print("[state]D*Y_g:", D*self.y_m[:, g])
    #                 print("[state]before lag:", (-a_m_f - a_m_b + self.theta_m * self.theta_m.T + self.lambda_1_l2*np.eye(self.K)).I * (
    #                     1 / 2 * b_m_f + 1 / 2 * b_m_b + self.theta_m * self.x_m[:, g]))
    #                 print("[state]after lag:", self.y_m[:, g])
                if not self.check_positive(new_y_m):
                    # print("[state] negative state!")
                    # print(new_y_m)
                    pass
                self.y_m[:, g] = new_y_m
            else:
                self.y_m[: ,g] = new_y_m
            g += 1

    def update_y_idv(self):
        g = 0
        while g < self.G:
            if g == 0:
                self.y_m[:, g] = (self.theta_m * self.theta_m.T + self.lambda_m.T * self.lambda_m + (1 + self.lambda_1_l2) * np.eye(self.K)).I \
                                 * (self.lambda_m.T * self.y_m[:, g + 1] + self.theta_m * self.x_m[:, g])
            elif g == self.G - 1:
                self.y_m[:, g] = (self.theta_m * self.theta_m.T + self.lambda_m.T * self.lambda_m + (1 + self.lambda_1_l2) * np.eye(self.K)).I \
                                 * (self.lambda_m * self.y_m[:, g - 1] + self.theta_m * self.x_m[:, g])
            else:
                self.y_m[:, g] = (self.theta_m * self.theta_m.T + self.lambda_m.T * self.lambda_m + (1 + self.lambda_1_l2)  * np.eye(self.K)).I \
                                 * (self.lambda_m * self.y_m[:, g - 1] + self.lambda_m.T * self.y_m[:,g + 1] + self.theta_m * self.x_m[:, g])
            g += 1

    def check_positive(self, a):
        m,n = a.shape
        for i in range(m):
            for j in range(n):
                if a[i,j] < 0:
                    return False
        return True
    
    def update_theta(self, lhs_m, rhs_m):
        """
        this is a copy from ssmKalman, which is using different symbol
        here new_z_m(E by K) is the transpose of theta(K by E)
        """
        t = 0
        # FIXIME: fix the symbol
        self.P = self.E
        self.M = self.K
        self.T = self.G
        self.z_m = self.theta_m.T.copy()
        new_z_m = (lhs_m * rhs_m).T
        z_m_noLag = new_z_m.copy()
        if self.positive_em:
            # print('rhs_m', rhs_m)
            new_z_m = np.asmatrix(np.zeros((self.M, self.P))).T
            for p in range(self.P):
                old_z = self.z_m[p, :].T
                lag_lambda_m = np.asmatrix(np.zeros((1, self.M)))
                # get active set
                active = set()
                for _ in range(self.M):
                    # print('[z]:before update\n', old_z)
                    # print('[z]:update without lag\n', z_m_noLag[p, :])
                    # print('[z]:lag term\n', lag_lambda_m[p,:].T)
                    new_z = lhs_m * (rhs_m[:,p] + 1/2*lag_lambda_m[0,:].T)
                    # print('[z]:after lag\n', new_z)

                    # print('after round up\n', new_z)
                    v = new_z - old_z
                    min_k = 1
                    min_m = 0
                    found = False
                    for m in range(self.M):
                        if new_z[m,0] < 0:
                            if v[m, 0] < 0:
                                k = old_z[m, 0] / (-v[m, 0])
                                if k < min_k:
                                    min_k = k
                                    min_m = m
                                    found = True
                    old_z += min_k*v
                    for m in range(self.M):
                        if -self.precision <= old_z[m,0] <= self.precision:
                            old_z[m,0] = 0
                    # print('[z]:updated z', old_z[:, 0])
                    
                    if found:
                        active.add(min_m)
                    # print(active)
                    # print('[z]:active set', active)
                    # set lag term to 0
                    done = False
                    while not done:
                        lag_lambda_m = np.asmatrix(np.zeros((1, self.M)))
                        
                        if len(active) > 0:
                            active_len = len(active)
                            lag_lambda = np.asmatrix(np.zeros((1, active_len)))
                            lag_lhs = np.asmatrix(np.zeros((active_len, active_len)))
                            lag_lhs_bar = np.asmatrix(np.zeros((active_len, self.M - active_len)))
                            lag_rhs = np.asmatrix(np.zeros((active_len, 1)))
                            lag_rhs_bar = np.asmatrix(np.zeros((self.M - active_len, 1)))
                            for r, val in enumerate(active):
                                counter = 0
                                bar_counter = 0
                                for c in range(self.M):
                                    if c in active:
                                        lag_lhs[r, counter] = lhs_m[val, c]
                                        counter += 1
                                    else:
                                        lag_lhs_bar[r, bar_counter] = lhs_m[val, c]
                                        bar_counter += 1
                            counter = 0
                            bar_counter = 0
                            for c in range(self.M):
                                if c in active:
                                    lag_rhs[counter, 0] = rhs_m[c, p]
                                    counter += 1
                                else:
                                    lag_rhs_bar[bar_counter, 0] = rhs_m[c, p]
                                    bar_counter += 1
                            lag_lambda = (-2 * np.linalg.pinv(lag_lhs) * lag_lhs_bar * lag_rhs_bar - 2 * lag_rhs).T
                            # print('lag_lambda\n', lag_lambda)
                            # print('lag_lhs\n', lag_lhs)
                            # print('lag_lhs_bar\n', lag_lhs_bar)
                            # print('lag_rhs\n', lag_rhs)
                            # print('lag_rhs_bar\n', lag_rhs_bar)
                            expected_theta = lag_lhs * (lag_rhs + 1/2*lag_lambda.T) + lag_lhs_bar * lag_rhs_bar
                            # print('expected theta\n', expected_theta) 
                            
                            active_remove = set()
                            # print("[Neg lag]", lag_lambda)
                            for c, val in enumerate(active):
                                if lag_lambda[:, c] < 0:
                                    print("[Neg lag]", lag_lambda[:, c])
                                    active_remove.add(val)
                            active = active.difference(active_remove)
                            if len(active_remove) == 0:
                                done = True
                                for c, val in enumerate(active):
                                    lag_lambda_m[:, val] = lag_lambda[:, c]
                        else:
                            break
                if not self.check_positive(old_z):
                    print("[Warning]: get negative z value")
                    print(old_z)
                    # raise Exception('[z]:Non positive state')
                new_z_m[p, :] = old_z[:, 0].T
        self.theta_m = new_z_m.T.copy()   

    def update_lambda(self, lhs_m, rhs_m):
        self.lambda_m = np.asmatrix(np.linalg.lstsq(lhs_m + self.lambda_3_l2 * np.eye(self.K), rhs_m, rcond=1)[0]).T

    def update_state(self):
        # initialzie the self.message_dic for a_m,b_m
        self.message_dic["a_m_f"] = []
        self.message_dic["b_m_f"] = []
        self.message_dic["a_m_b"] = []
        self.message_dic["b_m_b"] = []
        self.message_dic["c_m_f"] = []
        self.message_dic["c_m_b"] = []
        # use forward_backward message to update y
        self.forward()
        self.backward()
        if self.message_passing:
            self.update_y()
        else:
            self.update_y_idv()

    def optimization(self, iteration=10, verbose=False):
        self.message_dic = {"a_m_f": [], "b_m_f": [], "c_m_f": [], "c_m_b": [], "a_m_b": [], "b_m_b": []}
        i = 0
        self.error_m.append(self.total_error())
        if verbose:
            print("[Error]: intial error", self.total_error())
        start_time = time.time()
        while i < iteration:
            self.update_state()
            if verbose:
                print("[Error] after update state:", self.total_error())
            
            lambda_lhs, lambda_rhs = self.get_lambda_msg()
            self.update_lambda(lambda_lhs, lambda_rhs)
            if verbose:
                print("[Error] after update tr:", self.total_error())

            
            theta_lhs, theta_rhs = self.get_theta_msg()
            self.update_theta(theta_lhs, theta_rhs)
            if verbose:
                print("[Error] after update em:", self.total_error())

            self.error_m.append(self.total_error())
            self.opt_time_m.append(time.time() - start_time)
            i += 1
        self.error_begin = self.error_m[0]
        self.error_end = self.error_m[-1]

if __name__ == "__main__":
    seed = 0
    M = 3
    T = 100
    P = 5
    lambda_3_l2 = 0.001
    lambda_2_l2 = 0.001
    lambda_2_l1 = 100
    lambda_1_l2 = 0.001
    # experiment
    model = ssm(seed=seed, E=P, G=T, K=M, lambda_1_l2=lambda_1_l2, 
        lambda_2_l1=lambda_2_l1, lambda_2_l2=lambda_2_l2, lambda_3_l2=lambda_3_l2, message_passing=True, sumone_state=False, positive_state=True, positive_em=True)

    test_iteration = 10
    model.optimization(test_iteration, verbose=True)
    model.print_error()
    print("errSSM:\n", model.error_m)
    print("emSSM:\n", model.theta_m)
    print("trSSM:\n", model.lambda_m)
    # for g in range(model.G):
    #     print("pos{}:\n".format(g), model.y_m[:,g].T)
    #     print("sum{}:".format(np.sum(model.y_m[:,g])))
        
    # plot
    # import time
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt


    # plot, predict vs real trans value
    fig, axarr = plt.subplots()
    axarr.plot(model.error_m, label='normal')
    axarr.legend(loc='upper right')
    plt.xlabel('iteration')
    plt.ylabel('error')
    plt.show()
