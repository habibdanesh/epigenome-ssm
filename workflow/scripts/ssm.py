import os
import sys
import copy
import time
import numpy as np

class ssm(object):
    def __init__(self, seed=0, 
                E=5, G=100, K=3, 
                lambda_1_l2=0.1, lambda_2_l1=0., lambda_2_l2=0.1, lambda_3_l2=0.1, 
                positive_state=True, sumone_state=False, positive_em=True, message_passing=True, 
                n_threads=1, verbose=False):
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
        self.seed = seed
        self.E = E
        self.G = G
        self.K = K
        self.lambda_1_l2 = lambda_1_l2
        self.lambda_2_l1 = lambda_2_l1
        self.lambda_2_l2 = lambda_2_l2
        self.lambda_3_l2 = lambda_3_l2
        self.positive_state = positive_state
        self.sumone_state = sumone_state
        self.positive_em = positive_em
        self.message_passing = message_passing
        self.n_threads = n_threads
        self.verbose = verbose
        self.precision = sys.float_info.epsilon
        
        self.x_m = np.asmatrix(np.random.rand(self.E, self.G))
        self.y_m = np.asmatrix(np.random.dirichlet(np.ones(self.K), size=self.G)).T
        self.theta_m = np.asmatrix(np.random.rand(self.K, self.E))
        self.theta_m_transpose = self.theta_m.T
        self.lambda_m = np.asmatrix(np.eye(self.K))
        self.lambda_m_transpose = self.lambda_m.T
        self.error_m = []
        self.opt_time_m = []
        self.message_dic = {"a_m_f": [], "b_m_f": [], "c_m_f": [], "c_m_b": [], "a_m_b": [], "b_m_b": []}
        
    def initialize(self):
        self.y_m = np.asmatrix(np.random.dirichlet(np.ones(self.K), size=self.G)).T
        g = 0
        while g < self.G:
            self.x_m[:, g] = (
                np.matrix(np.random.multivariate_normal((self.theta_m_transpose * self.y_m[:, g]).flatten().tolist()[0], np.eye(self.E)))).T
            g += 1
        self.message_dic = {"a_m_f": [], "b_m_f": [], "c_m_f": [], "c_m_b": [], "a_m_b": [], "b_m_b": []}

    def re_init(self, seed):
        np.random.seed(seed)
        self.seed = seed
        self.y_m = np.asmatrix(np.random.dirichlet(np.ones(self.K), size=self.G)).T
        self.theta_m = np.asmatrix(np.random.rand(self.K, self.E))
        self.theta_m_transpose = self.theta_m.T
        self.lambda_m = np.asmatrix(np.eye(self.K))
        self.lambda_m_transpose = self.lambda_m.T

    def set_x(self, x):
        self.x_m = copy.deepcopy(x)

    def set_y(self, y):
        self.y_m = copy.deepcopy(y)

    def set_theta(self, theta_m):
        self.theta_m = copy.deepcopy(theta_m)
        self.theta_m_transpose = self.theta_m.T

    def set_lambda(self, lambda_m):
        self.lambda_m = copy.deepcopy(lambda_m)
        self.lambda_m_transpose = self.lambda_m.T

    def total_error(self):
        t_error, e_error, r_error = self.get_error()
        return t_error + e_error + r_error

    def get_error(self):
        g = 0
        t_error = 0 # error comes from transition
        e_error = 0 # error comes from emission
        r_error = 0 # error comes from regularization
        while g < self.G:
            e_error += np.sum((self.x_m[:, g] - self.theta_m_transpose * self.y_m[:, g]).T \
                     * (self.x_m[:, g] - self.theta_m_transpose * self.y_m[:, g]))
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
        J = np.sign(self.theta_m)
        J[np.where(J==0)] = 1
        rhs_m -= J * self.lambda_2_l1 / 2
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
        K_eye = np.eye(self.K)
        a_m_initial = self.theta_m * self.theta_m_transpose + self.lambda_1_l2 * K_eye
        a_m = a_m_initial
        b_m = (-2 * self.x_m[:, g].T * self.theta_m_transpose).T
        self.message_dic["a_m_f"].append(copy.deepcopy(a_m))
        self.message_dic["b_m_f"].append(copy.deepcopy(b_m))
        lambda_m_square = self.lambda_m_transpose * self.lambda_m
        while g <= self.G - 2:
            # update
            g += 1
            b_m_transpose = b_m.T
            g_m = (lambda_m_square + a_m.T).I
            g_m_transpose = g_m.T
            expression_1 = g_m * self.lambda_m_transpose
            expression_2 = K_eye - self.lambda_m * expression_1
            expression_3 = g_m_transpose * a_m * expression_1
            a_m_next = a_m_initial + (expression_2).T * (expression_2) \
                       + self.lambda_m * expression_3
            b_m_next = (b_m_transpose * g_m_transpose * self.lambda_m_transpose * (expression_2) \
                        - b_m_transpose * expression_3 \
                        + b_m_transpose * expression_1 \
                        - 2 * self.x_m[:, g].T * self.theta_m_transpose).T
            a_m = a_m_next
            b_m = b_m_next
            self.message_dic["a_m_f"].append(copy.deepcopy(a_m))
            self.message_dic["b_m_f"].append(copy.deepcopy(b_m))

    def backward(self):
        g = self.G - 1
        # initialize
        K_eye = np.eye(self.K)
        a_m_initial = self.theta_m * self.theta_m_transpose + self.lambda_1_l2 * K_eye
        a_m = a_m_initial
        b_m = (-2 * self.x_m[:, g].T * self.theta_m_transpose).T
        self.message_dic["a_m_b"].insert(0, copy.deepcopy(a_m))
        self.message_dic["b_m_b"].insert(0, copy.deepcopy(b_m))
        while g > 0:
            g -= 1
            b_m_transpose = b_m.T
            g_m = (K_eye + a_m.T).I
            g_m_transpose = g_m.T
            expression_1 = g_m * self.lambda_m
            expression_2 = expression_1 - self.lambda_m
            expression_3 = a_m * expression_1
            expression_4 = b_m_transpose * g_m_transpose
            a_m_next = a_m_initial + (expression_2).T * (expression_2) \
                       + self.lambda_m_transpose * g_m_transpose * expression_3
            b_m_next = (- 2 * self.x_m[:, g].T * self.theta_m_transpose \
                        + b_m_transpose * expression_1 \
                        - expression_4 * expression_3 \
                        - expression_4 * (expression_2)).T
            a_m = a_m_next
            b_m = b_m_next
            self.message_dic["a_m_b"].insert(0, copy.deepcopy(a_m))
            self.message_dic["b_m_b"].insert(0, copy.deepcopy(b_m))
    
    def update_y(self):
        expression_1 = self.theta_m * self.theta_m_transpose + self.lambda_1_l2 * np.eye(self.K)
        g = 0
        while g < self.G:
            a_m_f = self.message_dic["a_m_f"][g]
            b_m_f = self.message_dic["b_m_f"][g]
            a_m_b = self.message_dic["a_m_b"][g]
            b_m_b = self.message_dic["b_m_b"][g]
            lhs_m = (-a_m_f - a_m_b + expression_1)
            lhs_m_inverse = lhs_m.I
            rhs_m = (1 / 2 * b_m_f + 1 / 2 * b_m_b + self.theta_m * self.x_m[:, g])
            if self.sumone_state:
                D = np.asmatrix(np.zeros((1, self.K)))
                D_transpose = D.T
                d = np.asmatrix(np.zeros((1, 1)))
                # sum one constraint
                for k in range(self.K):
                    D[0, k] = 1
                d[0, 0] = 1
                expression_2 = D * lhs_m_inverse
                lagrange_term = (expression_2 * D_transpose).I * (-d + expression_2 * rhs_m)
                new_y_m = lhs_m_inverse * (-D_transpose * lagrange_term + rhs_m)
            else:
                new_y_m = lhs_m_inverse * rhs_m
            if self.positive_state:
                active_set = set()
                for _ in range(self.K):    
                    if not np.any(new_y_m < 0): # if all positive
                        break
                    D = np.asmatrix(np.zeros((self.K, self.K)))
                    D_transpose = D.T
                    d = np.asmatrix(np.zeros((self.K, 1)))
                    for k in range(self.K):
                        if new_y_m[k, 0] <= self.precision:
                            active_set.add(k)
                    for k in active_set:        
                        D[k, k] = 1
                        d[k, 0] = 0
                    expression_2 = D * lhs_m_inverse
                    lagrange_term = np.linalg.pinv(expression_2 * D_transpose) * (-d + expression_2 * rhs_m)
                    new_y_m = lhs_m_inverse * (-D_transpose * lagrange_term + rhs_m)
                if np.any(new_y_m < 0): # if any negative
                    pass
                self.y_m[:, g] = new_y_m
            else:
                self.y_m[: ,g] = new_y_m
            g += 1

    def update_y_idv(self):
        g = 0
        while g < self.G:
            if g == 0:
                self.y_m[:, g] = (self.theta_m * self.theta_m_transpose + self.lambda_m_transpose * self.lambda_m + (1 + self.lambda_1_l2) * np.eye(self.K)).I \
                                 * (self.lambda_m_transpose * self.y_m[:, g + 1] + self.theta_m * self.x_m[:, g])
            elif g == self.G - 1:
                self.y_m[:, g] = (self.theta_m * self.theta_m_transpose + self.lambda_m_transpose * self.lambda_m + (1 + self.lambda_1_l2) * np.eye(self.K)).I \
                                 * (self.lambda_m * self.y_m[:, g - 1] + self.theta_m * self.x_m[:, g])
            else:
                self.y_m[:, g] = (self.theta_m * self.theta_m_transpose + self.lambda_m_transpose * self.lambda_m + (1 + self.lambda_1_l2)  * np.eye(self.K)).I \
                                 * (self.lambda_m * self.y_m[:, g - 1] + self.lambda_m_transpose * self.y_m[:,g + 1] + self.theta_m * self.x_m[:, g])
            g += 1
    
    def update_theta(self, lhs_m, rhs_m):
        """
        this is a copy from ssmKalman, which is using different symbol
        here new_z_m(E by K) is the transpose of theta(K by E)
        """
        # FIXIME: fix the symbol
        self.z_m = self.theta_m_transpose.copy()
        new_z_m = (lhs_m * rhs_m).T
        if self.positive_em:
            new_z_m = np.asmatrix(np.zeros((self.K, self.E))).T
            for p in range(self.E):
                old_z = self.z_m[p, :].T
                lag_lambda_m = np.asmatrix(np.zeros((1, self.K)))
                # get active set
                active = set()
                for _ in range(self.K):
                    new_z = lhs_m * (rhs_m[:,p] + 1/2*lag_lambda_m[0,:].T)
                    v = new_z - old_z
                    min_k = 1
                    min_m = 0
                    found = False
                    for m in range(self.K):
                        if new_z[m,0] < 0:
                            if v[m, 0] < 0:
                                k = old_z[m, 0] / (-v[m, 0])
                                if k < min_k:
                                    min_k = k
                                    min_m = m
                                    found = True
                    old_z += min_k*v
                    for m in range(self.K):
                        if -self.precision <= old_z[m,0] <= self.precision:
                            old_z[m,0] = 0                    
                    if found:
                        active.add(min_m)
                    # set lag term to 0
                    done = False
                    while not done:
                        lag_lambda_m = np.asmatrix(np.zeros((1, self.K)))
                        if len(active) > 0:
                            active_len = len(active)
                            lag_lambda = np.asmatrix(np.zeros((1, active_len)))
                            lag_lhs = np.asmatrix(np.zeros((active_len, active_len)))
                            lag_lhs_bar = np.asmatrix(np.zeros((active_len, self.K - active_len)))
                            lag_rhs = np.asmatrix(np.zeros((active_len, 1)))
                            lag_rhs_bar = np.asmatrix(np.zeros((self.K - active_len, 1)))
                            for r, val in enumerate(active):
                                counter = 0
                                bar_counter = 0
                                for c in range(self.K):
                                    if c in active:
                                        lag_lhs[r, counter] = lhs_m[val, c]
                                        counter += 1
                                    else:
                                        lag_lhs_bar[r, bar_counter] = lhs_m[val, c]
                                        bar_counter += 1
                            counter = 0
                            bar_counter = 0
                            for c in range(self.K):
                                if c in active:
                                    lag_rhs[counter, 0] = rhs_m[c, p]
                                    counter += 1
                                else:
                                    lag_rhs_bar[bar_counter, 0] = rhs_m[c, p]
                                    bar_counter += 1
                            lag_lambda = (-2 * np.linalg.pinv(lag_lhs) * lag_lhs_bar * lag_rhs_bar - 2 * lag_rhs).T
                            active_remove = set()
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
                if np.any(old_z < 0):
                    print("[Warning]: get negative z value")
                    print(old_z)
                new_z_m[p, :] = old_z[:, 0].T
        self.theta_m = new_z_m.T.copy()
        self.theta_m_transpose = self.theta_m.T

    def update_lambda(self, lhs_m, rhs_m):
        self.lambda_m = np.asmatrix(np.linalg.lstsq(lhs_m + self.lambda_3_l2 * np.eye(self.K), rhs_m, rcond=1)[0]).T
        self.lambda_m_transpose = self.lambda_m.T

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

    def optimization(self, iteration=10):
        self.message_dic = {"a_m_f": [], "b_m_f": [], "c_m_f": [], "c_m_b": [], "a_m_b": [], "b_m_b": []}
        i = 0
        self.error_m.append(self.total_error())
        if self.verbose:
            print("[Error]: intial error", self.total_error())
        start_time = time.time()
        while i < iteration:
            self.update_state()
            if self.verbose:
                print("[Error] after update state:", self.total_error())
            
            lambda_lhs, lambda_rhs = self.get_lambda_msg()
            self.update_lambda(lambda_lhs, lambda_rhs)
            if self.verbose:
                print("[Error] after update tr:", self.total_error())

            
            theta_lhs, theta_rhs = self.get_theta_msg()
            self.update_theta(theta_lhs, theta_rhs)
            if self.verbose:
                print("[Error] after update em:", self.total_error())

            self.error_m.append(self.total_error())
            self.opt_time_m.append(time.time() - start_time)
            i += 1
        self.error_begin = self.error_m[0]
        self.error_end = self.error_m[-1]

if __name__ == "__main__":
    seed = 0
    E = 30
    G = 10000
    K = 5
    LAMBDA_1_L2 = .1
    LAMBDA_2_L1 = .0
    LAMBDA_2_L2 = .1
    LAMBDA_3_L2 = .1
    ### run ssm
    model = ssm(seed=seed, E=E, G=G, K=K, 
                lambda_1_l2=LAMBDA_1_L2, lambda_2_l1=LAMBDA_2_L1, lambda_2_l2=LAMBDA_2_L2, lambda_3_l2=LAMBDA_3_L2, 
                positive_state=True, sumone_state=False, positive_em=True, message_passing=True,
                n_threads=1, verbose=False)
    test_iteration = 10
    model.optimization(test_iteration)
    model.print_error()
    print("errSSM:\n", model.error_m)
    print("emSSM:\n", model.theta_m)
    print("trSSM:\n", model.lambda_m)
    ### plot
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    ## error vs iteration
    err_plot = plt.figure().gca()
    err_k = model.error_m
    err_k = err_k[1:] # to get rid of the initial error which might be too high
    # normalize error values by the max error
    err_plot.plot(range(1, len(err_k)+1), err_k, label="K={}".format(K))
    plt.ylabel("Negative log-likelihood")
    plt.xlabel("Iteration")
    err_plot.xaxis.set_major_locator(MaxNLocator(integer=True))
    err_plot.legend(loc="best")
    plt.savefig("error_vs_iteration.pdf")
    ## error vs time
    plt.clf()
    err_plot = plt.figure().gca()
    err_k = model.error_m
    opt_time_k = model.opt_time_m
    # convert time from second to minute and then hour
    opt_time_k = [t / 60.0 for t in opt_time_k]
    opt_time_k = [t / 60.0 for t in opt_time_k]
    err_k = err_k[1:] # to get rid of the initial error which might be too high
    # normalize error values by the max error
    err_plot.plot(opt_time_k, err_k, label="K={}".format(K))
    plt.ylabel("Negative log-likelihood")
    plt.xlabel("Time (h)")
    err_plot.xaxis.set_major_locator(MaxNLocator(integer=True))
    err_plot.legend(loc="best")
    plt.savefig("error_vs_time.pdf")