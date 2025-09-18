# scipy, py_ecc, galois, sklearn
import logging

from scipy.interpolate import lagrange
import galois
# from tate_bilinear_pairing import eta, ecc, f3m
from py_ecc.bn128 import bn128_curve, bn128_pairing

from py_ecc.fields import (
    bn128_FQ as FQ,
    bn128_FQ2 as FQ2,
    bn128_FQ12 as FQ12,
    bn128_FQP as FQP,
)
import random
import math
import numpy as np
# from poly_utils import PrimeField
from scipy.interpolate import lagrange
from fractions import Fraction
from sklearn.preprocessing import normalize
import time
import pickle


class PrimeField():
    def __init__(self, modulus):
        assert pow(2, modulus, modulus) == 2
        self.modulus = modulus

    def add(self, x, y):
        return (x+y) % self.modulus

    def sub(self, x, y):
        return (x-y) % self.modulus

    def mul(self, x, y):
        return (x*y) % self.modulus

    def exp(self, x, p):
        if p < 0:
            return self.inv(pow(x, -p, self.modulus))
        else:
            return pow(x, p, self.modulus)
    
    def toField(self, x):
        return x % self.modulus
    
    def neg(self, x):
        return -x % self.modulus
    
    
    def moduloMultiplication(self, a, b):
        
        mod = self.modulus
        res = 0; 
        a = a % mod;
        b = b % mod;
        
        while (b):
            if (b & 1):
                res = (res + a) % mod;
            a = (2 * a) % mod;
            
            b >>= 1; 

        return res;

    # Modular inverse using the extended Euclidean algorithm
    def inv(self, a):
        if a == 0:
            return 0
        lm, hm = 1, 0
        low, high = a % self.modulus, self.modulus
        while low > 1:
            r = high//low
            nm, new = hm-lm*r, high-low*r
            lm, low, hm, high = nm, new, lm, low
        return lm % self.modulus

    def multi_inv(self, values):
        partials = [1]
        for i in range(len(values)):
            partials.append(self.mul(partials[-1], values[i] or 1))
        inv = self.inv(partials[-1])
        outputs = [0] * len(values)
        for i in range(len(values), 0, -1):
            outputs[i-1] = self.mul(partials[i-1], inv) if values[i-1] else 0
            inv = self.mul(inv, values[i-1] or 1)
        return outputs

    def div(self, x, y):
        return self.mul(x, self.inv(y))

    # Evaluate a polynomial at a point. The poly must have consecutive orders starting from "init_order".
    def eval_poly_at(self, p, x, init_order = 0):
        y = 0
        power_of_x = 1
        for i, p_coeff in enumerate(p):
            new_y = self.add(self.moduloMultiplication(power_of_x, p_coeff), y)
            y = new_y
            power_of_x = self.mul(power_of_x, x) 
            
        change_order = self.exp(x, int(init_order))

        return self.moduloMultiplication(y, change_order)
    
    def eval_poly_Y(self, coeffs, order, y):
        new_coeffs = []
        for poly in coeffs:
            new_coeffs.append([[self.eval_poly_at(poly[0], y, poly[1])], 0])
        return new_coeffs, order
    
    def eval_poly_X(self, coeffs, init_order, x):
        y = [[0],0]
        power_of_x = 1
        for i, p_coeff in enumerate(coeffs):
            part_poly = self.mul_by_const(p_coeff, power_of_x)
            y = self.add_polys(y[0], part_poly[0], y[1], part_poly[1])
            power_of_x = (power_of_x * x) % self.modulus
        if init_order < 0:
            change_order = self.inv(x ** -init_order)
        else:
            change_order = x ** init_order % self.modulus
        return self.mul_by_const(y, change_order)
        
    # Arithmetic for polynomials
    def add_polys(self, a, b, init_order_a = 0, init_order_b = 0):
        init_order_result = min(init_order_a, init_order_b)
        return ([((a[i-init_order_a+init_order_result] if (i-init_order_a+init_order_result < len(a) and i-init_order_a+init_order_result >= 0) else 0) + 
                  (b[i-init_order_b+init_order_result] if (i-init_order_b+init_order_result < len(b) and i-init_order_b+init_order_result >= 0) else 0))
                % self.modulus for i in range(max(len(a)+init_order_a, len(b)+init_order_b)-init_order_result)], init_order_result)
    
    def sub_polys(self, a, b, init_order_a = 0, init_order_b = 0):
        neg_b = self.mul_by_const([b,init_order_b], -1)
        return self.add_polys(a, neg_b[0], init_order_a, neg_b[1])
        
        # Arithmetic for polynomials
    def add_polys_bivar(self, a, b, init_order_a = 0, init_order_b = 0):
        init_order_result = min(init_order_a, init_order_b)
        result_poly = []
        for i in range(max(len(a)+init_order_a, len(b)+init_order_b)-init_order_result):
            poly_a = a[i-init_order_a+init_order_result] if (i-init_order_a+init_order_result < len(a) and i-init_order_a+init_order_result >= 0) else [[0], 0]
            poly_b = b[i-init_order_b+init_order_result] if (i-init_order_b+init_order_result < len(b) and i-init_order_b+init_order_result >= 0) else [[0], 0]
            result_poly.append(field.add_polys(poly_a[0], poly_b[0], poly_a[1], poly_b[1]))
        return result_poly, init_order_result

    def mul_by_const(self, a, c):
        return [(x*c) % self.modulus for x in a[0]], a[1]
    
    def mul_polys(self, a, b, init_order_a = 0, init_order_b = 0):
        o = [0] * (len(a) + len(b) - 1)
        for i, aval in enumerate(a):
            for j, bval in enumerate(b):
                addon = self.moduloMultiplication(a[i], b[j])
                o[i+j] = (o[i+j] + addon) % self.modulus
        return [x % self.modulus for x in o], init_order_a + init_order_b
    
    def mul_polys_bivar(self, a, b, init_order_a = 0, init_order_b = 0):
        o = [[[0], 0]] * (len(a) + len(b) - 1)
        for i, aval in enumerate(a):
            for j, bval in enumerate(b):
                mul_result = self.mul_polys(a[i][0], b[j][0], a[i][1], b[j][1])
                o[i+j] = self.add_polys(o[i+j][0], mul_result[0], o[i+j][1], mul_result[1])

        return o, init_order_a + init_order_b
    
    def div_polys(self, a, b, init_order_a = 0, init_order_b = 0):
        assert len(a) >= len(b)
        a = [x for x in a]
        o = []
        apos = len(a) - 1
        bpos = len(b) - 1
        diff = apos - bpos
        while diff >= 0:
            quot = self.div(a[apos], b[bpos])
            o.insert(0, quot)
            for i in range(bpos, -1, -1):
                a[diff+i] -= self.moduloMultiplication(b[i], quot)
            apos -= 1
            diff -= 1
        return [x % self.modulus for x in o], init_order_a - init_order_b

    # Build a polynomial from a few coefficients, together with init_order
    def sparse(self, coeff_dict):
        degree = max(coeff_dict.keys()) - min(coeff_dict.keys())
        o = [0] * (degree + 1)
        for k, v in coeff_dict.items():
            o[k - min(coeff_dict.keys())] = v % self.modulus
        return (o, min(coeff_dict.keys()))
    
    def sparse_bivar(self, coeff_dict):
        degree = max(coeff_dict.keys()) - min(coeff_dict.keys())
        o = [[[0], 0]] * (degree + 1)
#         print(o)
        for k, v in coeff_dict.items():
            o[k - min(coeff_dict.keys())] = v
        return (o, min(coeff_dict.keys()))
    
    def lagrange(self, xs, ys):
        fn = None
        for i, y in enumerate(ys):
            xlist = xs[:i]+xs[i+1:]
            xpoly = [list(np.array(np.poly1d(xlist, True).c)[::-1].astype(int)), 0]
            denominator = 1
            for x in xlist:
                denominator *= xs[i] - x
                
            pi = self.mul_by_const(xpoly, self.inv(denominator))
            pi = self.mul_by_const(pi, y)
            
            if fn is None:
                fn = pi
            else:
                fn = self.add_polys(fn[0], pi[0], fn[1], pi[1])
                
        return fn
    
    def isZero(self, a):
        return np.all(np.array(a) == 0)
    
    def dimension_change(self, a, init_order):
        newPoly = []
        for i in a:
            newPoly.append(i[0][0])
        return newPoly, init_order

def decision_tree(x_p, s_j, t_j, l_i, k, h, output):

    #np.set_printoptions(threshold=np.inf)
    n = len(x_p)
    no_of_interal_nodes = 2**h - 1
    
    v_j_p = [[0]*n for _ in range(no_of_interal_nodes)]

    for j in range(no_of_interal_nodes):
        v_j_p[j][s_j[j] - 1] = 1
        
    # Repeat the array
    x_j_p = [x_p for _ in range(no_of_interal_nodes)]

    x_s_j_p = np.multiply(x_j_p, v_j_p)
    
    x_s_j = [sum(array) for array in x_s_j_p]
   
    v_j_p_m = [[] for _ in range(no_of_interal_nodes)]
    v_j_p_m_prime = [[] for _ in range(no_of_interal_nodes)]
    omega_j_p = [np.zeros_like(x_p) for _ in range(no_of_interal_nodes)]
    for j in range(no_of_interal_nodes):
        for p in range(n):
            omega_j_p[j][p] = s_j[j] - p - 1 #starts from 0
            if omega_j_p[j][p] > 0: #bit decomposition of v_j_p_m, assign 0s to v_j_p_m_b_prime
                decomposition = np.array(fixed_length_decomposition((s_j[j] - p), k))
                decomposition = np.squeeze(decomposition)  # Convert to 1D array if needed
                v_j_p_m[j] = np.concatenate((v_j_p_m[j], decomposition)).astype(int)
                v_j_p_m_prime[j] = np.concatenate((v_j_p_m_prime[j], np.zeros(k))).astype(int)
            elif omega_j_p[j][p] < 0:
                decomposition_prime = np.array(fixed_length_decomposition((s_j[j] - p - 2), k))
                decomposition_prime = np.squeeze(decomposition_prime)  # Convert to 1D array if needed
                v_j_p_m_prime[j] = np.concatenate((v_j_p_m_prime[j], decomposition_prime)).astype(int)
                v_j_p_m[j] = np.concatenate((v_j_p_m[j], np.zeros(k))).astype(int)
            else:
                v_j_p_m_prime[j] = np.concatenate((v_j_p_m_prime[j], np.zeros(k))).astype(int)
                v_j_p_m[j] = np.concatenate((v_j_p_m[j], np.zeros(k))).astype(int)

    v_j_p_m_b = np.zeros_like(v_j_p_m) #v_j_p_m_b - 2**(m - 1)
    v_j_p_m_b_prime = np.zeros_like(v_j_p_m_prime) 
    sigma_v_j_p_m = np.zeros_like(v_j_p)
    sigma_v_j_p_m_prime = np.zeros_like(v_j_p)

    for j in range(no_of_interal_nodes):
        for p in range(n):
            for m in range(k):
                v_j_p_m_b[j][p*k + m] = v_j_p_m[j][p*k + m] - 2**(m)
                v_j_p_m_b_prime[j][p*k + m] = v_j_p_m_prime[j][p*k + m] + 2**(m)
                sigma_v_j_p_m[j][p] = sigma_v_j_p_m[j][p] + v_j_p_m[j][p*k + m]
                sigma_v_j_p_m_prime[j][p] = sigma_v_j_p_m_prime[j][p] + v_j_p_m_prime[j][p*k + m]
    
    theta_j_p = np.multiply((sigma_v_j_p_m - omega_j_p - 1),(sigma_v_j_p_m_prime - omega_j_p + 1))

    w_j = x_s_j - t_j
    c_j = [[0] for _ in range(no_of_interal_nodes)]
    b_j_m = [[] for _ in range(no_of_interal_nodes)]
    b_j_m_prime = [[] for _ in range(no_of_interal_nodes)]
    for j in range(no_of_interal_nodes):
        if w_j[j] > 0:
            b_j_m[j] = np.zeros(k).astype(int)
            decomposition = fixed_length_decomposition(w_j[j] - 1,k)
            decomposition = np.squeeze(decomposition)
            b_j_m_prime[j].extend(decomposition)
            c_j[j] = np.array([1])
        else: 
            b_j_m_prime[j] = np.zeros(k).astype(int)
            decomposition2 = fixed_length_decomposition(w_j[j],k)
            decomposition2 = np.squeeze(decomposition2)
            b_j_m[j].extend(decomposition2)
    
    b_j_m_b = np.zeros_like(b_j_m)
    b_j_m_prime_b = np.zeros_like(b_j_m_prime)
    for j in range(no_of_interal_nodes):
        for m in range(k):
            b_j_m_b[j][m] = b_j_m[j][m] + 2**m
            b_j_m_prime_b[j][m] = b_j_m_prime[j][m] - 2**m

    sigma_b_j_m = [sum(array) for array in b_j_m]
    sigma_b_j_m_prime = [sum(array) for array in b_j_m_prime]

 
    beta_i_a = [[] for _ in range(2**h,2**(h+1))]
    for i in range(2**h,2**(h+1)): #stops before 2**(h+1) - 1
        beta_i_a[i-2**h] = fixed_length_binary_decomposition(i,h+1)

    two_power_a_minus_1 = np.zeros_like(beta_i_a)       
        
    for i in range(2**h,2**(h+1)):
        for alpha in range(h + 1):
            two_power_a_minus_1[i - 2**h][alpha] = 2**alpha
            
    u_i_a = beta_i_a * two_power_a_minus_1


    c_theta_i_d = [[0]*h for _ in range(2**h,2**(h+1))]
    theta_i_d = np.zeros_like(c_theta_i_d)
    one_minus_c_a = np.zeros_like(c_theta_i_d)
    for i in range(2**h,2**(h+1)):
        for d in range(h):
            for f in range(d + 1):
                value = beta_i_a[i - 2**h][h - f] * 2**(d - f)
                theta_i_d[i - 2**h][d] = theta_i_d[i - 2**h][d] + value
                c_theta_i_d[i - 2**h][d] = c_theta_i_d[i - 2**h][d] + value
            c_theta_i_d[i - 2**h][d] = c_j[c_theta_i_d[i - 2**h][d] - 1][0]
            one_minus_c_a[i - 2**h][d] = 1 - c_theta_i_d[i - 2**h][d]

    beta_i_h_minus_d_plus_one = np.zeros_like(c_theta_i_d)
    one_minus_beta_i_h_minus_d_plus_one = np.zeros_like(beta_i_h_minus_d_plus_one)
    for i in range(2**h,2**(h+1)):
        for d in range(h):
            beta_i_h_minus_d_plus_one[i - 2**h][d] = beta_i_a[i - 2**h][h - d - 1]
            one_minus_beta_i_h_minus_d_plus_one[i - 2**h][d] = 1 - beta_i_a[i - 2**h][h - d - 1]
    
    zeta_i_d = c_theta_i_d * beta_i_h_minus_d_plus_one
    zeta_i_d_prime = one_minus_c_a * one_minus_beta_i_h_minus_d_plus_one
    
    gamma_i_d = zeta_i_d + zeta_i_d_prime
    z_i_d_l = [[0]*(h + 1) for _ in range(2**h,2**(h+1))]
    z_i_d_minus_one = np.zeros_like(gamma_i_d)
    z_i_d = np.zeros_like(gamma_i_d)
    for i in range(2**h,2**(h+1)):
        for d in range(1, h + 1):
            z_i_d_l[i - 2**h][0] = 1
            z_i_d_l[i - 2**h][d] = gamma_i_d[i - 2**h][d - 1] * z_i_d_l[i - 2**h][d - 1]
            z_i_d_minus_one[i - 2**h][d - 1] = z_i_d_l[i - 2**h][d - 1]
            z_i_d[i - 2**h][d - 1] = z_i_d_l[i - 2**h][d]

    ###
    # Repeat c_j_i_d
    c_j_i_d = [[c_j for _ in range(h)] for _ in range(2**h,2**(h+1))]
    mu_j_i_d = [[[0 for _ in range(2**h - 1)] for _ in range(h)] for _ in range(2**h,2**(h+1))]#j,d,i

    mu_j_i_d_m = [[[[] for _ in range(2**h - 1)] for _ in range(h)] for _ in range(2**h,2**(h+1))]

    mu_j_i_d_m_prime = [[[[] for _ in range(2**h - 1)] for _ in range(h)] for _ in range(2**h,2**(h+1))]
    
    Omega_j_i_d = np.zeros_like(mu_j_i_d)
    for j in range(1,2**h):
        for i in range(2**h,2**(h+1)):
            for d in range(h):
                Omega_j_i_d[i - 2**h][d][j-1] = theta_i_d[i - 2**h][d] - j                 
                if Omega_j_i_d[i - 2**h][d][j-1] > 0:
                    decomposition = np.array(fixed_length_decomposition((theta_i_d[i - 2**h][d] - j + 1), k))
                    decomposition = np.squeeze(decomposition)  # Convert to 1D array if needed
                    mu_j_i_d_m[i - 2**h][d][j-1] = np.concatenate((mu_j_i_d_m[i - 2**h][d][j-1], decomposition)).astype(int)
                    mu_j_i_d_m_prime[i - 2**h][d][j-1] = np.concatenate((mu_j_i_d_m_prime[i - 2**h][d][j-1], np.zeros(k))).astype(int)
                elif Omega_j_i_d[i - 2**h][d][j-1] < 0:
                    decomposition = np.array(fixed_length_decomposition((theta_i_d[i - 2**h][d] - j - 1), k))
                    decomposition = np.squeeze(decomposition)  # Convert to 1D array if needed
                    mu_j_i_d_m_prime[i - 2**h][d][j-1] = np.concatenate((mu_j_i_d_m_prime[i - 2**h][d][j-1], decomposition)).astype(int)
                    mu_j_i_d_m[i - 2**h][d][j-1] = np.concatenate((mu_j_i_d_m[i - 2**h][d][j-1], np.zeros(k))).astype(int)
                else: #=0
                    mu_j_i_d[i - 2**h][d][j-1] = 1
                    mu_j_i_d_m[i - 2**h][d][j-1] = np.concatenate((mu_j_i_d_m[i - 2**h][d][j-1], np.zeros(k))).astype(int)
                    mu_j_i_d_m_prime[i - 2**h][d][j-1] = np.concatenate((mu_j_i_d_m_prime[i - 2**h][d][j-1], np.zeros(k))).astype(int)
                    
    mu_j_i_d_m_b = np.zeros_like(mu_j_i_d_m)
    mu_j_i_d_m_b_prime = np.zeros_like(mu_j_i_d_m_prime)
    for j in range(1,2**h):
        for i in range(2**h,2**(h+1)):
            for d in range(h):            
                for m in range(1,k+1):
                    mu_j_i_d_m_b[i - 2**h][d][j-1][m-1] = mu_j_i_d_m[i - 2**h][d][j-1][m-1] - 2**(m - 1)
                    mu_j_i_d_m_b_prime[i - 2**h][d][j-1][m-1] = mu_j_i_d_m_prime[i - 2**h][d][j-1][m-1] + 2**(m - 1)
                    
    #flattened_list = [item for sublist1 in three_d_array for sublist2 in sublist1 for item in sublist2]

    c_theta_j_i_d = np.array(c_j_i_d).squeeze(-1) * np.array(mu_j_i_d)

    #c_theta_i_d_two = np.sum(c_theta_j_i_d, axis=2)
    sigma_mu_j_i_d_m = np.sum(mu_j_i_d_m, axis=3)
    sigma_mu_j_i_d_m_prime = np.sum(mu_j_i_d_m_prime, axis=3)
    eta_j_i_d = (sigma_mu_j_i_d_m - Omega_j_i_d - 1)*(sigma_mu_j_i_d_m_prime - Omega_j_i_d + 1)

    
    z_i = [[0] for _ in range(2**h,2**(h+1))]
    for i in range(2**h,2**(h+1)):
        z_i[i - 2**h] = z_i_d[i - 2**h][h - 1] 
    
    epsilon_i = z_i * l_i

    #flattening lists
    x_j_p = np.array([item for sublist in x_j_p for item in sublist])
    v_j_p = np.array([item for sublist in v_j_p for item in sublist])
    v_j_p_m = np.array([item for sublist in v_j_p_m for item in sublist])
    v_j_p_m_prime = np.array([item for sublist in v_j_p_m_prime for item in sublist])
    omega_j_p = np.array([item for sublist in omega_j_p for item in sublist])
    sigma_v_j_p_m = np.array([item for sublist in sigma_v_j_p_m for item in sublist])
    theta_j_p = np.array([item for sublist in theta_j_p for item in sublist])
    b_j_m = np.array([item for sublist in b_j_m for item in sublist])
    c_j = np.array([item for sublist in c_j for item in sublist])
    b_j_m_prime = np.array([item for sublist in b_j_m_prime for item in sublist])
    beta_i_a = np.array([item for sublist in beta_i_a for item in sublist])
    c_theta_i_d = np.array([item for sublist in c_theta_i_d for item in sublist])
    gamma_i_d = np.array([item for sublist in gamma_i_d for item in sublist])
    v_j_p_m_b = np.array([item for sublist in v_j_p_m_b for item in sublist])
    v_j_p_m_b_prime = np.array([item for sublist in v_j_p_m_b_prime for item in sublist])
    sigma_v_j_p_m_prime = np.array([item for sublist in sigma_v_j_p_m_prime for item in sublist])
    b_j_m_b = np.array([item for sublist in b_j_m_b for item in sublist])
    b_j_m_prime_b = np.array([item for sublist in b_j_m_prime_b for item in sublist])
    two_power_a_minus_1 = np.array([item for sublist in two_power_a_minus_1 for item in sublist])
    beta_i_h_minus_d_plus_one = np.array([item for sublist in beta_i_h_minus_d_plus_one for item in sublist])
    z_i_d_minus_one = np.array([item for sublist in z_i_d_minus_one for item in sublist])
    x_s_j_p = np.array([item for sublist in x_s_j_p for item in sublist])
    u_i_a = np.array([item for sublist in u_i_a for item in sublist])
    zeta_i_d = np.array([item for sublist in zeta_i_d for item in sublist])
    zeta_i_d_prime = np.array([item for sublist in zeta_i_d_prime for item in sublist])
    c_j_i_d = np.array([item for sublist1 in c_j_i_d for sublist2 in sublist1 for sublist3 in sublist2 for item in sublist3])
    mu_j_i_d = np.array([item for sublist1 in mu_j_i_d for sublist2 in sublist1 for item in sublist2])
    c_theta_j_i_d = np.array([item for sublist1 in c_theta_j_i_d for sublist2 in sublist1 for item in sublist2])
    mu_j_i_d_m = np.array([item for sublist1 in mu_j_i_d_m for sublist2 in sublist1 for sublist3 in sublist2 for item in sublist3])
    mu_j_i_d_m_b = np.array([item for sublist1 in mu_j_i_d_m_b for sublist2 in sublist1 for sublist3 in sublist2 for item in sublist3])
    mu_j_i_d_m_prime = np.array([item for sublist1 in mu_j_i_d_m_prime for sublist2 in sublist1 for sublist3 in sublist2 for item in sublist3])
    mu_j_i_d_m_b_prime = np.array([item for sublist1 in mu_j_i_d_m_b_prime for sublist2 in sublist1 for sublist3 in sublist2 for item in sublist3])
    Omega_j_i_d = np.array([item for sublist1 in Omega_j_i_d for sublist2 in sublist1 for item in sublist2])
    sigma_mu_j_i_d_m = np.array([item for sublist1 in sigma_mu_j_i_d_m for sublist2 in sublist1 for item in sublist2])
    sigma_mu_j_i_d_m_prime = np.array([item for sublist1 in sigma_mu_j_i_d_m_prime for sublist2 in sublist1 for item in sublist2])
    eta_j_i_d = np.array([item for sublist1 in eta_j_i_d for sublist2 in sublist1 for item in sublist2])
    z_i_d = np.array([item for sublist in z_i_d for item in sublist])
    theta_i_d = np.array([item for sublist in theta_i_d for item in sublist])

    a_upper = np.concatenate( [x_j_p,
                                v_j_p,
                                v_j_p_m,
                                v_j_p_m_prime,
                                omega_j_p,
                                (sigma_v_j_p_m - omega_j_p - 1),
                                theta_j_p,
                                b_j_m,
                                c_j,
                                (sigma_b_j_m - w_j),
                                b_j_m_prime,
                                (sigma_b_j_m_prime - w_j + 1),
                                c_theta_i_d,
                                (1 - c_theta_i_d),
                                gamma_i_d,
                                c_j_i_d,
                                mu_j_i_d,
                                mu_j_i_d_m,
                                mu_j_i_d_m_prime,
                                Omega_j_i_d,
                                (sigma_mu_j_i_d_m - Omega_j_i_d - 1),
                                eta_j_i_d,
                                l_i
                               ]
                            ).astype(int)
    
    b_upper = np.concatenate( [v_j_p,
                               (1 - v_j_p),
                                v_j_p_m_b,
                                v_j_p_m_b_prime,
                                v_j_p,
                                (sigma_v_j_p_m_prime - omega_j_p + 1),
                                (1 - v_j_p),
                                b_j_m_b,
                                (1 - c_j),
                                (1 - c_j),
                                b_j_m_prime_b,
                                c_j,
                                beta_i_h_minus_d_plus_one,
                                (1 - beta_i_h_minus_d_plus_one),
                                z_i_d_minus_one,
                                mu_j_i_d,
                                (1 - mu_j_i_d),
                                mu_j_i_d_m_b,
                                mu_j_i_d_m_b_prime,
                                mu_j_i_d,
                                (sigma_mu_j_i_d_m_prime - Omega_j_i_d + 1),
                                (1 - mu_j_i_d),
                                z_i
                               ]
                            ).astype(int)
    
    c_upper = np.concatenate( [x_s_j_p,
                                np.zeros_like(v_j_p),
                                np.zeros_like(v_j_p_m),
                                np.zeros_like(v_j_p_m_prime),
                                np.zeros_like(omega_j_p),
                                theta_j_p,
                                np.zeros_like(theta_j_p),
                                np.zeros_like(b_j_m),
                                np.zeros_like(c_j),
                                np.zeros_like(c_j),
                                np.zeros_like(b_j_m_prime),
                                np.zeros_like(c_j),
                                zeta_i_d,
                                zeta_i_d_prime,
                                z_i_d,
                                c_theta_j_i_d,
                                np.zeros_like(mu_j_i_d),
                                np.zeros_like(mu_j_i_d_m),
                                np.zeros_like(mu_j_i_d_m_prime),
                                np.zeros_like(mu_j_i_d),
                                eta_j_i_d,
                                np.zeros_like(mu_j_i_d),
                                epsilon_i
                                ]
                            ).astype(int)
    
    # Ext input
    y = np.array([output])
    a_lower = np.concatenate([x_p, s_j, t_j, l_i,y]).astype(int)
    b_lower = np.zeros_like(a_lower).astype(int)
    c_lower = np.zeros_like(a_lower).astype(int)
        
    a = np.concatenate([a_upper, a_lower])
    b = np.concatenate([b_upper, b_lower])
    c = np.concatenate([c_upper, c_lower])

    # linear constraints. k_q is a list of scalar, u_q, v_q, w_q are list of numpy arrays
    k_q = []
    u_q = []
    v_q = []
    w_q = []
    
    
    #\sum_{p = 1}^{n} x_{[s_{[j,p]}]})_{j = 1}^{2^h - 1} = (x_{[s_{[j]}]})_{j = 1}^{2^h - 1}
    for j in range(2**h-1):
        k_2 = 1
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        #t_j
        u[(2**h - 1) * (5*n + 2*n*k + 2*k + 4 + (2**h)*5*h + (2**h)*2*h*k) + 2**h*(3*h+1) + n + j] = -2
        #sigma_b_j_m - w_j
        u[(2**h - 1) * (5*n + 2*n*k + k + 1) + j] = 1
        #sigma_b_j_m' - w_j + 1
        u[(2**h - 1) * (5*n + 2*n*k + 2*k + 2) + j] = 1
        for p in range(n):
            # x_s_j_p
            w[j*n+p] = 2      
        for m in range(k):
            # b_j_m
            u[(2**h - 1) * (5*n + 2*n*k) + j*k + m] = -1
            # b_j_m_prime
            u[(2**h - 1) * (5*n + 2*n*k + k + 2) + j*k + m] = -1

        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)


    #\(\sum_{j = 1}^{2^h-1} c_{[\theta_{[i,d,j]}]})_{i = 2^h, d = 1}^{2^{h+1} - 1, h} = (c_{[\theta_{[i,d]}]})_{i = 2^h, d = 1}^{2^{h+1} - 1, h}
    for ind in range(len(c_theta_i_d)):
        k_2 = 0
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        u[(2**h - 1) * (5*n + 2*n*k + 2*k + 3)+ ind] = 1
        for j in range(2**h-1):
            w[(2**h - 1) * (5*n + 2*n*k + 2*k + 3) + (2**h)*3*h + ind*(2**h-1) + j] = -1
 
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
    
    
    #(z_{[i,h]})_{i = 2^h}^{2^{h + 1} - 1} = (z_{[i]})_{i = 2^h}^{2^{h + 1} - 1}
    for i in range(2**h,2**(h+1)):
        k_2 = 0
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        v[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*5*h + (2**h)*2*h*k) + (2**h)*3*h+i-2**h] = 1 
        w[(2**h - 1) * (5*n + 2*n*k + 2*k + 3) + (2**h)*2*h + (i-2**h+1)*h - 1] = -1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
    
    
    # \sum_{i = 2^h}^{2^{h + 1} - 1} \epsilon_{[i]} = y
    k_2 = 0
    u = np.zeros_like(a)
    v = np.zeros_like(b)
    w = np.zeros_like(c)
    
    u[(2**h - 1) * (5*n + 2*n*k + 2*k + 5 + (2**h)*5*h + (2**h)*2*h*k) + 2**h*(3*h+2) + n] = 1 
    w[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*5*h + (2**h)*2*h*k) + (2**h)*3*h:(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*5*h + (2**h)*2*h*k) + 2**h*(3*h+1)] = -1
    
    k_q.append(k_2)
    u_q.append(u)
    v_q.append(v)
    w_q.append(w)
    
    
    # (s_{[j]} - p = \omega_{[j,p]})_{j = 1, p = 1}^{2^h - 1, n}
    for j in range(2**h-1):
        for p in range(n):
            k_2 = p+1
            u = np.zeros_like(a)
            v = np.zeros_like(b)
            w = np.zeros_like(c)
            #omega_j_p
            u[(2**h - 1) * (2*n + 2*n*k) + j*n + p] = -1 
            #l_i
            u[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*5*h + (2**h)*2*h*k) + 2**h*(3*h+1) + n + j] = 1            
            k_q.append(k_2)
            u_q.append(u)
            v_q.append(v)
            w_q.append(w)
            

    # (\theta_{[i,d]} - j = \Omega_{[j,i,d]})_{j = 1, i = 2^h, d = 1}^{2^h-1, 2^{h+1}-1, h}               
    for iandd in range(2**h*h):
        for j in range(2**h-1):    
            k_2 = theta_i_d[iandd] - j - 1
            u = np.zeros_like(a)
            v = np.zeros_like(b)
            w = np.zeros_like(c)
            # Omega_j_i_d
            u[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*2*h + (2**h)*2*h*k) + (2**h)*3*h + iandd*(2**h-1) + j] = 1
            k_q.append(k_2)
            u_q.append(u)
            v_q.append(v)
            w_q.append(w)

    
    #(\zeta{[i,d]} + \zeta'{[i,d]} = \gamma_{[i,d]})_{i = 2^h, d = 1}^{2^{h+1}-1, h}
    for ind in range(len(zeta_i_d)):
        k_2 = 0
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        
        u[(2**h - 1) * (5*n + 2*n*k + 2*k + 3) + (2**h)*2*h + ind] = 1 
        w[(2**h - 1) * (5*n + 2*n*k + 2*k + 3) + ind] = -1
        w[(2**h - 1) * (5*n + 2*n*k + 2*k + 3) + (2**h)*h + ind] = -1
        
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
        
    
    # (x_{[p]})_{p = 1}^{n} \) match \( (x_{[p]})_{j = 1, p = 1}^{2^h - 1, n}
    for j in range(2**h-1):
        for ind in range(len(x_p)):
            k_2 = 0
            u = np.zeros_like(a)
            v = np.zeros_like(b)
            w = np.zeros_like(c)
            u[j*len(x_p)+ind] = 1
            u[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*5*h + (2**h)*2*h*k) + 2**h*(3*h+1) + ind] = -1
            k_q.append(k_2)
            u_q.append(u)
            v_q.append(v)
            w_q.append(w)

    
    # (l_{[i]})_{i = 2^h}^{2^{h +1} - 1} \) match \( (l_{[i]})_{i = 2^h}^{2^{h +1} - 1} 
    for ind in range(len(l_i)):
        k_2 = 0
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        u[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*5*h + (2**h)*2*h*k) + (2**h)*3*h + ind] = 1
        u[(2**h - 1) * (5*n + 2*n*k + 2*k + 5 + (2**h)*5*h + (2**h)*2*h*k) + 2**h*(3*h+1) + n + ind] = -1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
    
    
    # (c_{[j]})_{j = 1}^{2^h - 1}\) match \((c_{[j]})_{j = 1, i=2^h, d = 1 }^{2^h - 1, 2^{h+1} - 1, h}
    for iandd in range((2**(h+1)-2**h)*h):
        for ind in range(len(c_j)):
            k_2 = 0
            u = np.zeros_like(a)
            v = np.zeros_like(b)
            w = np.zeros_like(c)
            #c_j
            u[(2**h - 1) * (5*n + 2*n*k + k) + ind] = 1
            # c_j_i_d
            u[(2**h - 1) * (5*n + 2*n*k + 2*k + 3) + (2**h)*3*h + iandd*len(c_j) + ind] = -1
            k_q.append(k_2)
            u_q.append(u)
            v_q.append(v)
            w_q.append(w)
    
    
    ### (\beta{[i,h-d+1]})_{i = 2^h, d = 1}^{2^{h +1} - 1, h} \) is correct
    for i in range(2**h,2**(h+1)):   
        decomposition = np.array(fixed_length_binary_decomposition((i), k))
        for d in range(h):
            k_2 = decomposition[h-d-1]
            u = np.zeros_like(a)
            v = np.zeros_like(b)
            w = np.zeros_like(c)
            v[(2**h - 1) * (5*n + 2*n*k + 2*k + 3) + (i-2**h)*h + d] = 1
            k_q.append(k_2)
            u_q.append(u)
            v_q.append(v)
            w_q.append(w)


    # Copies of \( (v_{[j,p]})_{j = 1, p = 1}^{2^h -1 , n} \)
    for ind in range(len(v_j_p)):
        k_2 = 0
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        u[(2**h-1)*n + ind] = 1
        v[ind] = -1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
    
    
    # Copies of \( (v_{[j,p]})_{j = 1, p = 1}^{2^h -1 , n} \)
    for ind in range(len(v_j_p)):
        k_2 = 0
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        v[ind] = -1
        v[(2**h - 1) * (2*n + 2*n*k) + ind] = 1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
    
    
    # Copies of \( (1 - v_{[j,p]})_{j = 1, p = 1}^{2^h -1 , n} \)
    for ind in range(len(v_j_p)):
        k_2 = 1
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        v[ind] = 1
        v[(2**h-1)*n + ind] = 1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
    
    
    # Copies of \( (1 - v_{[j,p]})_{j = 1, p = 1}^{2^h -1 , n} \)
    for ind in range(len(v_j_p)):
        k_2 = 0
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        v[(2**h-1)*n + ind] = -1
        v[(2**h - 1) * (4*n + 2*n*k) + ind] = 1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)

    
    for jandp in range((2**h-1)*n):
        for m in range(k):
            k_2 = 2**m
            u = np.zeros_like(a)
            v = np.zeros_like(b)
            w = np.zeros_like(c)
            u[(2**h-1)*2*n + jandp*k + m] = 1
            v[(2**h-1)*2*n + jandp*k + m] = -1
            k_q.append(k_2)
            u_q.append(u)
            v_q.append(v)
            w_q.append(w)

    
    
    #  Copies of \((v'_{[j,p,m]} + 2^{m-1})_{j = 1, p = 1, m = 1}^{2^h -1, n, k}\)
    for jandp in range((2**h-1)*n):
        for m in range(k):
            k_2 = -2**m
            u = np.zeros_like(a)
            v = np.zeros_like(b)
            w = np.zeros_like(c)
            u[(2**h - 1) * (2*n + n*k) + jandp*k + m] = 1
            v[(2**h - 1) * (2*n + n*k) + jandp*k + m] = -1
            k_q.append(k_2)
            u_q.append(u)
            v_q.append(v)
            w_q.append(w)
    
    
    # Copies of \((\Theta{[j,p]})_{j = 1, p = 1}^{2^h - 1, n}\) 
    for ind in range(len(theta_j_p)):
        k_2 = 0
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        u[(2**h - 1) * (4*n + 2*n*k) + ind] = 1
        w[(2**h - 1) * (3*n + 2*n*k) + ind] = -1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
    
    
    # Copies of \((\sum_{m=1}^{k} v_{[j,p,m]} - \omega_{[j,p]} - 1)_{j = 1, p = 1}^{2^h - 1, n}\)
    for jandp in range((2**h-1)*n):
        k_2 = -1
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        for m in range(k):
            u[(2**h-1)*2*n + jandp*k + m] = -1
            
        u[(2**h - 1) * (3*n + 2*n*k) + jandp] = 1            
        u[(2**h - 1) * (2*n + 2*n*k) + jandp] = 1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)

    
    # Copies of \((\sum_{m=1}^{k} v'_{[j,p,m]} - \omega_{[j,p]} + 1)_{j = 1, p = 1}^{2^h - 1, n}\)
    for jandp in range((2**h-1)*n):
        k_2 = 1
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        for m in range(k):
            u[(2**h - 1) * (2*n + n*k) + jandp*k + m] = -1
            
        v[(2**h - 1) * (3*n + 2*n*k) + jandp] = 1            
        u[(2**h - 1) * (2*n + 2*n*k) + jandp] = 1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
    
    
    # Copies of \((b_{[j,m]} + 2^{m - 1})_{j = 1, m = 1}^{2^h - 1, k}\)
    for j in range(2**h-1):
        for m in range(k):
            k_2 = 2**m
            u = np.zeros_like(a)
            v = np.zeros_like(b)
            w = np.zeros_like(c)
            v[(2**h - 1) * (5*n + 2*n*k) + j*k+m] = 1            
            u[(2**h - 1) * (5*n + 2*n*k) + j*k+m] = -1
            k_q.append(k_2)
            u_q.append(u)
            v_q.append(v)
            w_q.append(w)
            
    
    # Copies of \((c_{[j]})_{j = 1}^{2^h - 1}\)
    for ind in range(len(c_j)):
        k_2 = 0
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        u[(2**h - 1) * (5*n + 2*n*k + k) + ind] = 1
        v[(2**h - 1) * (5*n + 2*n*k + 2*k + 2) + ind] = -1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
    
    
    # Copies of \((1 - c_{[j]})_{j = 1}^{2^h - 1}\) 
    for ind in range(len(c_j)):
        k_2 = 1
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        u[(2**h - 1) * (5*n + 2*n*k + k) + ind] = 1
        v[(2**h - 1) * (5*n + 2*n*k + k) + ind] = 1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
    
    
    # Copies of \((1 - c_{[j]})_{j = 1}^{2^h - 1}\) 
    for ind in range(len(c_j)):
        k_2 = 0
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        v[(2**h - 1) * (5*n + 2*n*k + k) + ind] = 1
        v[(2**h - 1) * (5*n + 2*n*k + k + 1) + ind] = -1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
    
    
    # Copies of \((b'_{[j,m]} - 2^{m - 1})_{j = 1, m = 1}^{2^h - 1, k}\)
    for j in range(2**h-1):
        for m in range(k):
            k_2 = 2**m
            u = np.zeros_like(a)
            v = np.zeros_like(b)
            w = np.zeros_like(c)
            u[(2**h - 1) * (5*n + 2*n*k + k + 2) + j*k+m] = 1            
            v[(2**h - 1) * (5*n + 2*n*k + k + 2) + j*k+m] = -1
            k_q.append(k_2)
            u_q.append(u)
            v_q.append(v)
            w_q.append(w)
    
    
    # Copies of \((1 - c_{[\theta_{[i,d]}]})_{i = 2^h, d = 1}^{2^{h+1} - 1, h}\)
    for ind in range(len(c_theta_i_d)):
        k_2 = 1
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        u[(2**h - 1) * (5*n + 2*n*k + 2*k + 3) + ind] = 1
        u[(2**h - 1) * (5*n + 2*n*k + 2*k + 3) + (2**h)*h + ind] = 1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
    
    
    # Copies of \((1 - \beta_{[i,h - d + 1]})_{i = 2^h, d = 1}^{2^{h+1} - 1, h}\)
    for ind in range(len(beta_i_h_minus_d_plus_one)):
        k_2 = 1
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        v[(2**h - 1) * (5*n + 2*n*k + 2*k + 3) + ind] = 1
        v[(2**h - 1) * (5*n + 2*n*k + 2*k + 3) + (2**h)*h + ind] = 1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
    
    
    # (z_{[i,0]})_{i = 2^h}^{2^{h+1} - 1} = 1\) 
    for i in range(2**h,2**(h+1)):
        k_2 = 1
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        v[(2**h - 1) * (5*n + 2*n*k + 2*k + 3) + (2**h)*2*h + (i-2**h)*h] = 1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
    
    
    # Copies of \((z_{[i,d]})_{i = 2^h, d = 1}^{2^{h+1} - 1, h}\)
    for i in range(2**h,2**(h+1)):
        for d in range(h-1):
            k_2 = 0
            u = np.zeros_like(a)
            v = np.zeros_like(b)
            w = np.zeros_like(c)
            #zid-1
            v[(2**h - 1) * (5*n + 2*n*k + 2*k + 3) + (2**h)*2*h + (i-2**h)*h + d + 1] = 1
            #zid
            w[(2**h - 1) * (5*n + 2*n*k + 2*k + 3) + (2**h)*2*h + (i-2**h)*h + d] = -1
            k_q.append(k_2)
            u_q.append(u)
            v_q.append(v)
            w_q.append(w)
            
    
    # Copies of \((\mu_{[j,i,d]})_{j = 1, i = 2^h, d = 1}^{2^h-1, 2^{h+1} - 1, h}\)
    for ind in range(len(mu_j_i_d)):
        k_2 = 0
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        v[(2**h - 1) * (5*n + 2*n*k + 2*k + 3) + (2**h)*3*h + ind] = 1
        u[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*h) + (2**h)*3*h + ind] = -1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
        
    
    # Copies of \((\mu_{[j,i,d]})_{j = 1, i = 2^h, d = 1}^{2^h-1, 2^{h+1} - 1, h}\)
    for ind in range(len(mu_j_i_d)):
        k_2 = 0
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        v[(2**h - 1) * (5*n + 2*n*k + 2*k + 3) + (2**h)*3*h + ind] = 1
        v[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*2*h + (2**h)*2*h*k) + (2**h)*3*h + ind] = -1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
    
    
    # Copies of \((1 - \mu_{[j,i,d]})_{j = 1, i = 2^h, d = 1}^{2^h-1, 2^{h+1} - 1, h}\)
    for ind in range(len(mu_j_i_d)):
        k_2 = 1
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        v[(2**h - 1) * (5*n + 2*n*k + 2*k + 3) + (2**h)*3*h + ind] = 1
        v[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*h) + (2**h)*3*h + ind] = 1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
    
    
    # Copies of \((1 - \mu_{[j,i,d]})_{j = 1, i = 2^h, d = 1}^{2^h-1, 2^{h+1} - 1, h}\)
    for ind in range(len(mu_j_i_d)):
        k_2 = 0
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        v[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*h) + (2**h)*3*h + ind] = 1
        v[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*4*h + (2**h)*2*h*k) + (2**h)*3*h + ind] = -1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
    
    
    #  Copies of \((\mu_{[j,i,d,m]} - 2^{m-1})_{j = 1, i = 2^h, d = 1, m =1}^{2^h-1, 2^{h+1} - 1, h, k}\)
    for jid in range((2**h-1)*(2**h)*h):
        for m in range(k):
            k_2 = 2**m
            u = np.zeros_like(a)
            v = np.zeros_like(b)
            w = np.zeros_like(c)
            u[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*2*h) + (2**h)*3*h+ jid*k+m] = 1            
            v[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*2*h) + (2**h)*3*h + jid*k+m] = -1
            k_q.append(k_2)
            u_q.append(u)
            v_q.append(v)
            w_q.append(w)
    
    
    # Copies of \((\mu'_{[j,i,d,m]} + 2^{m-1})_{j = 1, i = 2^h, d = 1, m =1}^{2^h-1, 2^{h+1} - 1, h, k}\)
    for jid in range((2**h-1)*(2**h)*h):
        for m in range(k):
            k_2 = -2**m
            u = np.zeros_like(a)
            v = np.zeros_like(b)
            w = np.zeros_like(c)
            u[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*2*h + (2**h)*h*k) + (2**h)*3*h + jid*k+m] = 1            
            v[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*2*h + (2**h)*h*k) + (2**h)*3*h + jid*k+m] = -1
            k_q.append(k_2)
            u_q.append(u)
            v_q.append(v)
            w_q.append(w)
    
    
    # Copies of (\sum_{m=1}^{k} \mu_{[j,i,d,m]} - \Omega_{[j,i,d]} - 1)_{j = 1, i = 2^h, d = 1}^{2^h-1, 2^{h+1} - 1, h} 
    for ind in range((2**h-1)*(2**h)*h):
        k_2 = -1
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        #Omega_j_i_d
        u[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*2*h + (2**h)*2*h*k) + (2**h)*3*h + ind] = 1 
        # sigmamujid
        u[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*3*h + (2**h)*2*h*k) + (2**h)*3*h + ind] = 1
        for m in range(k):
            # summujidm
            u[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*2*h) + (2**h)*3*h + ind*k + m] = -1

    # Copies of (\sum_{m=1}^{k} \mu'_{[j,i,d,m]} - \Omega_{[j,i,d]} + 1)_{j = 1, i = 2^h, d = 1}^{2^h-1, 2^{h+1} - 1, h} \)
    for ind in range((2**h-1)*(2**h)*h):
        k_2 = 1
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        #Omega_j_i_d
        u[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*2*h + (2**h)*2*h*k) + (2**h)*3*h + ind] = 1 
        # sigmamujid'
        v[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*3*h + (2**h)*2*h*k) + (2**h)*3*h + ind] = 1
        for m in range(k):
            # summujidm'
            u[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*2*h + (2**h)*h*k) + (2**h)*3*h + ind*k + m] = -1
    
    # Copies of \((\eta_{[j,i,d]})_{j = 1, i = 2^h, d = 1}^{2^h-1, 2^{h+1} - 1, h}\) 
    for ind in range(len(eta_j_i_d)):
        k_2 = 0
        u = np.zeros_like(a)
        v = np.zeros_like(b)
        w = np.zeros_like(c)
        w[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*3*h + (2**h)*2*h*k) + (2**h)*3*h + ind] = 1
        u[(2**h - 1) * (5*n + 2*n*k + 2*k + 3 + (2**h)*4*h + (2**h)*2*h*k) + (2**h)*3*h + ind] = -1
        k_q.append(k_2)
        u_q.append(u)
        v_q.append(v)
        w_q.append(w)
   
    
    print(f"abc length: {np.concatenate([a_upper, a_lower]).shape}")
    print(f"k length: {np.array(k_q).shape}")  
    
    np.savetxt("../../input/aL.txt", a.astype(int), delimiter=' ', newline=" ", fmt="%0d")
    np.savetxt("../../input/aO.txt", c.astype(int), delimiter=' ', newline=" ", fmt="%0d")
    np.savetxt("../../input/aR.txt", b.astype(int), delimiter=' ', newline=" ", fmt="%0d")
    np.savetxt("../../input/cs.txt", np.array(k_q).astype(int), delimiter=' ', newline=" ", fmt="%0d")
    np.savetxt("../../input/wL.txt", np.array(u_q).reshape((-1)).astype(int), delimiter=' ', newline=" ", fmt="%0d")
    np.savetxt("../../input/wO.txt", np.array(w_q).reshape((-1)).astype(int), delimiter=' ', newline=" ", fmt="%0d")
    np.savetxt("../../input/wR.txt", np.array(v_q).reshape((-1)).astype(int), delimiter=' ', newline=" ", fmt="%0d")

    return (a, 
            b, 
            c, 
            np.array(k_q), np.array(u_q), np.array(v_q), np.array(w_q))
            
def fixed_length_decomposition(n, k):
    if abs(n) >= 2 ** k:
        raise ValueError("Absolute value of n is too large to be represented with k bits.")

    result = [0] * k  
    abs_n = abs(n)
    for i in range(k - 1, -1, -1):
        coefficient = abs_n // (2 ** i)
        abs_n -= coefficient * (2 ** i)
        result[i] = coefficient * (2 ** i) if n >= 0 else -coefficient * (2 ** i)
    return result
    
def fixed_length_binary_decomposition(n, k):
    if abs(n) >= 2 ** k:
        raise ValueError("Absolute value of n is too large to be represented with k bits.")

    result = [0] * k  
    abs_n = abs(n)
    for i in range(k - 1, -1, -1):
        coefficient = abs_n // (2 ** i)
        abs_n -= coefficient * (2 ** i)
        result[i] = coefficient if n >= 0 else -coefficient
    return result


order = bn128_curve.curve_order
field = PrimeField(order)
size = 4
srsX = 12
srsAlpha = 10


def rPoly(aL, aR, aO, n):
    list_of_coeff = np.concatenate([np.flip(aO), np.flip(aR), aL])
    list_of_power = np.concatenate([np.arange(-2*n, 0), np.arange(1,n+1)])
    list_of_bi_coeff = []
    for i in range(len(list_of_coeff)):
        dummy_dict = {}
        dummy_dict[list_of_power[i]] = list_of_coeff[i]
        list_of_bi_coeff.append(field.sparse(dummy_dict))
    return field.sparse_bivar(dict(zip(list_of_power, list_of_bi_coeff)))


def sPoly(u,v,w, n, q):
    uiYs = []
    viYs = []
    wiYs = []
    for i in range(n):
        uiYs.insert(0,field.sparse(dict(zip(np.arange(n+1, n+q+1), u[:,i]))))
        viYs.append(field.sparse(dict(zip(np.arange(n+1, n+q+1), v[:,i]))))
        wiPart1 = field.sparse(dict(zip(np.arange(n+1, n+q+1), w[:,i])))
        wiPart2 = field.sparse(dict(zip([-i-1, i+1], [-1, -1])))
        wiYs.append(field.add_polys(wiPart1[0], wiPart2[0], wiPart1[1], wiPart2[1]))
    return field.sparse_bivar(dict(zip(np.concatenate([np.arange(-n, 0), np.arange(1,2*n+1)]), np.concatenate([uiYs, viYs, wiYs], dtype=object))))

def kPoly(k, n, q):
    return [[field.mul_by_const(field.sparse(dict(zip(np.arange(n+1, n+q+1), k))), -1)], 0]

class KZGCommitment():
    def __init__(self, n, srsX, srsAlpha, field):
        self.G1 = bn128_curve.G1
        self.G2 = bn128_curve.G2
        self.srsD = n * 6
        self.gNegativeX = [bn128_curve.multiply(bn128_curve.G1, field.exp(srsX, -i)) for i in range(1,self.srsD)]
        self.gPositiveX = [bn128_curve.multiply(bn128_curve.G1, field.exp(srsX, i)) for i in range(0,self.srsD)]
        self.hPositiveX = [bn128_curve.multiply(bn128_curve.G2, field.exp(srsX, i)) for i in range(0,2)]
        self.field = field

    def commit(self, p, init_order):
        c = None
        for i in range(len(p)):
            if init_order + i < 0:
                if c is None:
                    c = bn128_curve.multiply(self.gNegativeX[abs(init_order + i)-1], p[i])
                else:
                    c = bn128_curve.add(bn128_curve.multiply(self.gNegativeX[abs(init_order + i)-1], p[i]), c)
            else:
                if c is None:
                    c = bn128_curve.multiply(self.gPositiveX[init_order + i], p[i])
                else:
                    c = bn128_curve.add(bn128_curve.multiply(self.gPositiveX[init_order + i], p[i]), c)
        return c
    
    def openC(self, c, z, p, init_order):
        fz = self.field.eval_poly_at(p, z, init_order)
        dummy_dict = {}
        dummy_dict[0] = fz
        dummy_poly = self.field.mul_by_const(self.field.sparse(dummy_dict), -1)
        numerator = self.field.add_polys(p, dummy_poly[0], init_order, dummy_poly[1])
        dummy_dict = {}
        dummy_dict[0] = -z
        dummy_dict[1] = 1
        denominator = self.field.sparse(dummy_dict)
        qx = self.field.div_polys(numerator[0], denominator[0], numerator[1], denominator[1])
        return self.commit(qx[0], qx[1]), fz
    
    def verify(self, c, z, fz, w):
        leftleft = bn128_curve.add(c, bn128_curve.multiply(self.G1, self.field.neg(fz)))
        leftright = self.hPositiveX[0]
        rightleft = w
        rightright = bn128_curve.add(self.hPositiveX[1], bn128_curve.multiply(self.G2, self.field.neg(z)))
        
        e_left = bn128_pairing.pairing(leftright, leftleft)
        e_right = bn128_pairing.pairing(rightright, rightleft)
        return e_left == e_right
    
def load_srs(length, srsX):
    gNegativeXf = np.load('gNegativeX.npy', allow_pickle=True)
    gNegativeX = gNegativeXf[0:length, :]
    gPositiveXf = np.load('gPositiveX.npy', allow_pickle=True)
    gPositiveX = gPositiveXf[0:length, :]
    hNegativeXf = np.load('hNegativeX.npy', allow_pickle=True)
    hNegativeX = np.char.split(hNegativeXf[0:length, :].astype(str), ",") 
    hPositiveXf = np.load('hPositiveX.npy', allow_pickle=True)
    hPositiveX = np.char.split(hPositiveXf[0:length, :].astype(str), ",")
    gNegativeAlphaXf = np.load('gNegativeAlphaX.npy', allow_pickle=True)
    gNegativeAlphaX = gNegativeAlphaXf[0:length, :]
    gPositiveAlphaXf = np.load('gPositiveAlphaX.npy', allow_pickle=True)
    gPositiveAlphaX = gPositiveAlphaXf[0:length, :]
    hNegativeAlphaXf = np.load('hNegativeAlphaX.npy', allow_pickle=True)
    hNegativeAlphaX = np.char.split(hNegativeAlphaXf[0:length, :].astype(str), ",")
    hPositiveAlphaXf = np.load('hPositiveAlphaX.npy', allow_pickle=True)
    hPositiveAlphaX = np.char.split(hPositiveAlphaXf[0:length, :].astype(str), ",")
    
    for i in range(hNegativeX.shape[0]):
        hNegativeX[i][0][0] = int(hNegativeX[i][0][0].strip("()"))
        hNegativeX[i][0][1] = int(hNegativeX[i][0][1].strip("()"))
        hNegativeX[i][0] = FQ2(hNegativeX[i][0])
        hPositiveX[i][0][0] = int(hPositiveX[i][0][0].strip("()"))
        hPositiveX[i][0][1] = int(hPositiveX[i][0][1].strip("()"))
        hPositiveX[i][0] = FQ2(hPositiveX[i][0])
        hNegativeAlphaX[i][0][0] = int(hNegativeAlphaX[i][0][0].strip("()"))
        hNegativeAlphaX[i][0][1] = int(hNegativeAlphaX[i][0][1].strip("()"))
        hNegativeAlphaX[i][0] = FQ2(hNegativeAlphaX[i][0])
        hPositiveAlphaX[i][0][0] = int(hPositiveAlphaX[i][0][0].strip("()"))
        hPositiveAlphaX[i][0][1] = int(hPositiveAlphaX[i][0][1].strip("()"))
        hPositiveAlphaX[i][0] = FQ2(hPositiveAlphaX[i][0])
        
        hNegativeX[i][1][0] = int(hNegativeX[i][1][0].strip("()"))
        hNegativeX[i][1][1] = int(hNegativeX[i][1][1].strip("()"))
        hNegativeX[i][1] = FQ2(hNegativeX[i][1])
        hPositiveX[i][1][0] = int(hPositiveX[i][1][0].strip("()"))
        hPositiveX[i][1][1] = int(hPositiveX[i][1][1].strip("()"))
        hPositiveX[i][1] = FQ2(hPositiveX[i][1])
        hNegativeAlphaX[i][1][0] = int(hNegativeAlphaX[i][1][0].strip("()"))
        hNegativeAlphaX[i][1][1] = int(hNegativeAlphaX[i][1][1].strip("()"))
        hNegativeAlphaX[i][1] = FQ2(hNegativeAlphaX[i][1])
        hPositiveAlphaX[i][1][0] = int(hPositiveAlphaX[i][1][0].strip("()"))
        hPositiveAlphaX[i][1][1] = int(hPositiveAlphaX[i][1][1].strip("()"))
        hPositiveAlphaX[i][1] = FQ2(hPositiveAlphaX[i][1])
        
        hNegativeX[i] = tuple(hNegativeX[i])
        hPositiveX[i] = tuple(hPositiveX[i])
        hNegativeAlphaX[i] = tuple(hNegativeAlphaX[i])
        hPositiveAlphaX[i] = tuple(hPositiveAlphaX[i])
        
    
    return gNegativeX, gPositiveX, hPositiveX, hNegativeX, gNegativeAlphaX, gPositiveAlphaX, hPositiveAlphaX, hNegativeAlphaX


class KZGBatchCommitment():
    def __init__(self, n, srsX, srsAlpha, field):
        self.G1 = bn128_curve.G1
        self.G2 = bn128_curve.G2
        self.srsD = n * 6
        srss = load_srs(self.srsD, srsX)
        self.gNegativeX = srss[0]
        self.gPositiveX = srss[1]
        self.hPositiveX = srss[2]
        self.hNegativeX = srss[3]
        self.gNegativeAlphaX = srss[4]
        self.gPositiveAlphaX = srss[5]
        self.hPositiveAlphaX = srss[6]
        self.hNegativeAlphaX = srss[7]

        self.field = field
        self.srsX = srsX

    def commit(self, list_of_p, list_of_init_order, g_max):
        list_of_c = []
        for j, p in enumerate(list_of_p):
            c = None
            init_order = list_of_init_order[j]
            for i in range(len(p)):
                index = init_order + i + self.srsD - g_max
                if index < 0:
                    if c is None:
                        c = bn128_curve.multiply(self.gNegativeX[abs(index)-1], p[i])
                    else:
                        c = bn128_curve.add(bn128_curve.multiply(self.gNegativeX[abs(index)-1], p[i]), c)
                else:
                    if c is None:
                        c = bn128_curve.multiply(self.gPositiveX[index], p[i])
                    else:
                        c = bn128_curve.add(bn128_curve.multiply(self.gPositiveX[index], p[i]), c)
            list_of_c.append(c)
            
        return list_of_c
    
    def commita(self, list_of_p, list_of_init_order, g_max):
        list_of_c = []
        for j, p in enumerate(list_of_p):
            c = None
            init_order = list_of_init_order[j]
            for i in range(len(p)):
                index = init_order + i
                if index < 0:
                    if c is None:
                        c = bn128_curve.multiply(self.gNegativeAlphaX[abs(index)-1], p[i])
                    else:
                        c = bn128_curve.add(bn128_curve.multiply(self.gNegativeAlphaX[abs(index)-1], p[i]), c)
                elif index == 0:
                    print(f"zero term with coefficient {p[i]}")
                else:
                    if c is None:
                        c = bn128_curve.multiply(self.gPositiveAlphaX[index-1], p[i])
                    else:
                        c = bn128_curve.add(bn128_curve.multiply(self.gPositiveAlphaX[index-1], p[i]), c)
            list_of_c.append(c)
            
        return list_of_c
    
    
    def openC(self, list_of_c, list_of_z_for_p, list_of_p, list_of_init_order, g_max):
        
        # evaluate all polynomials at all points and then interpolations
        list_of_fz = []
        list_of_all_z = list(set().union(*list_of_z_for_p))
        list_of_gamma = []
        zTx = list(np.array(np.poly1d(list_of_all_z, True).c)[::-1].astype(int))
#         print(zTx)
        list_of_zSix = []
        fx = [[0],0]
        fxz = [[0],0]
        beta = 1
        rand_z = 1
        for i, p in enumerate(list_of_p):
            list_of_fz_p = []
            # compute all p(z) for z \in S_i
            for z in list_of_z_for_p[i]:
                list_of_fz_p.append(self.field.eval_poly_at(p, z, list_of_init_order[i]))
                
            list_of_fz.append(list_of_fz_p)           
            # compute r_i
            gamma_p = self.field.lagrange(list_of_z_for_p[i], list_of_fz_p)
            # compute z_{T\Si}
            SExSi = list(set(list_of_all_z) ^ set(list_of_z_for_p[i]))
            zSix = list(np.array(np.poly1d(SExSi, True).c)[::-1].astype(int))
            # f_i(x)-r_i(x)
            sub_poly = self.field.sub_polys(p, gamma_p[0], list_of_init_order[i], gamma_p[1])
            # z_{T\Si}*[f_i(x)-r_i(x)]
            part_poly = self.field.mul_polys(zSix, sub_poly[0], 0, sub_poly[1])
            # beta^i * z_{T\Si}*[f_i(x)-r_i(x)]
            part_poly = self.field.mul_by_const(part_poly, self.field.exp(beta, i))
            # summation
            fx = self.field.add_polys(fx[0], part_poly[0], fx[1], part_poly[1])
            # r_i(z)
            eval_gamma = [[self.field.eval_poly_at(gamma_p[0], rand_z, gamma_p[1])],0]
            list_of_gamma.append(eval_gamma[0][0])
            # compute z_{T\Si}(z)
            eval_zSix = [[self.field.eval_poly_at(zSix, rand_z, 0)],0]
            list_of_zSix.append(eval_zSix[0][0])
#             print(eval_zSix)
            # f_i(x)-r_i(z)
            sub_poly = self.field.sub_polys(p, eval_gamma[0], list_of_init_order[i], eval_gamma[1])
            # z_{T\Si}(z)*[f_i(x)-r_i(z)]
            part_poly = self.field.mul_polys(eval_zSix[0], sub_poly[0], 0, sub_poly[1])
            # beta^i * z_{T\Si}(z)*[f_i(x)-r_i(z)]
            part_poly = self.field.mul_by_const(part_poly, self.field.exp(beta, i))
            fxz = self.field.add_polys(fxz[0], part_poly[0], fxz[1], part_poly[1])
        
        # f/Z_T
        px = self.field.div_polys(fx[0], zTx, fx[1], 0)
        test = 6
        
        # Z_T(z)
        eval_zt = self.field.eval_poly_at(zTx, rand_z, 0)
        # f/Z_T * Z_T(z)
        l_second_half = self.field.mul_by_const(px, eval_zt)
        # L = beta^i * z_{T\Si}(z)*[f_i(x)-r_i(z)] - f/Z_T * Z_T(z)
        lx = self.field.sub_polys(fxz[0], l_second_half[0], fxz[1], l_second_half[1])
#         print(self.field.eval_poly_at(lx[0], 13295, lx[1])/(13295-rand_z))
        # L / (x-z)
        lx = self.field.div_polys(lx[0], [-rand_z, 1], lx[1], 0)
        
        w = self.commit([px[0]], [px[1]], g_max)
        wDash = self.commit([lx[0]], [lx[1]], g_max)  
        
        return w[0], wDash[0], rand_z, beta, list_of_fz, list_of_gamma, list_of_zSix, eval_zt, g_max
        
    def verify(self, list_of_c, w, wDash, rand_z, beta, list_of_fz, list_of_gamma, list_of_zSix, eval_zt, g_max):
        
        cm_poly = None
        ga_poly = 0
        for i, c in enumerate(list_of_c):
            gamma_p = list_of_gamma[i]
            zSix = list_of_zSix[i]
            if cm_poly == None:
                cm_poly = bn128_curve.multiply(c, self.field.mul(zSix, self.field.exp(beta, i)))
            else:
                cm_poly = bn128_curve.add(bn128_curve.multiply(c, self.field.mul(zSix, self.field.exp(beta, i))), cm_poly)
            ga_poly = self.field.add(self.field.mul(gamma_p, self.field.mul(zSix, self.field.exp(beta, i))), ga_poly)

        w_poly = bn128_curve.multiply(w, eval_zt)

        non_a_part = bn128_curve.neg(bn128_curve.add(w_poly, bn128_curve.multiply(bn128_curve.G1, ga_poly)))   
        non_a_part_w_a = bn128_curve.multiply(bn128_curve.add(non_a_part, bn128_curve.multiply(wDash, rand_z)), srsAlpha)
        left1right = self.hPositiveX[self.srsD - g_max].astype(FQ2)
        left2left = bn128_curve.add(cm_poly, non_a_part_w_a)
        rightleft = wDash
        rightright = self.hPositiveAlphaX[1 + self.srsD - g_max].astype(FQ2)
        e_left = bn128_pairing.pairing(left1right, left2left)
        e_right = bn128_pairing.pairing(rightright, rightleft)
        return e_left == e_right
    

def parse_number(str_in):
    i = 0
    while str_in[i].isdigit() or str_in[i] == "-":
        i += 1

    return str_in[:i], i

def parse_term(str_in):
    str_in = str_in.strip()
    order = int(str_in.split(",")[0][1:])
    coeff = int(str_in.split(")")[0].split("P")[1])
    if len(str_in.split(")")[1].strip()) == 0:
        pow = 0
    else:
        if len(str_in.split(")")[1].strip().split("^")) == 1:
            pow = 1
        else:
            pow = int(str_in.split(")")[1].split("^")[1])

    return order, coeff, pow

def make_poly(coeffs):
    init_order = coeffs[-1][-1]
    coeff_dict = {}
    for coef in coeffs:
        coeff_dict[coef[0]+init_order] = coef[1]

    return field.sparse(coeff_dict)

def make_poly_bivar(coeffs):
    init_order = coeffs[-1][-1]
    coeff_dict = {}
    for coef in coeffs:
        coeff_dict[coef[0]+init_order] = coef[1]

    return field.sparse_bivar(coeff_dict)
    
def read_poly(string_in):
    new_poly = [[], 0]
    tokens = string_in.split("+")
    list_of_tokens = []
    main_poly = []
    for token in tokens:
        if token.count("(") == 2 and token.count(")") == 1:
            list_of_tokens = []
            tos = token.split("(")
            order = int(tos[1][:-1])
            coeff = parse_term("("+tos[2])
            list_of_tokens.append(order)
            list_of_tokens.append(coeff)
        
        elif token.count("(") == 1 and token.count(")") == 1:
            list_of_tokens.append(parse_term(token))

        elif token.count("(") == 1 and token.count(")") == 2:
            tos = token.split(")")
            coeff = parse_term(tos[0]+")"+tos[1])
            list_of_tokens.append(coeff)
            if len(tos) == 2:
                pow = 0
            elif len(tos[2].split("^")) == 1:
                pow = 1
            else:
                pow = int(tos[2].split("^")[1])

            list_of_tokens.append(pow)

            main_poly.append((list_of_tokens[0], make_poly(list_of_tokens[1:-1]), list_of_tokens[-1]))

        elif token.count("(") == 2 and token.count(")") == 2:
            tos = token.split("(")
            order = int(tos[1][:-1])
            toss = tos[2].split(")")
            coeff = parse_term("("+toss[0]+")"+toss[1])

            if len(toss) == 2:
                pow = 0
            elif len(toss[2].split("^")) == 1:
                pow = 1
            else:
                pow = int(toss[2].split("^")[1])

            main_poly.append((order, make_poly([coeff]), pow))

    return make_poly_bivar(main_poly)


def sonic_experiment(size, aL, aR, aO, k, u, v, w, n, q, save=False, load=False):
    
    cmScheme = KZGCommitment(n, srsX, srsAlpha, field)
    
    if save:
        with open(f"../../output/polys_dt.txt", "r") as f:
            polys = f.readline()
            polys = polys.split("=")[1:]
            last = polys[-1]
            polys = [x.rsplit(',', 1)[0].strip() for x in polys[:-1]]
            polys.append(last[:-1])

        sXY = read_poly(polys[0])
        rXY = read_poly(polys[1])
        tXY = read_poly(polys[-2])

        neg_kXY = kPoly(k, n, q)
        with open(f"sXY_{size}.txt", 'wb') as f:
            pickle.dump(sXY, f)

        with open(f"rXY_{size}.txt", 'wb') as f:
            pickle.dump(rXY, f)

        with open(f"tXY_{size}.txt", 'wb') as f:
            pickle.dump(tXY, f)
                
        with open(f"neg_kXY_{size}.txt", 'wb') as f:
            pickle.dump(neg_kXY, f)

        

    if load:
        with open(f"sXY_{size}.txt", 'rb') as f:
            sXY = pickle.load(f)

        with open(f"rXY_{size}.txt", 'rb') as f:
            rXY = pickle.dump(f)

        with open(f"tXY_{size}.txt", 'rb') as f:
            tXY = pickle.dump(f)

        with open(f"neg_kXY_{size}.txt", 'rb') as f:
            neg_kXY = pickle.dump(f)

    st = time.process_time()
    
    
    if save:
        rX1 = field.eval_poly_Y(rXY[0], rXY[1], 1)
        rX1 = field.dimension_change(rX1[0], rX1[1])

        r1Raw = [rX1[0][:len(rX1)-n//2], rX1[1]]
        r1Local = field.sub_polys(rX1[0], r1Raw[0], rX1[1], r1Raw[1])
        print(f"N: {n}")
        r1Local = field.mul_polys(r1Local[0], [1], r1Local[1], n)

        with open(f"rX1_{size}.txt", 'wb') as f:
            pickle.dump(rX1, f)

        with open(f"r1Raw_{size}.txt", 'wb') as f:
            pickle.dump(r1Raw, f)

        with open(f"r1Local_{size}.txt", 'wb') as f:
            pickle.dump(r1Local, f)


    if load:
        with open(f"rX1_{size}.txt", 'rb') as f:
            rX1 = pickle.dump(f)

        with open(f"r1Raw_{size}.txt", 'b') as f:
            r1Raw = pickle.dump(f)

        with open(f"r1Local_{size}.txt", 'rb') as f:
            r1Local = pickle.dump(f)
        
    et = time.process_time()
    res = et - st
    print('CPU Execution time-poly construction:', res, 'seconds')


    st = time.process_time()

    g_max = 6 * n

    s1Y = field.eval_poly_X(sXY[0], sXY[1], 1)
    kY = field.mul_by_const(neg_kXY[0][0], -1)
    kY

    commitS1Y = cmScheme.commit(s1Y[0], s1Y[1])
    commitK = cmScheme.commit(kY[0], kY[1])

    # commit R
    commitR = cmScheme.commit(rX1[0], rX1[1])

    commitR

    commitRRaw = cmScheme.commit(r1Raw[0], r1Raw[1])
    commitRLocal = cmScheme.commit(r1Local[0], r1Local[1])


    # y
    y = 3
    #commit T
    tXy = field.eval_poly_Y(tXY[0], tXY[1], y)
    tXy = field.dimension_change(tXy[0], tXy[1])
    commitT = cmScheme.commit(tXy[0], tXy[1])
    

    #commit Sx
    sXy = field.eval_poly_Y(sXY[0], sXY[1], y)
    sXy = field.dimension_change(sXy[0], sXy[1])
    commitSx = cmScheme.commit(sXy[0], sXy[1])
    

    # z
    z = 2

    #open all together

    list_of_c = [commitR, commitRRaw, commitRLocal, commitT, commitS1Y, commitSx, commitK] 
    list_of_z_for_p = [[z, z*y], [z], [z], [z], [y], [z, 1], [y]]
    list_of_p = [rX1[0], r1Raw[0], r1Local[0], tXy[0], s1Y[0], sXy[0], kY[0]]
    list_of_init_order = [rX1[1], r1Raw[1], r1Local[1], tXy[1], s1Y[1], sXy[1], kY[1]]

    opens = []
    for i, com in enumerate(list_of_c):
        z_for_p = list_of_z_for_p[i]
        for z in z_for_p:
            opens.append(cmScheme.openC(com, z, list_of_p[i], 
                                        list_of_init_order[i]))

    et = time.process_time()
    
    res = et - st
    print('CPU Execution time - proof generation:', res, 'seconds')

    pi_D = opens[2][0]
    pi_R_tilde = opens[3][0]
    pi_R1 = opens[0][0]
    pi_R2 = opens[1][0]
    pi_T = opens[4][0]
    pi_K = opens[8][0]
    pi_S_x1 = opens[6][0]
    pi_S_x2 = opens[7][0]
    pi_S_y = opens[5][0]
    r_1 = opens[0][1]
    r_tilde = opens[3][1]
    t = opens[4][1]
    k = opens[8][1]
    s_tilde = opens[6][1]
    r_2 = opens[1][1]
    s_1_tilde = opens[7][1]
    s_2_tilde = opens[5][1]
    # poly commitments; Fs
    D = commitRRaw
    R_tilde = commitR
    R = commitRLocal
    T = commitT
    K = commitK
    S_x = commitSx
    S_y = commitS1Y

    d = opens[2][1]
    
    
    # verify
        
    print(f"\
        G1Point pi_D = {pi_D} \n \
        G1Point pi_R_tilde = {pi_R_tilde} \n \
        G1Point pi_R1 = {pi_R1} \n \
        G1Point pi_R2 = {pi_R2} \n \
        G1Point pi_T = {pi_T} \n \
        G1Point pi_K = {pi_K} \n \
        G1Point pi_S_x1 = {pi_S_x1} \n \
        G1Point pi_S_x2 = {pi_S_x2} \n \
        G1Point pi_S_y = {pi_S_y} \n \
        uint256 r_1 = {r_1} \n \
        uint256 r_tilde = {r_tilde} \n \
        uint256 t = {t} \n \
        uint256 k = {k} \n \
        uint256 s_tilde = {s_tilde} \n \
        uint256 r_2 = {r_2} \n \
        uint256 s_1_tilde = {s_1_tilde} \n \
        uint256 s_2_tilde = {s_2_tilde} \n \
        G1Point D = {D} \n \
        G1Point R_tilde = {R_tilde} \n \
        G1Point R = {R} \n \
        G1Point T = {T} \n \
        G1Point K = {K} \n \
        G1Point S_x = {S_x} \n \
        G1Point S_y = {S_y} \n \
        uint256 d = {d}")


    cmVerify = True
    j = 0
    for i, com in enumerate(list_of_c):
       z_for_p = list_of_z_for_p[i]
       for z in z_for_p:
           cmVerify = cmVerify and cmScheme.verify(com, 
                                        z, 
                                        opens[j][1], 
                                       opens[j][0])
           j += 1
           
    return cmVerify and t == field.sub(field.mul(r_1, field.add(r_2, s_tilde)), k) and s_1_tilde == s_2_tilde

# Change the data for DTs with different heights:
def setup_data():
    #232
    aL, aR, aO, k, u, v, w = decision_tree(np.array([1, 1, 1, 1, 1, 5000, 6982, 548, 421, 3, 4]), np.array([10]), np.array([2]), np.array([0, 1]), 5, 1, 1)
    assignment = [aL, aR, aO]
    circuit = [u,v,w,k]

    assert (aL @ u.T + aR @ v.T + aO @ w.T == k).all()

    n = aL.shape[0]
    q = k.shape[0]
    
    print("constraints generation completed")

    return aL, aR, aO, k, u, v, w, n, q

if __name__ == "__main__":
    aL, aR, aO, k, u, v, w, n, q = setup_data()

    result = sonic_experiment(4, aL, aR, aO, k, u, v, w, n, q, True)
    print(result)