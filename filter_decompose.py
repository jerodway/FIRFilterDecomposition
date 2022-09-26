import scipy.signal
import scipy.special
import numpy as np
import cmath
import math

import matplotlib.pyplot as plt

DBL_EPSILON = 2.2204460492503131e-16
DBL_MAX = 1.7976931348623157e308
MULLER_ITTERMAX = 150
MULLER_FVALUE = 1e36
MULLER_NOISESTART = DBL_EPSILON * 1e2
MULLER_NOISEMAX = 5
MULLER_MAXDIST = 1e3
MULLER_KITERMAX = 1e3
MULLER_BOUND4 = math.sqrt(DBL_MAX) / 1e4
MULLER_BOUND6 = math.log10(MULLER_BOUND4) - 4
CONVERGENCE = 100
EPSILON_TO_USE = 5 * DBL_EPSILON


class Poly:
    def __init__(self, fil):
        self.num_coefficients = len(fil)
        self.current_degree = self.num_coefficients - 1
        self.original_degree = self.current_degree
        self.degree_reduced = 0

        self.original_coeffs = list(fil)
        self.num_coefficients = len(self.original_coeffs)
        self.num_distinct = len(set(self.original_coeffs))
        self.coeffs = [complex(x, 0.0) for x in list(fil)]

        self.roots = [complex(0.0, 0.0)] * self.original_degree
        self.polynomial_reduced = [
            complex(
                0.0,
                0.0,
            )
        ] * (self.num_coefficients)
        self.psave = self.polynomial_reduced  # .copy()

        # initialization for muller from l764
        self.x0 = complex(0.0, 1.0)
        self.x1 = complex(0.0, -1.0)
        self.x2 = complex(1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0))
        self.h1 = self.x1 - self.x0
        self.h2 = self.x2 - self.x1
        self.q2 = self.h2 / self.h1

        self.f0 = None
        self.f1 = None
        self.f2 = None

        self.xb = self.x2
        self.epsilon = 5 * DBL_EPSILON
        self.iter = 0

        self.flag = (
            True  # not sure why this is necessary, but flag for real coefficients?
        )

        self.gain = None

    def get_roots(self):
        # self.current_degree -= 1 # reduce for indexing from zero
        self.degree_reduced = self.current_degree

        diff = 0  # number of roots at zero

        self.poly_check()

        diff = self.current_degree - self.degree_reduced

        # pointer should change?
        self.current_degree = self.degree_reduced

        # check for linear or quadratic, which can be solved directly l438

        # convert to monic polynomial
        self.monic()
        self.psave = self.coeffs.copy()

        # prepare for input of Muller
        for i in range(
            0, self.current_degree + 1
        ):  # from 0 to current_degree + 1 because these are the coefficients
            self.polynomial_reduced[i] = self.coeffs[i]

        max_error = 0.0

        while self.degree_reduced > 2:
            new_solution = self.muller()

            self.roots[self.degree_reduced - 1], new_error = self.newton(new_solution)

            if new_error > max_error:
                max_error = new_error

            reduction = self.polynomial_deflate(self.flag)
            # self.pred = self.pred + reduction
            for i in range(0, len(self.polynomial_reduced)):
                self.psave[i] = self.polynomial_reduced[i]

            self.polynomial_reduced = self.polynomial_reduced[reduction:]
            self.degree_reduced = self.degree_reduced - reduction
        self.lin_or_quadratic()
        # extras if degree_reduced == 2
        if self.degree_reduced == 2:
            if abs(self.roots[1]) <= 1:
                self.roots[1], new_error = self.newton(self.roots[1])
                if new_error > max_error:
                    max_error = new_error
            if abs(self.roots[0]) <= 1:
                self.roots[0], new_error = self.newton(self.roots[0])

        if max_error > 9e-5:
            print("ROOT FINDING FAILED")


    def poly_check(self):
        i = -1

        if self.current_degree < 0:
            return 1

        for j in range(0, self.current_degree + 1):
            if abs(self.coeffs[j]) != 0.0:
                i = j

        # if polynomial is null
        if i == -1:
            return 2

        # if polynomial is all zeros
        if i == 0:
            return 3

        # get the new order of the polynomial
        self.current_degree = i

        # check for number of zeros in input polynomial
        zero_count = 0
        # not_found = True
        for i in range(0, self.current_degree):
            if abs(self.coeffs[i]) == 0.0:
                zero_count += 1
            else:
                break
        if (
            zero_count == 0
        ):  # case with no zeros just floating around, degree isn't reduced at all
            self.degree_reduced = self.current_degree
            return 0
        else:
            for j in range(0, zero_count + 1):
                self.roots[self.current_degree - j - 1] = complex(0.0, 0.0)
        self.degree_reduced = self.current_degree - zero_count

    # computes a monic polynomial from the original
    def monic(self):
        factor = 1.0 / abs(self.coeffs[self.current_degree])
        if factor != 1.0:
            self.coeffs = [x * factor for x in self.coeffs]
            # for i in range(0, self.current_degree+1):
            #    self.coeffs[i] *= factor

    def muller(self):
        f2absq = MULLER_FVALUE
        f2absqb = MULLER_FVALUE
        f1absq = None

        second_iter = 0  # second iteration for when root is too bad?
        noise = 0  # noise counter
        rootd = False  # ?

        # apparently we have to do this every time
        self.x0 = complex(0.0, 1.0)
        self.x1 = complex(0.0, -1.0)
        self.x2 = complex(1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0))
        self.h1 = self.x1 - self.x0
        self.h2 = self.x2 - self.x1
        self.q2 = self.h2 / self.h1
        self.xb = self.x2
        self.iter = 0

        # use Horner's Method, get f0=P(x0), f1=P(x1), f2=P(x2)
        self.f0 = self.muller_f_value(self.degree_reduced, self.x0)
        self.f1 = self.muller_f_value(self.degree_reduced, self.x1)
        self.f2 = self.muller_f_value(self.degree_reduced, self.x2)

        while second_iter <= 2:
            while (
                self.iter < MULLER_ITTERMAX
                and not rootd
                and noise <= MULLER_NOISEMAX
            ):
                self.muller_root_of_parabola()

                self.x0 = self.x1
                self.x1 = self.x2
                h2_abs = abs(self.h2)

                h2_abs = self.iteration_equation(h2_abs)

                self.f0 = self.f1
                self.f1 = self.f2
                f1absq = f2absq

                f2absq = self.compute_function(f1absq, f2absq, EPSILON_TO_USE)

                # check if the new x2 is good enough , these two checks are necessary
                self.xb, f2absqb, noise, rootd = self.check_x_value(
                    self.xb, f1absq, f2absq, f2absqb, rootd, EPSILON_TO_USE, noise
                )
                if (
                    abs((abs(self.xb) - abs(self.x2)) / abs(self.xb))
                    < MULLER_NOISESTART
                ):
                    noise += 1
            second_iter += 1
            second_iter, rootd, noise = self.root_check(
                f2absqb, second_iter, rootd, noise, self.xb
            )
        return self.xb

    def muller_f_value(self, n: int, x0: complex):
        f = self.polynomial_reduced[n]
        index = list(range(0, n))
        index.reverse()
        for i in index:
            temp = f * x0
            f = temp + self.polynomial_reduced[i]
        return f

    # root of muller's parabola l778
    def muller_root_of_parabola(self):
        # A2 = q2(f2 - (1 + q2)* f1 + f0q2)
        # B2 = q2[q2(f0 - f1) + 2(f2 - f1)] + (f2 - f1)
        # C2 = (1 + q2) * f[2]

        a_2 = self.q2 * (self.f2 - (1.0 + self.q2) * self.f1 + self.f0 * self.q2)
        b_2 = self.q2 * (self.q2 * (self.f0 - self.f1) + 2.0 * (self.f2 - self.f1)) + (
            self.f2 - self.f1
        )
        c_2 = (1.0 + self.q2) * self.f2

        # discr = B2^2 - 4A2C2
        discriminant = b_2 ** 2 - 4.0 * a_2 * c_2

        # denominators of q2
        n_1 = b_2 - cmath.sqrt(discriminant)
        n_2 = b_2 + cmath.sqrt(discriminant)

        # choose denominator with largest modulus
        if abs(n_1) > abs(n_2) and abs(n_1) > DBL_EPSILON:
            self.q2 = -2.0 * c_2 / n_1
        elif abs(n_2) > DBL_EPSILON:
            self.q2 = -2.0 * c_2 / n_2
        else:
            self.q2 = complex(math.cos(self.iter), math.sin(self.iter))

    # muller's iteration equation l801
    # main iteration equation: x2 = h2*q2 + x2
    def iteration_equation(self, h2_abs):
        self.h2 = self.h2 * self.q2
        new_h2_abs = abs(self.h2)
        if new_h2_abs > (h2_abs * MULLER_MAXDIST):
            temp = MULLER_MAXDIST / new_h2_abs
            self.h2 = self.h2 * temp
            self.q2 = self.q2 * temp

        self.x2 = self.x2 + self.h2
        return new_h2_abs

    # muller's compute function l886
    # compute P(x2) and make some checks
    def compute_function(self, f1absq, f2absq, epsilon):
        overflow = True
        while overflow:
            overflow = False

            self.suppress_overflow()

            # calculate new value => result in f2
            self.f2 = self.muller_f_value(self.degree_reduced, self.x2)

            f2absq = self.too_big_function_values(f2absq)

            self.iter += 1

            # Muller's modification to improve convergence
            overflow = self.convergence_check(overflow, f1absq, f2absq, epsilon)

        return f2absq

    # suppress overflow? l842
    def suppress_overflow(self):
        loop = True
        k_iter = 0
        while loop:
            loop = False
            abs_x2 = abs(self.x2)

            if (
                abs_x2 > 1.0
                and abs(self.degree_reduced * math.log10(abs_x2)) > MULLER_BOUND6
            ):
                k_iter += 1
                if k_iter < MULLER_KITERMAX:
                    # halve distance between new and old x2
                    self.h2 = 0.5 * self.h2
                    self.q2 = 0.5 * self.q2
                    self.x2 = self.x2 = self.h2
                    loop = True
                else:
                    k_iter = 0

    def too_big_function_values(self, f2absq):
        if abs(self.f2.real) + abs(self.f2.imag) > MULLER_BOUND4:
            f2absq = abs(self.f2.real) + abs(self.f2.imag)
        else:
            f2absq = self.f2.real ** 2 + self.f2.imag ** 2
        return f2absq

    # Muller's modification to improve convergence l872
    def convergence_check(self, overflow, f1absq, f2absq, epsilon):

        if (
            f2absq > (CONVERGENCE * f1absq)
            and abs(self.q2) > epsilon
            and self.iter < MULLER_ITTERMAX
        ):
            self.q2 = 0.5 * self.q2
            self.h2 = 0.5 * self.h2
            self.x2 = self.x2 - self.h2
            overflow = True

        return overflow

    def check_x_value(
        self, xb: complex, f1absq, f2absq, f2absqb, rootd, epsilon, noise
    ):
        BOUND1 = 1.01
        BOUND2 = 0.99
        BOUND3 = 0.01

        if f2absq <= (BOUND1 * f1absq) and f2absq >= (BOUND2 * f1absq):
            if abs(self.h2) < BOUND3:
                self.q2 = self.q2 * 2.0
                self.h2 = self.h2 * 2.0
            else:
                self.q2 = complex(math.cos(self.iter), math.sin(self.iter))
                self.h2 = self.h2 * self.q2
        elif f2absq < f2absqb:
            f2absqb = f2absq
            xb = self.x2
            noise = 0
            if (
                math.sqrt(f2absq) < epsilon
                and abs((self.x2 - self.x1) / self.x2) < epsilon
            ):
                rootd = True

        return xb, f2absqb, noise, rootd

    def root_check(self, f2absqb, second_iter, rootd, noise, xb):
        BOUND7 = 1e-5
        df = None
        if second_iter == 1 and f2absqb > 0:
            df = self.muller_f_value_2()
            if (abs(self.f2) / (abs(df) * abs(xb))) > BOUND7:
                self.x0 = complex(
                    -1.0 / math.sqrt(2), 1.0 / math.sqrt(2)
                )  # start second iteration with new initial estimates
                self.x1 = complex(1.0 / math.sqrt(2), -1.0 / math.sqrt(2))
                self.x2 = complex(-1.0 / math.sqrt(2), -1.0 / math.sqrt(2))
                self.f0 = self.muller_f_value(self.degree_reduced, self.x0)
                self.f1 = self.muller_f_value(self.degree_reduced, self.x1)
                self.f2 = self.muller_f_value(self.degree_reduced, self.x2)
                self.iter = 0
                second_iter += 1
                rootd = False
                noise = 0

        return second_iter, rootd, noise

    def muller_f_value_2(self):
        # this seems to only ever be called for f2 and xb
        self.f2 = self.psave[self.degree_reduced]
        df = complex(0.0, 0.0)

        index = list(range(0, self.degree_reduced))
        index.reverse()
        for i in index:
            df = df * self.xb + self.f2
            self.f2 = self.f2 * self.xb + self.psave[i]

        return df

    # deflate the polynomial
    def polynomial_deflate(self, flag):

        x0 = self.roots[self.degree_reduced - 1]
        if x0.imag != 0.0:  # x0 is complex
            flag = False

        if not flag:  # if x0 is complex, then there are two roots
            a = 2 * x0.real
            b = -(x0.real ** 2 + x0.imag ** 2)
            self.roots[self.degree_reduced - 2] = x0.conjugate()
            self.horncd(a, b)
            return 2
        else:
            self.hornc(x0, flag)
            return 1

    # Horner method to deflate two roots l595
    def horncd(self, a, b):
        self.polynomial_reduced[self.degree_reduced - 1] = complex(
            self.polynomial_reduced[self.degree_reduced - 1].real
            + self.polynomial_reduced[self.degree_reduced].real * a,
            self.polynomial_reduced[self.degree_reduced - 1].imag,
        )
        index = list(range(2, self.degree_reduced - 1))
        index.reverse()
        for i in index:
            self.polynomial_reduced[i] = complex(
                self.polynomial_reduced[i].real
                + (a * self.polynomial_reduced[i + 1].real + b * self.polynomial_reduced[i + 2].real),
                self.polynomial_reduced[i].imag,
            )

    # Horner method to deflate one root l580
    def hornc(self, x0, flag):
        index = list(range(1, self.degree_reduced))
        index.reverse()

        if flag:  # real coefficients
            for i in index:
                self.polynomial_reduced[i] = complex(
                    self.polynomial_reduced[i].real + (x0.real * self.polynomial_reduced[i + 1].real),
                    self.polynomial_reduced[i].imag,
                )
        else:  # complex coefficients
            for i in index:
                temp = self.polynomial_reduced[i + 1] * x0
                self.polynomial_reduced[i] = temp + self.polynomial_reduced[i]

    # newtons method
    def newton(self, new_solution):
        new_error = None
        polished_root = None

        ITERMAX_1 = 20
        BOUND = math.sqrt(DBL_EPSILON)
        NOISEMAX = 5
        FACTOR = 5
        FVALUE = 1e36

        fabsmin = FVALUE
        noise = 0

        x0 = new_solution
        xmin = new_solution
        dx = complex(1.0, 0.0)
        new_error = abs(dx)
        # self.iter = 0
        for i in range(0, ITERMAX_1):
            f, df = self.f_value_1(x0)
            if abs(f) < fabsmin:
                xmin = x0
                fabsmin = abs(f)
                noise = 0
            if abs(df) != 0.:
                dxh = f / df
                if abs(dxh) < new_error * FACTOR:
                    dx = dxh
                    new_error = abs(dx)
            if abs(xmin) != 0.0:
                if new_error / abs(xmin) < DBL_EPSILON or noise == NOISEMAX:
                    if noise == NOISEMAX:
                        print("cutting out on noisemax")
                    if abs(xmin.imag) < BOUND and self.flag:
                        xmin = complex(
                            xmin.real, 0.0
                        )  # if the imaginary part is super small, just zero it out

                    new_error = new_error / abs(xmin)
                    return xmin, new_error

            x0 = x0 - dx
            noise = noise + 1

        # we iterated through and didn't find a root that met our bounds, so return what we got anyway
        if abs(xmin.imag) < BOUND and self.flag:
            xmin = complex(xmin.real, 0.0)
        if abs(xmin) != 0.0:
            new_error = new_error / abs(xmin)
        return xmin, new_error

    # calculate roots of quadratic
    def quadratic(self):
        discriminant = self.polynomial_reduced[1] * self.polynomial_reduced[1] - 4.0 * self.polynomial_reduced[2] * self.polynomial_reduced[0]
        z1 = -self.polynomial_reduced[1] + cmath.sqrt(discriminant)
        z2 = -self.polynomial_reduced[1] - cmath.sqrt(discriminant)
        n = 2.0 * self.polynomial_reduced[2]
        self.roots[0] = z1 / n
        self.roots[1] = z2 / n

    def lin_or_quadratic(self):
        if self.degree_reduced == 1:
            self.roots[0] = -self.polynomial_reduced[0] / self.polynomial_reduced[1]
        elif self.degree_reduced == 2:
            self.quadratic()

    # f value calculation in newton's method l
    def f_value_1(self, x0):

        f = self.coeffs[self.current_degree]
        df = complex(0.0, 0.0)
        index = list(range(0, self.current_degree))
        index.reverse()
        for i in index:
            df = df * x0 + f
            f = f * x0 + self.coeffs[i]

        return f, df

    def build_filters(self):
        inside_unit_circle = []
        outside_unit_circle = []
        on_unit_circle = []
        on_axis = []
        on_one = []

        h_0 = [1.] * self.original_degree
        h_1 = [1.] * self.original_degree
        h_2 = [1.] * self.original_degree
        h_3 = [0.] * self.original_degree
        h_4 = [0.] * self.original_degree

        temp_0 = [0.] * self.original_degree
        temp_1 = [0.] * self.original_degree
        temp_2 = [0.] * self.original_degree
        temp_3 = [0.] * self.original_degree
        temp_4 = [0.] * self.original_degree

        h_m = complex(1., 1.)

        count = 0  # filter count, I think
        k = 0  # I really have no idea
        types = []

        for root in self.roots:

            if abs(root) < (1 - 1e-8) and root.imag > 0:
                inside_unit_circle.append(root)

            if (abs(root) - 1) > 1e-8 and root.imag > 0:
                outside_unit_circle.append(root)

            if 1 - 1e-8 < abs(root) < 1 + 1e-8 and abs(root.imag) > 0:
                on_unit_circle.append(root)

            if abs(root.imag) == 0. and (abs(root.real) > 1 + 1e-8 or abs(root.real) < 1 - 1e-8):
                on_axis.append(root)

            if abs(root.imag) == 0. and (abs(root.real) < 1 + 1e-8 and abs(root.real) > 1 - 1e-8):
                on_one.append(root)

        num_inside_unit_circle = len(inside_unit_circle)
        num_outside_unit_circle = len(outside_unit_circle)
        num_on_unit_circle = len(on_unit_circle)
        num_on_axis = len(on_axis)
        num_on_one = len(on_one)

        for i in range(0, num_on_one):
            if on_one[i].real - 1 > -1e-8 and on_one[i].real - 1 < 1e-8:
                h_1[i] = -1.
            temp_0[count] = h_0[i]
            temp_1[count] = h_1[i]
            temp_2[count] = 0.
            temp_3[count] = 0.
            temp_4[count] = 0.
            count += 1
            types.append('on_one')

        for i in range(0, num_on_axis):

            if abs(on_axis[i].real) > 1:
                k += 1
            elif abs(on_axis[i].real) < 1:
                if abs(on_axis[i].real) > 0. or abs(on_axis[i].imag) > 0.:
                    h_1[i + num_on_one - k] = -(on_axis[i].real + 1 / on_axis[i].real)
                    temp_0[count] = h_0[i + num_on_one -k]
                    temp_1[count] = h_1[i + num_on_one -k]
                    temp_2[count] = h_2[i + num_on_one]
                    temp_3[count] = 0.
                    temp_4[count] = 0.
                    count += 1
                    types.append('on_axis')

        for i in range(0, num_on_unit_circle):

            if on_unit_circle[i].imag < 0:
                k += 1
            elif on_unit_circle[i].imag > 0:
                h_1[i + num_on_one + num_on_axis - k] = -2 * on_unit_circle[i].real / abs(on_unit_circle[i])
                temp_0[count] = h_0[i + num_on_one + num_on_axis - k]
                temp_1[count] = h_1[i + num_on_one + num_on_axis - k]
                temp_2[count] = h_2[i + num_on_one + num_on_axis]
                temp_3[count] = 0.
                temp_4[count] = 0.
                count += 1
                types.append('on_unit_circle')

        for i in range(0, num_inside_unit_circle):
            r = abs(inside_unit_circle[i])

            h_2[i + num_on_one + num_on_axis + num_on_unit_circle - k] = (math.pow(r, 2.) + 1 / math.pow(r, 2.)) + math.pow(2 * inside_unit_circle[i].real / abs(inside_unit_circle[i]), 2.)

            h_1[i + num_on_one + num_on_axis + num_on_unit_circle - k] = -2. * (r + 1 / r) * inside_unit_circle[i].real / abs(inside_unit_circle[i])

            temp_0[count] = h_0[i + num_on_one + num_on_axis + num_on_unit_circle - k]
            temp_1[count] = h_1[i + num_on_one + num_on_axis + num_on_unit_circle - k]
            temp_2[count] = h_2[i + num_on_one + num_on_axis + num_on_unit_circle - k]
            temp_3[count] = h_1[i + num_on_one + num_on_axis + num_on_unit_circle - k]
            temp_4[count] = 1.
            count += 1
            types.append('inside_unit_circle')

        decomposed_filters = []
        for i in range(0, count):
            loc_filt = [temp_0[i], temp_1[i], temp_2[i], temp_3[i], temp_4[i]]
            decomposed_filters.append(loc_filt)

        for i in range(0, count):
            print(f"{i}\t{decomposed_filters[i]}\t{types[i]}")


    def calculate_gain(self):
        NPLOT = 2048

        upper = 0.
        lower = 0.

        H = complex(1., 0.)

        for k in range(0, NPLOT):
            omega = math.pi * float(k) / float(NPLOT)

            w_re = math.cos(omega)
            w_im = -1 * math.sin(omega)

            t_re = math.cos(2.*omega)
            t_im = -1 * math.sin(2.*omega)

            s_t_re = 1.
            s_t_im = 0.
            s_h_re = self.original_coeffs[0]
            s_h_im = 0.

            #for m in range(0, self.num_coefficients-1):
            m = 0
            while m < self.num_coefficients-1:

                if self.roots[m].imag == 0.:
                    temp = complex(1 - self.roots[m].real * w_re, self.roots[m].real * w_im)
                else:
                    temp_temp = math.pow(abs(self.roots[m]), 2.)

                    temp = complex(1. - (2 * self.roots[m].real * w_re) + (temp_temp * t_re), temp_temp * t_im - 2 * self.roots[m].real * w_im)
                    m += 1

                H = H * temp
                m += 1

            for i in range(1, self.num_coefficients):
                s_temp = s_t_re * w_re - s_t_im * w_im
                s_t_im = s_t_re * w_im + s_t_im * w_re
                s_t_re = s_temp

                s_h_re += self.original_coeffs[i] * s_t_re
                s_h_im += self.original_coeffs[i] * s_t_im

            s_h_mag = math.sqrt(math.pow(s_h_re, 2.) + math.pow(s_h_im, 2.))
            h_mag = abs(H)
            H = complex(1., 0.)
            upper += s_h_mag * h_mag
            lower += math.pow(h_mag, 2.)

        self.gain = upper / lower
        print(f"gain = {self.gain}")





if __name__ == "__main__":
    regular_filter = [
        0.0003616671417606181,
        -0.00079401999138141,
        -0.0006358046027427011,
        0.005099253312893104,
        -0.00868131808575515,
        0.0025415463320707884,
        0.017470296034074767,
        -0.03822079444428064,
        0.029733487197886137,
        0.03308035489798998,
        -0.14192083414359752,
        0.24841654629966603,
        0.707099240102832,
        0.24841654629966603,
        -0.14192083414359752,
        0.03308035489798998,
        0.029733487197886137,
        -0.03822079444428064,
        0.017470296034074767,
        0.0025415463320707884,
        -0.00868131808575515,
        0.005099253312893104,
        -0.0006358046027427011,
        -0.00079401999138141,
        0.0003616671417606181,
    ]

    decomposed = Poly(np.asarray(regular_filter))
    decomposed.get_roots()
    decomposed.build_filters()
    decomposed.calculate_gain()

    pass
