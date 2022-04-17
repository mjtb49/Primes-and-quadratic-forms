import mpmath
import math


def conv_to_super(n):
    if n == 1 or n == 0:
        return ""

    up = list(map(chr, [8304, 185, 178, 179, 8308, 8309, 8310, 8311, 8312, 8313]))
    st = ""
    for i in str(n):
        st += up[int(i)]
    return st


class Polynomial:
    def __init__(self, coefficients):
        self.coefficients = coefficients

    def __mul__(self, other):
        deg1 = len(self.coefficients)
        deg2 = len(other.coefficients)
        new_thing = [0] * (deg1 + deg2 - 1)
        for i in range(deg1):
            for j in range(deg2):
                new_thing[i + j] += self.coefficients[i] * other.coefficients[j]
        return Polynomial(new_thing)

    __rmul__ = __mul__

    def __str__(self):
        st = ""
        for n in range(len(self.coefficients) - 1, -1, -1):
            if self.coefficients[n] != 1 and self.coefficients[n] != 0:
                if self.coefficients[n] == -1:
                    st += "-"
                else:
                    st += str(self.coefficients[n]).split('.')[0]
            if n > 0 and self.coefficients[n] != 0:
                st += "x" + conv_to_super(n) + " + "
        return st


class QuadraticForm:
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c
        self.disc = b * b - 4 * a * c

    def value(self, x, y):
        return self.a * x * x + self.b * x * y + self.c * y * y

    def is_reduced_pos_def(self):
        flag = True
        if abs(self.b) == self.a or self.a == self.c:
            flag = self.b >= 0
        return self.is_primitive() and self.is_positive_definite() and abs(self.b) <= self.a <= self.c and flag

    def is_positive_definite(self):
        return self.a > 0 and self.disc < 0

    def is_primitive(self):
        return abs(math.gcd(self.a, self.b, self.c)) == 1

    def get_ideal_class_rep(self):
        return (-self.b + mpmath.sqrt(self.disc)) / (2 * self.a)

    def get_genus_representative(self):
        for x in range(0, -self.disc):
            for y in range(0, -self.disc):
                k = self.value(x, y) % -self.disc
                if abs(math.gcd(k, -self.disc)) == 1:
                    return k

    def get_genus_for_principal(self):
        genus = set()
        if self.a == 1 and self.b == 0:
            for x in range(2 * self.c):
                k = self.value(x, 0) % -self.disc
                if math.gcd(k, -self.disc) == 1:
                    genus.add(k)
                k = self.value(x, 1) % -self.disc
                if math.gcd(k, -self.disc) == 1:
                    genus.add(k)
        return genus

    def genus_size(self):
        forms = narrow_class_group(self.disc)
        genus = self.get_genus_for_principal()
        forms_found = 0
        for form in forms:
            if form.get_genus_representative() in genus:
                forms_found += 1
        return forms_found

    def is_convenient(self):
        return self.genus_size() == 1

    def print_residue_string(self):
        residues = []
        if self.disc < 0:
            for i in range(1, -self.disc + 1):
                if math.gcd(i, -self.disc) == 1 and jacobi_symbol(self.disc, i) == 1:
                    residues.append(i)
        print(f"Primes p not dividing {-self.disc} must leave residue mod {-self.disc} in {residues}")

    def print_genus_string(self):        
        print(f"Primes p not dividing {-self.disc} must leave residue mod {-self.disc} in "
              f"{sorted(list(self.get_genus_for_principal()))}")

    def __str__(self):
        coefficients = [self.c, self.b, self.a]
        st = ""
        for n in range(len(coefficients) - 1, -1, -1):
            if coefficients[n] != 1 and coefficients[n] != 0:
                if coefficients[n] == -1:
                    st += "-"
                else:
                    st += str(coefficients[n])
            if n > 0 and coefficients[n] != 0:
                st += "x" + conv_to_super(n)
            if 2 - n > 0 and coefficients[n] != 0:
                st += "y" + conv_to_super(2 - n)
            if n > 0 and coefficients[n] != 0:
                st += " + "
        return st


def jacobi_symbol(n, k):
    if k < 0 or k % 2 == 0:
        print("Bad args for jacobi")
    n = n % k
    t = 1
    while n != 0:
        while n % 2 == 0:
            n = n // 2
            r = k % 8
            if r == 3 or r == 5:
                t = -t
        n, k = k, n
        if n % 4 == 3 and k % 4 == 3:
            t = -t
        n = n % k
    if k == 1:
        return t
    return 0


def narrow_class_group(discriminant):
    a = 1
    reduced_forms = []
    while a * a <= -discriminant / 3:
        for b in range(-a, a + 1):
            c = (b * b - discriminant) // (4 * a)
            form = QuadraticForm(a, b, c)
            if form.disc == discriminant and form.is_reduced_pos_def():
                reduced_forms.append(form)
        a += 1
    return reduced_forms


def narrow_class_number(discriminant):
    return len(narrow_class_group(discriminant))


def make_integer_poly(poly: Polynomial):
    c = poly.coefficients
    clean = [0] * len(c)
    for i in range(len(c)):
        number = c[i]
        if abs(mpmath.im(number)) > 0.00000000001:
            print("Large imaginary part")
        clean[i] = mpmath.nint(mpmath.re(number))
        if abs(c[i] - clean[i]) > 0.00000000001:
            print("Nonintegral real part")
            print(i)
            print(c[i])
            print(clean[i])
    return Polynomial(clean)


def make_min_poly_hilbert_class(form_class_group):
    result = Polynomial([1])
    for g in form_class_group:
        result = result * Polynomial([-1728 * mpmath.kleinj(g.get_ideal_class_rep()), 1])
    return result


def print_prime_condition(n):
    form = QuadraticForm(1, 0, n)
    form_class_group = narrow_class_group(form.disc)
    class_number = len(form_class_group)
    if class_number == 1:
        print(str(form) + " is convenient and class number 1")
        form.print_residue_string()
    elif form.is_convenient():
        print(f"{form} is convenient of class number {class_number}")
        form.print_genus_string()
    else:
        print(f"{form} is class number {class_number} with genus size {form.genus_size()}")
        result = make_min_poly_hilbert_class(form_class_group)
        form.print_genus_string()
        print(f"Furthermore, there must exist solutions mod p to {make_integer_poly(result)}")


def main():
    mpmath.mp.dps = 1000
    max_class = 1
    for n in range(1, 2000):
        print("---------------------------------------------------")
        print_prime_condition(n)
        # form = QuadraticForm(1, 0, n)
        # group = narrow_class_group(form.disc)
        # # if len(group) > max_class:
        # max_class = len(group)
        # print(form)
        # result = Polynomial([1])
        # for g in group:
        #     result = result * Polynomial([-1728 * mpmath.kleinj(g.get_ideal_class_rep()), 1])
        # mpmath.mp.pretty = True
        # make_integer_poly(result)
        # form.print_residue_string()
        # print("")


if __name__ == '__main__':
    main()
