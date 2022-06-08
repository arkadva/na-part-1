import math
import sympy as sp

z = sp.symbols('z')


def f(function, val):
    return function.subs([(z, val)])

class Result:
    def __init__(self, iterations, x):
        self.iterations = iterations
        self.x = x


def derivative(polynomial):
    return polynomial.diff(z)


def bisection_method(polynomial, start_point, end_point, epsilon):
    a = start_point
    b = end_point
    c = (a + b) / 2
    max_iterations = math.ceil(-math.log10(epsilon / (b - a))/math.log10(2))

    for i in range(max_iterations + 1):
        c = (a + b) / 2
        if abs(f(polynomial, c)) < epsilon or (b - a)/2 < epsilon:
            return Result(i + 1, c)
        if f(polynomial, a) * f(polynomial, c) < 0:
            b = c
        else:
            a = c
    return None


def newton_method(polynomial, start_point, end_point, epsilon):
    polynomial_derivative = derivative(polynomial)
    x = (start_point + end_point)/2
    print(f'x = ({start_point} + {end_point})/2')
    for i in range(30):
        print(f'x = {x} - {f(polynomial, x)}/{f(polynomial_derivative, x)} = {x - f(polynomial, x)/f(polynomial_derivative, x)}')
        x = x - f(polynomial, x)/f(polynomial_derivative, x)
        if abs(f(polynomial, x)) < epsilon:
            return Result(i + 1, x)
    return None


def secant_method(polynomial, start_point, end_point, epsilon):
    print(f'a = {start_point}')
    a = start_point
    print(f'b = {end_point}')
    b = end_point
    i = 0
    while abs(f(polynomial, b) - f(polynomial, a)) > epsilon:
        print(f'a = {b}')
        print(f'b = ({a}*{f(polynomial, b)} - {b}*{f(polynomial, a)})/({f(polynomial, b)} - {f(polynomial, a)}) = {(a*f(polynomial, b) - b*f(polynomial, a))/(f(polynomial, b) - f(polynomial, a))}')
        print()
        b, a = (a*f(polynomial, b) - b*f(polynomial, a))/(f(polynomial, b) - f(polynomial, a)), b
        if abs(f(polynomial, b)) < epsilon:
            return Result(i + 1, b)
        i += 1
    return None


def find_max(function):
    function_d = function.diff(z)
    res = sp.solve(function_d, z)
    for result in res:
        if f(function_d.diff(z), result) < 0:
            return result
    return None


def find_max_val(function):
    max_x = find_max(function)
    if max_x is not None:
        return f(function, max_x)
    else:
        return 1

def simpson_method(polynomial, partition, start, end):
    integration_sum = 0
    step = (end - start)/partition
    p_x = start

    coefficient = [4, 2]
    coefficient_pointer = 0

    print(f'int_sum = {f(polynomial, p_x)}')
    integration_sum += f(polynomial, p_x)
    p_x += step

    for i in range(1, partition):
        print(f'int_sum = {integration_sum} + {coefficient[coefficient_pointer] * f(polynomial, p_x)}')
        integration_sum += coefficient[coefficient_pointer] * f(polynomial, p_x)
        coefficient_pointer ^= 1
        p_x += step
    print(f'int_sum = {integration_sum} + {f(polynomial, p_x)} = {integration_sum}')
    integration_sum += f(polynomial, p_x)
    print(f'Final result = {integration_sum} * {step}/3 = {step/3 * integration_sum}\n')
    deg_4_derivative = polynomial.diff(z).diff(z).diff(z).diff(z)
    tolerance = step**5/90 * partition * f(deg_4_derivative, find_max_val(polynomial))/2
    print("Tolerance with respect to the chosen partition:", tolerance)
    return step/3 * integration_sum

def main():
    methods = [bisection_method, newton_method, secant_method]
    def_epsilon = 0.00001

    polynomial = z**2 * sp.E**(-z**2 + 5*z - 3) * (z**2 + 3*z - 5)
    polynomial_derivative = derivative(polynomial)

    start_point = 0
    end_point = 1.5
    step = 0.1

    integral_start = 0.5
    integral_end = 1

    choice = int(input("1. Bisection method\n2. Newton–Raphson method\n3. Secant method\n4. Simpson's Integral\n5. Exit\nInput: ")) - 1
    while 0 <= choice <= 3:
        if choice == 3:
            print("Simpson's method:", simpson_method(polynomial, 10, integral_start, integral_end))
            choice = int(input("1. Bisection method\n2. Newton–Raphson method\n3. Secant method\n4. Simpson's Integral\n5. Exit\nInput: ")) - 1
            continue
        x = start_point
        results = []
        while x + step <= end_point:
            if f(polynomial, x) * f(polynomial, x + step) < 0:
                print("- Found possible root, f({}) * f({}) < 0".format(x, x + step))
                result = methods[choice](polynomial, x - step, x + step, def_epsilon)
                if result is not None:
                    print("- Result x = {} found, it took {} iterations.".format(result.x, result.iterations))
                    results.append(result.x)
                else:
                    print("Result not found, method failed.")
            elif f(polynomial_derivative, x) * f(polynomial_derivative, x + step) < 0:
                print("- Found possible root using the derivative, f'({}) * f'({}) < 0".format(x, x + step))
                result = methods[choice](polynomial_derivative, x - step, x + step, def_epsilon)
                if result is not None:
                    print("- Result x = {} found, it took {} iterations.".format(result.x, result.iterations), end=' ')
                    if abs(f(polynomial, result.x)) > def_epsilon:
                        print("But it's invalid as it's not a root for the original function.")
                    else:
                        results.append(result.x)
                        print()
                else:
                    print("Result not found, method failed.")
            x += step
        if len(results):
            print("\nFinal results:")
            for x in results:
                print("x = {}".format(x))
            print()
        choice = int(input("1. Bisection method\n2. Newton–Raphson method\n3. Secant method\n4. Simpson's Integral\n5. Exit\nInput: ")) - 1


if __name__ == '__main__':
    main()