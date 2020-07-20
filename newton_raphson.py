from scipy import optimize
import matplotlib.pyplot as plt

def func(x) :
    return x + 1

a = optimize.newton(func, 1.5)

print(a)