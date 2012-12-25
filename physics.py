"""
Fundamental constants in the KSP universe (at least in 0.18.1)
"""
import math

# Standard gravity.
g0 = 9.80665

def quadratic(a, b, c):
    """
    Return both real roots of ax^2 + bx + c = 0, with the positive (or less
    negative) one first.
    DomainError if they are not real.
    """
    # (-b +/- sqrt(b^2 - 4ac)) / 2a
    sqrtdiscriminant = math.sqrt(b*b - 4*a*c)
    a2 = 2 * a
    return ( (-b + sqrtdiscriminant) / a2, (-b - sqrtdiscriminant) / a2)

def L2(vector):
    return math.sqrt(sum(x*x for x in vector))

def L2_2(vector):
    return sum(x*x for x in vector)
