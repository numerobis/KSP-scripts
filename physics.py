# KSP Physics basics.
# Copyright 2012 Benoit Hudson
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


"""
Fundamental constants in the KSP universe (at least in 0.18.1), and various math bits.
"""
import math

# Standard gravity.
g0 = 9.81 # closer to 9.8065 in the SI world.

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
