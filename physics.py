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
from __future__ import division
import math

# Standard gravity.
g0 = 9.82 # closer to 9.8065 in the SI world, but it's 9.82 in KSP

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

def vectoradd(a, b):
    return [ ax + bx for (ax,bx) in zip(a,b) ]

def vectorscale(s, v):
    return [ s * vx for vx in v ]

class PiecewiseLinearCurve(object):
    def __init__(self, points):
        """
        Specify a sorted list of (t, value) pairs.
        """
        self.points = tuple(sorted(points, key = lambda x : x[0]))

    @staticmethod
    def fromfile(filename):
        pts = []
        with open(filename) as f:
            for line in f:
                (k, v) = line.split(' ')
                pts.append( (float(k),float(v)) )
        return PiecewiseLinearCurve(pts)

    def lookup(self, t):
        """
        Return the index of the first point whose parameter is greater than or
        equal to t.  Returns None if t is outside the bounds of the curve.
        """
        for (i, (k, v)) in enumerate(self.points):
            if k >= t:
                return i
        return None

    def evaluate(self, t):
        """
        Evaluate the curve at parameter position t.
        TODO: this is piecewise linear, whereas it's something else
        in-game.
        """
        i = self.lookup(t)
        if i is None:
            return self.points[-1][1]
        if i == 0:
            return self.points[0][1]

        # linear interpolation between i-1 and i
        (t0, v0) = self.points[i - 1]
        (t1, v1) = self.points[i]
        ratio = (t1 - t) / (t1 - t0)
        return ratio * v0 + (1-ratio) * v1

class AnimationCurve(object):
    """
    Attempted re-implementation of Unity's AnimationCurve.
    TODO: compute tangents automatically.
    """
    def __init__(self, data):
        """
        Data must be a list of quads (key, value, inTangent, outTangent)
        """
        for x in data: assert len(x) == 4
        # copy the data
        self._data = tuple(tuple(x) for x in sorted(data))

    def evaluate(self, key):
        # linear search for the first key bigger than the input
        def find():
            for i, (k,v,_,__) in enumerate(self._data):
                if key <= k: return i
            return -1
        i = find()

        # before the first, return the first value; similarly after the last
        if i == 0: return self._data[0][1]
        if i == -1: return self._data[-1][1]

        # cubic interpolation, formula from
        # http://answers.unity3d.com/questions/464782/t-is-the-math-behind-animationcurveevaluate.html
        t0, v0, _, t0out = self._data[i - 1]
        t1, v1, t1in, _ = self._data[i]
        if t0 == t1:
            return 0.5 * (v0 + v1)

        dt = t1 - t0
        t = (key - t0) / dt
        t2 = t * t
        t3 = t * t * t

        m0 = dt * t0out
        m1 = dt * t1in

        a = 2 * t3 - 3 * t2 + 1
        b = t3 - 2 * t2 + t
        c = t3 - t2
        d = -2 * t3 + 3 * t2

        return a * v0 + b * m0 + c * m1 + d * v1
