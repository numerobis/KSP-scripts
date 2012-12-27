# KSP Planets module
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

from __future__ import division
import math
from math import sqrt, cos, sin, exp, log, pi

from physics import g0, L2, quadratic

"""
This module provides information about the planets and other major bodies in
the KSP solar system.
"""

# The drag model in KSP 0.18.1 is
#       F_{drag} = gamma D v^2 m
# where D is atmospheric density, v is the speed, m the mass (as a standin
# for cross-section).  gamma combines a bunch of coefficients in one
# (including the 1/2 that would normally be there).  It is not actually
# a constant, but it varies little over reported terminal velocities and over
# the parts we actually use (though aero components have lower drag).
coefficientOfDrag = 0.2 # Assumed constant, not really the case
dragMultiplier = 0.008
kerbinSurfaceDensity = 1.2230948554874
gamma = 0.5 * coefficientOfDrag * dragMultiplier * kerbinSurfaceDensity


class planet(object):

    def __init__(self, name, gravityParam, radiusKm, siderealPeriod, datumPressure, scale):
        self.name = name
        self.datumPressure = datumPressure # atm
        self.scale = scale      # m (pressure falls by e every scale altitude)
        # Use siderealPeriod because those are more exact values
	# 2 * pi * r = circumference, divide by period to get m/s
        self.siderealRotationSpeed = 2 * pi * radiusKm * 1000 / siderealPeriod # m/s
        self.radius = radiusKm * 1000 # stored in m
        self.mu = gravityParam  # m^3/s^2

    def __str__(self): return self.name

    def gravity(self, altitude = 0):
        """
        Return the gravitational acceleration at a given altitude above the
        surface.

        altitude is in meters
        """
        r = self.radius + altitude
        return self.mu / (r*r)

    def orbitalVelocity(self, altitude):
        """
        Return the velocity required to maintain an orbit with the given
        altitude at the semimajor axis (e.g. a circular orbit at that
        altitude).

        altitude is in meters.
        """
        return sqrt(self.mu / (self.radius + altitude))

    def escapeVelocity(self, altitude):
        """
        At a given altitude, return the minimum velocity that will result in
        escape.

        altitude is in meters.
        """
        return sqrt(2 * self.mu / (self.radius + altitude))

    def determineOrbit(self, h, velocity):
        """
        Given distance from the planet's core, and current velocity vector,
        return the apsides in meters from the planet's core.

        If one of them is negative, the orbit is hyperbolic.

        h is in meters relative to the core
        velocity is in polar coordinates (m/s, radians) where the angle is
            relative to the tangent to the surface.
        """
        # http://forum.kerbalspaceprogram.com/showthread.php/36178-How-to-predict-my-apoapsis-given-altitude-and-velocity
        # Plugging in Horn Brain's derivation (which is the same as Jason
        # Patterson's, but seems a bit easier to follow):
        # We use two equations, twice each:
        # vis viva:     v^2 = mu (2/r - 1/a) where a is the semimajor axis
        # momentum:     p = r v cos(theta) where theta is the angle of v to
        #               the surface tangent.
        # Momentum is preserved, so p is the same value at all points of
        # the orbit.
        #
        # Steps:
        # 1. vis viva with r=h, v=v, solve for 1/a.
        # 2. momentum from the parameters: p = rv cos theta
        # 3. momentum at the apses to relate r_apsis and v_apsis
        # 4. vis viva with r being the apses.
        #    Substitue 1/a from before, and substitute v_apsis from
        #    momentum.
        # 5. Solve quadratic equation:
        #    0 = (1/a) r_apsis^2 - 2 r_apsis + p^2 / mu
        #
        # We get two solutions, which are the apoapsis and periapsis
        # (measured from the center of the planet), in that order.
        #
        (v, theta) = velocity
        ainv = (2 / h) - (v*v / self.mu)
        p = h * v * cos(theta)
        return quadratic(ainv, -2, p*p / self.mu)

    def determineOrbit2(self, position, velocity):
        """
        Given a position and velocity in some coordinate frame centered about
        the center of the planet (doesn't matter, as long as 1 unit is 1
        meter), return the apsides, in meters above the core.
        If one of them is negative, the orbit is hyperbolic.
        """

        # theta is the angle of the position above the horizontal.
        # phi is the angle of the velocity above the horizontal.
        # psi is the angle of the velocity above the surface tangent.
        # psi = (90 - theta) + phi
        # thetabar = 90-theta
        theta = math.atan2(position[1], position[0])
        phi   = math.atan2(velocity[1], velocity[0])
        thetabar = (math.pi/2) - theta
        psi = thetabar + phi

        h = L2(position)
        v = L2(velocity)
        return self.determineOrbit(h, (v, psi))

    def terminalVelocity(self, altitude):
        """
        Return the terminal velocity at a given altitude.
        This is the speed at which drag is as strong as gravity.

        altitude is in meters.
        """

        # v = sqrt(g/ (datumPressure * gamma * e^{-altitude/scale}))
        # Then pull the e term out, and remember that gravity changes
        # (slightly) with altitude.
        return exp(0.5 * altitude / self.scale) * sqrt(
                    self.gravity(altitude) / (gamma * self.datumPressure) )

    def drag(self, altitude, velocity):
        """
        Return the acceleration due to drag, in m/s^2.

        This is only possible since the drag force depends on mass;
        otherwise we'd need to know the mass of the vessel also.

        altitude is in meters, velocity in m/s
        """
        # F_{drag} = D v^2 m gamma, where D is atmospheric pressure, v is
        # the speed, m the mass (as a standin for cross-section), and gamma
        # combines a bunch of variables that are close enough to constant.
        # To get the acceleration, divide by mass, which cancels out m:
        # a_{drag} = gamma P v^2
        return gamma * self.pressure(altitude) * velocity * velocity

    def pressure(self, altitude):
        """
        Return the atmospheric pressure in Atm.

        altitude is in meters
        """
        if altitude is None or altitude >= self.topOfAtmosphere(): return 0
        return self.datumPressure * exp(-altitude / self.scale)

    def altitude(self, pressure):
        """
        Return the altitude at which the given atmospheric pressure prevails.

        pressure is in kerbin standard atmospheres
        """
        if pressure < 1e-6 * self.datumPressure:
            return None
        # p = datum * e^(-alt/scale)
        # scale * log(datum/p) = alt
        return self.scale * log(self.datumPressure / pressure)

    def topOfAtmosphere(self):
        """
        Return the altitude of the top of the atmosphere, defined
        as being the altitude where atmospheric pressure falls below
        one millionth the pressure at the datum.

        returns in meters
        """
        # self.altitude with pressure = self.datumPressure * 1e-6, inlined:
        # log(1e6) ~ 13.81551...
        return self.scale * 13.81551

kerbin = planet("Kerbin", 3531600000000, 600, 21600,  1,   5000)
eve    = planet("Eve",    8171730200000, 700, 80500,  5,   7000)
laythe = planet("Laythe", 1962000000000, 500, 52980.8790593796,  0.8, 4000)
jool   = planet("Jool",     2.82528E+14, 600, 36000, 15,   9000)

_planets = dict([ (p.name.lower(), p) for p in (kerbin, eve, laythe, jool) ])
def getPlanet(name):
    return _planets[name.lower()]
