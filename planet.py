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

    def __init__(self, name, gravityParam, SOI, radiusKm, siderealPeriod, datumPressure, scale):
        self.name = name
        self.mu = gravityParam          # m^3/s^2
        self.SOI = SOI                  # m
        self.radius = radiusKm * 1000   # stored in m, received in km

        self.siderealPeriod = siderealPeriod    # s
        self.siderealRotationSpeed = 2 * pi * self.radius / siderealPeriod # m/s
        self.datumPressure = datumPressure      # atm
        self.scale = scale      # m (pressure falls by e every scale altitude)

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
        escape from this gravity well.  Note that we can likely leave the
        sphere of influence with much less speed.

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

    def hohmann(self, a1, a2):
        """
        Given a circular orbit at altitude a1 and a target orbit at altitude
        a2, return the delta-V budget of the two burns required for a Hohmann
        transfer.

        a1 and a2 are in meters above the surface.
        returns a pair (dV1, dV2) corresponding to the two burns required.
            if a1 < a2, burn prograde; else burn retrograde
        """
        r1 = a1 + self.radius
        r2 = a2 + self.radius

        # source: wikipedia
        sqrt_r1 = math.sqrt(r1)
        sqrt_r2 = math.sqrt(r2)
        sqrt_2_sum = math.sqrt(2 / (r1 + r2))
        sqrt_mu = math.sqrt(self.mu)
        dV1 =   sqrt_mu / sqrt_r1 * (sqrt_r2 * sqrt_2_sum - 1)
        dV2 =   sqrt_mu / sqrt_r2 * (1 - sqrt_r1 * sqrt_2_sum)
        return (dV1, dV2)

    def bielliptic(self, a1, a2, aintermediate):
        """
        Given a circular orbit at altitude a1 and a target orbit at altitude a2,
        along with an intermediate altitude, return the delta-V budget of the
        three burns required for a bielliptic transfer.

        a1, a2, and aintermediate are in meters above the surface
        returns a triple (dV1, dV2, dV3) corresponding to the three burns required.
            if a1 < a2 then burn dV1 and dV2 prograde, and dV3 retrograde
            else burn dV1 prograde, and dV2 and dV3 retrograde
        """
        r0 = a1 + self.radius
        rb = aintermediate + self.radius
        rf = a2 + self.radius
        print ("bielliptic: start %g, intermediate %g, end %g"
            % (r0, rb, rf))

        # source: wikipedia and some simplification
        mu = self.mu
        two_r0rb = 2 / (r0 + rb)
        sqrt_2_rbrf = math.sqrt(2 * rf / (rb + rf))
        dV1 = math.sqrt(mu / r0) * (math.sqrt(rb * two_r0rb) - 1)
        dV2 = math.sqrt(mu / rb) * (sqrt_2_rbrf - math.sqrt(r0 * two_r0rb))
        dV3 = math.sqrt(mu / rf) * (1 - sqrt_2_rbrf)
        return (dV1, dV2, dV3)

    def soiBurn(self, altitude, soiSpeed):
        """
        If you are in a circular orbit at the given altitude and you want to burn enough
        to reach the given speed at the sphere of influence, how much should you burn now?

        Returns m/s of deltaV required.
        """
        # source: http://forum.kerbalspaceprogram.com/showthread.php/16511-Interplanetary-How-To-Guide
        r1 = altitude + self.radius
        r2 = self.SOI
        v2 = soiSpeed
        mu = self.mu
        # v = math.sqrt( (r1 * (r2 * v2 * v2 - 2 * mu) + 2 * r2 * mu) / (r1 * r2) )
        v = math.sqrt(v2 * v2 + 2 * mu * (r2 - r1) / (r1 * r2))
        return v - self.orbitalVelocity(altitude)


    def terminalVelocity(self, altitude):
        """
        Return the terminal velocity at a given altitude.
        This is the speed at which drag is as strong as gravity.

        altitude is in meters.
        """

        # v = sqrt(g/ (datumPressure * gamma * e^{-altitude/scale}))
        # Then pull the e term out, and remember that gravity changes
        # (slightly) with altitude.
        if altitude is None or altitude >= self.topOfAtmosphere(): return float("inf")
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
        if altitude is None or altitude >= self.topOfAtmosphere(): return 0
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
        if pressure < 1e-6 * self.datumPressure or self.datumPressure is 0:
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


# the sun has an infinite SOI; arbitrarily set it to 500x the orbit of Jool
sunSOI = 500 * 71950638386

planets = dict([ (p.name.lower(), p) for p in (
    #      name,     mu,             SOI,          radiusKm, sidereal, atm, scale
    planet("Kerbol", 1.172332794E+18,        sunSOI, 261600,  432000,   0,    0),
    planet("Moho",      245250003655,   11206449   ,    250, 1210000,   0,    0),
    planet("Eve",      8171730229211,   85109365   ,    700,   80500,   5, 7000),
    planet("Gilly",          8289450,     126123.27,     13,   28255,   0,    0),
    planet("Kerbin",      3.5316E+12,   84159286   ,    600,   21600,   1, 5000),
    planet("Mun",        65138397521,    2429559.1 ,    200,  138984,   0,    0),
    planet("Minmus",      1765800026,    2247428.4 ,     60,   40400,   0,    0),
    planet("Duna",      301363211975,   47921949   ,    320,   65518, 0.2, 3000),
    planet("Ike",        18568368573,    1049598.9 ,    130,   65518,   0,    0),
    planet("Dres",       21484488600,   32700000   ,    138,   34800,   0,    0),
    planet("Jool",   282528004209995, 2455985200   ,    600,   36000,  15, 9000),
    planet("Laythe",   1962000029236,    3723645.8 ,    500,   52981, 0.8, 4000),
    planet("Vall",      207481499474,    2406401.4 ,    300,  105962,   0,    0),
    planet("Tylo",     2825280042100,   10856518   ,    600,  211926,   0,    0),
    planet("Bop",         2486834944,     993002.8 ,     65,  544507,   0,    0),
    planet("Pol",          227905920,    1041613   ,     44,  901903,   0,    0),
    planet("Eeloo",      74410814527,  119087000   ,    210,   19460,   0,    0),
) ])
del sunSOI

def getPlanet(name):
    return planets[name.lower()]

# register the planets with the module, so users can write planet.kerbin or planet.getPlanet("kerbin") equivalently
for (name, p) in planets.iteritems():
    globals()[name] = p
