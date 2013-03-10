# KSP Lander designer module
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

import planet
import math
import engine

partdrag = 0.2
droguedrag = 100
chutedrag = 500

def landingspeed(m, chutemass, droguemass, planet = planet.kerbin, altitude = 0):
    """
    Assuming a lander landing in an atmosphere with the given mass (not
    including parachutes), parachutes, and drogues, what is the landing speed?
    """
    rho = planet.airdensity(altitude)
    g = planet.gravity(altitude)
    gforce = g * (m + droguemass + chutemass)
    area = partdrag * m + droguedrag * droguemass + chutedrag * chutemass
    return math.sqrt(gforce / (area * 0.008 * rho))


def chutemass(m, v, planet = planet.kerbin, altitude = 0):
    """
    Assuming a 0.2-drag lander, how many tonnes of parachutes 
    does it take to land at speed v (in m/s)?
    """
    rho = planet.airdensity(altitude)
    g = planet.gravity(altitude)
    Cdrhov2 = 0.008 * rho * v * v
    if partdrag * Cdrhov2 >= g:
        # We slow down below v even with no parachutes.
        return 0
    if chutedrag * Cdrhov2 <= g:
        # We can't slow down this much even with just parachutes!
        return float('inf')
    chuteratio = (partdrag * Cdrhov2 - g) / (g - chutedrag * Cdrhov2)
    return m * chuteratio

def powerburn(m, vland, Isp, chutemass, droguemass, 
        planet = planet.kerbin, altitude = 0):
    """
    Assuming a lander with the given mass, chutes, and drogues in the given
    atmospheric conditions, how much fuel and tank space is required to achieve
    the given speed?

    Returns the pair (mprop, mtank)

    Assumption:
    - The rocket engines have infinite thrust:
        - We stop instantaneously
        - We don't need additional engine mass
    - We are at terminal velocity and the target altitude upon turning on the
      rockets.
    """
    def calc(mprop, mtank):
        """
        Calculate the deltaV required to stop the spacecraft,
        given its propellant and tank space.
        Then check how much propellant and tank space we need to stop the
        spacecraft.
        """

        vterm = landingspeed(m + mprop + mtank, chutemass, droguemass, planet, altitude)
        dV = vterm - vland
        if dV <= 0:
            return (mprop, mtank)
        m0 = m + chutemass + droguemass
        return engine.burnMass(dV, Isp, m0)
    mproplast = 0
    (mprop, mtank) = calc(0, 0)
    while mprop > mproplast + 0.01:
        mproplast = mprop
        (mprop, mtank) = calc(mprop, mtank)
    return (mprop, mtank)

def semipowered(m, vland, Isp, planet = planet.kerbin, altitude = 0):
    """
    Suggests an "optimal" mass of parachutes and fuel to reduce landing speed
    to v.  Returns a pair (chutemass, fuelmass).

    Approximations:
    - The rocket engines have infinite thrust:
        - We stop instantaneously
        - We don't need additional engine mass
    - We are at terminal velocity upon turning on the rockets.
    - I try 0.1t of chutes, then 0.2, etc.  This means you can't have an odd
      number of radials.

    For safety, add ~100m-200m to the altitude compared to what you actually
    want to achieve.
    """
    # You can have x*0.1 + y*0.15 tonnes of chutes for any integer x,y,
    # equivalent to having z*0.05 tonnes for any integer z >= 2, or z=0.
    singlechutemass = 0.05
    minchutes = 2
    def masses(nchutes):
        chutemass = singlechutemass * nchutes
        (mprop, mtank) = powerburn(m, vland, Isp, chutemass, 0,
                planet=planet, altitude=altitude)
        return (chutemass, mtank, mprop)

    bestchutes = 0
    bestmass = sum(masses(0))
    for nchutes in xrange(minchutes, int(m / singlechutemass)):
        mass = sum(masses(nchutes))
        if mass < bestmass:
            bestchutes = nchutes
            bestmass = mass

    (mchutes, mtanks, mprop) = masses(bestchutes)
    return (mchutes, mprop)
