# KSP Engine descriptions.
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

from physics import g0

"""
This module encodes information about the engines available in KSP.

It also includes some functions related to the ideal rocket equation.

TODO: read this from the part.cfg files.
"""



###########################################################################
# This exception is thrown when we ask for the required propellant + tank mass,
# but it's infinite because of insufficient Isp.
class WeakEngineException(Exception):
    def __init__(self, Isp):
        self.Isp = Isp

###########################################################################

class engine(object):
    def __init__(self, name, IspAtm, IspVac, mass, thrust,
            vectoring = False, radial = False, large = False, fuel="liquid"):
        self.name = name
        self.IspAtm = IspAtm        # seconds
        self.IspVac = IspVac        # seconds
        self.fuel = fuel            # "liquid" or "solid"
        self.mass = mass            # tonnes
        self.thrust = thrust        # kN
        self.vectoring = vectoring  # true or false
        self.radial = radial        # true of false
        self.large = large          # true: 2m, false: 1m (can use bi/tricoupler)

    def __str__(self): return self.name

    def Isp(self, planet, altitude):
        # Assumption: Isp is in a linear correspondence with pressure,
        # clipped to 1 Atm (as determined by experiments on Kerbin and Eve).
        #
        # This is patently false for the jets.
        #
        if planet is None or altitude is None:
            return self.IspVac
        pressure = planet.pressure(altitude)
        if pressure > 1: pressure = 1
        Isp = pressure * self.IspAtm + (1.0 - pressure) * self.IspVac
        return Isp

    def __str__(self):
        return self.name

types = (
    # Bipropellent engines.
    engine("24-77",     250, 300,    0.09,     20, vectoring=True, radial=True),
    engine("48-7S",     300, 350,    0.1,      30, vectoring=True),
    engine("Aerospike", 388, 390,    1.5,     175),
    engine("KR-1x2",    290, 340,   10.0,    2000, vectoring=True, large=True),
    engine("KR-2L",     280, 380,    6.5,    2500, vectoring=True, large=True),
    engine("KS-25x4",   320, 360,    9.75,   3200, vectoring=True, large=True),
    engine("LV-1",      220, 290,    0.03,      4, vectoring=True),
    engine("LV-1R",     220, 290,    0.03,      4, vectoring=True, radial=True),
    engine("LV-909",    300, 390,    0.5,      50, vectoring=True),
    engine("LV-N",      220, 800,    2.25,     60, vectoring=True),
    engine("LV-T30",    320, 370,    1.25,    215),
    engine("LV-T45",    320, 370,    1.5,     200, vectoring=True),
    engine("Mainsail",  320, 360,    6,      1500, vectoring=True, large=True),
    engine("Mark 55",   320, 360,    0.9,     120, vectoring=True, radial=True),
    engine("Poodle",    270, 390,    2,       220, vectoring=True, large=True),
    engine("RAPIER",    320, 360,    1.2,     175, vectoring=True),
    engine("Skipper",   320, 370,    3,       650, vectoring=True, large=True),

    # Solids
    engine("Launch Escape System",
                        320, 360,    1,       750, fuel="solid"),
    engine("SRB-KD25k", 230, 250,    3,       650, fuel="solid"),
    engine("BACC",      230, 250,    1.5,     315, fuel="solid"),
    engine("RT-10",     225, 240,    0.5,     250, fuel="solid"),
    engine("Sepratron", 100, 100,    0.0125,   18, fuel="solid"),
)
_engineDict = dict((e.name.lower(), e) for e in types)
def getEngine(name):
    return _engineDict[name.lower()]


# To help the heuristics, choose the best possible Isp at a given altitude.
def maxIsp(planet, altitude):
    ispMaxEngine = max(types, key = lambda x: x.Isp(planet, altitude))
    return ispMaxEngine.Isp(planet, altitude)


# To help the heuristics, choose the best possible mass to achieve a given
# thrust.
maxThrustPerMassEngine = max(types, key =
        lambda x: 0 if x.thrust == 0 else x.thrust / x.mass)
def lightestEngineForThrust(thrust):
    num = thrust / maxThrustPerMassEngine.thrust
    return (maxThrustPerMassEngine, num)


##############################
## Tsiokolvsky rocket equation

beta = 8 # ratio of propellant mass : dry mass in the big stock tanks.

def deltaV(m1, m0, Isp):
    Ve = Isp * g0
    return Ve * math.log(m1 / m0)

def alpha(deltaV, Isp):
    # Corner cases: if deltaV is 0, alpha is 1: m1 = m0 obviously.
    # If Isp is 0, alpha is undefined; raise an exception
    if deltaV == 0: return 1
    if Isp == 0: raise WeakEngineException(Isp)

    return math.exp(deltaV / (Isp * g0))

def propellantMass(deltaV, Isp, m0):
    return m0 * (alpha(deltaV, Isp) - 1)

def propellantMassBurned(deltaV, Isp, m1):
    # m1 = m0 alpha
    # m0 = m1 / alpha
    # mprop = m1 - m0
    # mprop = m1 (1 - 1/alpha)
    return m1 * (1 - 1.0 / alpha(deltaV, Isp))

def payloadAllowance(deltaV, Isp, mprop):
    """
    Say we burn mprop tonnes of fuel in a burn that achieves deltaV m/s,
    using an engine of the given Isp.

    What's the payload?
    """
    # m1 = m0 alpha
    # mprop = m1 - m0
    # mprop = m0 alpha - m0
    # m0 = mprop / (alpha - 1)
    return mprop / (alpha(deltaV, Isp) - 1.0)

def burnMass(deltaV, Isp, m0, b = beta):
    """
    Return the mass of propellant and tanks that we'll need to burn.

    Assume tanks hold beta times their mass.  The assumption is false
    for some of the smallest tanks.

    Raise WeakEngineException if it is impossible to achieve the deltaV with
    that given Isp.  While the ideal rocket equation allows arbitrary deltaV,
    it doesn't take account of tank dry mass, which grows as propellant mass
    grows.

    deltaV: m/s
    Isp: s
    m0: tonnes, should include payload, engines, decouplers, but
        not the propellant and tanks we're using.

    return (propellant mass, tank mass)
    """
    # The amount of fuel we need is derived from the ideal rocket
    # equation.  Let alpha = e^{deltaV/Isp.g0}, beta = ratio of
    # propellant to dry mass of a tank.  Then:
    # tankMass = (alpha-1) * payload / (1 - alpha + beta)
    # Where the relevant payload here includes engines and decouplers.
    # Clearly, if (1-alpha+beta) <= 0 then we are in an impossible
    # state: this corresponds to needing infinite fuel.
    a = alpha(deltaV, Isp)
    if 1 - a + b <= 0: raise WeakEngineException(Isp)
    tankMass = m0 * (a - 1) / (1 - a + b)
    propMass = tankMass * b
    return (propMass, tankMass)



def burnTimeForPayload(deltaV, Isp, thrust, m0):
    """
    Return the time needed to perform the burn.
    m0 is the dry mass including tanks.

    Raise WeakEngineException if there is no thrust.
    """
    # If there's no deltaV, or no mass, we don't burn at all.
    if deltaV == 0 or m0 == 0: return 0

    # If there's no thrust but we want to move, problem...
    if Isp == 0 or thrust == 0:
        raise WeakEngineException(Isp)

    # The mass flow rate of an engine is thrust / (Isp * g0)
    # The mass we expel is from the ideal rocket equation.
    # The amount of time we burn is thus mprop / (thrust/(Isp*g0))
    mprop = propellantMass(deltaV, Isp, m0)
    return mprop * Isp * g0 / thrust

def burnTimeFromFull(deltaV, Isp, thrust, m1):
    a = alpha(deltaV, Isp)
    mprop = (1.0 - 1.0/a) * m1
    return mprop * Isp * g0 / thrust

def minThrustForBurnTime(deltaV, Isp, m0, time):
    """
    Return the minimum thrust necessary to achieve a burn in the time
    limit.
    m0 is the dry mass including tanks.
    """
    m1 = propellantMass(deltaV, Isp, m0)
    return m1 * Isp * g0 / time

def combineIsp(engines, planet, altitude):
    """
    Given a dictionary mapping engine -> count, or a list of pairs, compute the
    Isp of the system at the given altitude on the given planet.  Pass in None
    for the planet to get vacuum values.

    This is the weighted sum of the impulses of each type, with weights
    based off the relative mass flow rate of the engines.  Proving this
    requires re-deriving the ideal rocket equation, but with two engines,
    and generalizing in the obvious way.
    """
    def alpha(e, c):
        # no thrust => no contribution (even if Isp is zero)
        if c == 0 or e.thrust == 0: return 0
        else: return e.thrust * c / e.Isp(planet, altitude)

    try:
        totalthrust  = sum(e.thrust * c for (e,c) in engines.iteritems())
        totalweights = sum(alpha(e,c)   for (e,c) in engines.iteritems())
    except AttributeError:
        totalthrust  = sum(e.thrust * c for (e,c) in engines)
        totalweights = sum(alpha(e,c)   for (e,c) in engines)

    # no thrust at all => Isp may as well be zero
    if totalthrust == 0:
        return 0
    else:
        return totalthrust / totalweights


# combine 2x nuke and 1x sail
# print combineIsp( [(getEngine("mainsail"), 1), (getEngine("lv-n"), 2)], None, None)



###########################################################################
## Notes about experiments:
## I ran an experiment to see how Isp varied for the nuke as altitude
## changes.  Here are the results.
# Explanation for the systematic error, in part, is that the altimeter measures
# about 20m off from the engine -- below during launch on Kerbin, above while
# crashing on Eve.

# On Kerbin:
# alt  |    IspExperiment   | IspCalculated | Error
# kerbinaltIsp =
# (83, 229.6, 229.5, -0.05),
# (530, 279.7, 278.3, -1.37),
# (1020, 328.5, 327.0, -1.47),
# (2010, 413.5, 412.0, -1.51),
# (3020, 484.4, 483.0, -1.44),
# (4030, 542.3, 540.9, -1.35),
# (5034, 589.3, 588.1, -1.22),
# (6100, 629.9, 628.8, -1.13),
# (7700, 676.7, 675.7, -1.04),
# (9135, 707.6, 706.7, -0.92),
# (10767, 733.5, 732.7, -0.83),
# (13400, 760.7, 760.2, -0.47),
# (16874, 780.5, 780.1, -0.35),
# (23600, 795.0, 794.8, -0.17),
# (31147, 798.9, 798.9, -0.04),
# (35833, 799.6, 799.6, -0.05),
# (47552, 800.0, 800.0, -0.04),

# On Eve:
# evealtIsp = (
# (78213, 800.0, 800.0, -0.04),
# (76300, 799.9, 799.9, 0.05),
# (67000, 799.8, 799.8, -0.00),
# (61261, 799.5, 799.5, 0.04),
# (50040, 798.0, 797.7, -0.28),
# (39847, 790.1, 790.2, 0.12),
# (33970, 777.2, 777.4, 0.16),
# (28439, 749.6, 750.1, 0.51),
# (23500, 699.4, 699.0, -0.42),
# (18714, 599.0, 599.9, 0.86),
# (15946, 501.7, 502.8, 1.08),
# (13845, 397.5, 398.7, 1.24),
# (12351, 301.9, 303.3, 1.37),
# (11635, 248.3, 249.8, 1.48),
# (11271, 220.0, 220.4, 0.41),
# (9528, 220.0, 220.0, 0.00),
# (5100, 220.0, 220.0, 0.00),
# (4000, 220.0, 220.0, 0.00),
# )
