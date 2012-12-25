from __future__ import division
import math

from physics import g0

"""
This module encodes information about the engines available in KSP.

It also includes some functions related to the ideal rocket equation.
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
            vectoring = False, radial = False, large = False):
        self.name = name
        self.IspAtm = IspAtm        # seconds
        self.IspVac = IspVac        # seconds
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

# We have a none engine for fuel-only stages in asparagus staging.  It has no
# mass, no thrust, Isp doesn't matter.  Make it radial, so that we can have
# radial engines above and we needn't add towers.
noEngine = engine("none", 0, 0, 0, 0, radial=True)

types = (
    # One of the choices of engines is to have none...
    noEngine,

    # Jets.  TODO

    # Bipropellant engines, in order of Isp first, and thrust:mass ratio second
    engine("LV-N",      220, 800,    2.25,     60, vectoring=True),
    engine("Aerospike", 388, 390,    1.5,     175),
    engine("LV-909",    300, 390,    0.5,      50, vectoring=True),
    engine("Poodle",    270, 390,    2.5,     220, vectoring=True, large=True),
    engine("LV-T30",    320, 370,    1.25,    215),
    engine("LV-T45",    320, 370,    1.5,     200, vectoring=True),
    engine("Mainsail",  280, 330,    6,      1500, vectoring=True, large=True),
    engine("Mark 55",   290, 320,    0.9,     120, vectoring=True, radial=True),
    engine("24-77",     250, 300,    0.09,     20, vectoring=True, radial=True),
    engine("LV-1",      220, 290,    0.03,      1.5),

    # The ion engine is a bit off: we need a *lot* of power for it,
    # which starts to add more mass than just the 0.25.  Also, the
    # dry mass is much more than 1/4 of the propellant, and the number
    # of containers starts to matter.  So just ignore it.
    #engine("ion",      4200,4200,    0.25,      0.5),

    # Note: solid-fuel rockets aren't handled correctly yet.
    # engine("Sepratron", 100,    0.15,     20, solid=9),
    # engine("RT-10",     240,    0.5,     250, solid=433),
    # engine("BACC",      250,    1.75,    300, solid=850),
)

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
## Available tanks

#class tank(object):
#    def __init__(self, name, drymass, liquid, oxidizer):
#        """
#        Liquid and oxidizer in 5kg units (!?); drymass in tonnes.
#        """
#        self.name = name
#        self.drymass = drymass
#        self.propellant = drymass + (liquid+oxidizer) * 0.005
#        self.fullmass = drymass + self.propellant

## Store the tanks sorted by capacity.  But only the big tanks; the two
## small tanks I don't feel like handling.
#tanks = [
#    tank("Jumbo-64", 4, 2880, 3520),
#    tank("X200-32", 2, 1440, 1760),
#    tank("X200-16", 1, 720, 880),
#    tank("FL-T800 or X200-8", 0.5, 360, 440),
#    tank("FL-T400", 0.25, 180, 220),
#    tank("FL-T200", 0.125, 90, 110),
#    tank("Oscar-B", 0.015, 5.735, 7),
#    tank("ROUND-8", 0.025, 10, 12.2),
#]

##############################
## Tsiokolvsky rocket equation

beta = 8 # ratio of propellant mass : dry mass in the big stock tanks.

def alpha(deltaV, Isp):
    # Corner cases: if deltaV is 0, alpha is 1: m1 = m0 obviously.
    # If Isp is 0, alpha is undefined; raise an exception
    if deltaV == 0: return 1
    if Isp == 0: raise WeakEngineException(Isp)

    return math.exp(deltaV / (Isp * g0))

def propellantMass(deltaV, Isp, m0):
    return m0 * (alpha(deltaV, Isp) - 1)

def burnMass(deltaV, Isp, m0):
    """
    Return the mass of propellant and tanks that we'll need to burn.

    assuming tanks hold 8 times their mass.  The assumption is false
    for some of the smallest tanks.

    Return None if it is impossible to achieve the deltaV with that given
    Isp given the mass of the tanks.

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
    beta = 8 # TODO: handle the smaller, less efficient tanks too!
    if 1 - a + beta <= 0: raise WeakEngineException(Isp)
    tankMass = m0 * (a - 1) / (1 - a + beta)
    propMass = tankMass * beta
    return (propMass, tankMass)


def burnTime(deltaV, Isp, thrust, m0):
    """
    Return the time needed to perform the burn.
    m0 is the dry mass including tanks.
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
    Given a dictionary mapping engine -> count, compute the
    Isp of the system at the given altitude.

    This is the weighted sum of the impulses of each type, with weights
    based off the relative mass flow rate of the engines.
    """
    def alpha(e, c):
        # no thrust => no contribution (even if Isp is zero)
        if c == 0 or e.thrust == 0: return 0
        else: return e.thrust * c / e.Isp(planet, altitude)

    return (sum(e.thrust * c for (e,c) in engines.iteritems())
        /   sum(alpha(e,c)   for (e,c) in engines.iteritems()))

# combine 3x nuke and sail
# print combineIsp( [(engines[1], 1), (engines[2], 3)] )



###########################################################################
## Notes about experiments:
## I ran an experiment to see how Isp varied for the nuke as altitude
## changes.  Here are the results:

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
