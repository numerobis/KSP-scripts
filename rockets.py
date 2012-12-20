from __future__ import division # / means float div always

import math
from numbers import Number
import heapq
from LinkedList import LinkedList, cons, nil

###########################################################################
## Set these parameters, then run the script (a poor man's argument
## processing).

# Where do we start in Kerbin's atmosphere.  0 means the ground, n means n
# meters high, None means we start in space.
# We assume we start there along an optimal trajectory.
# If you can ride a jet to about 18km, you will save nearly 2500 m/s.
startAltitude = 18000

# How many burns do you need after circularizing, and how much acceleration
# in m/s^2 (aka thrust:mass ratio) do you need for each burn?  Here we go to
# the Mun and land.  We want to maintain 1m/s^2 for the tran-mun injection
# burn, and 3m/s^2 (nearly twice moon gravity) for the landing.  Then we do
# it again in reverse.  Don't worry about staging; the script (tries to)
# handle that.
# Each burn also has a payload.  The last payload has what we're left with when
# all has been consumed; if you land and don't mind leaving the lander legs behind,
# you can give that landing burn an extra payload.
class burn(object):
    def __init__(self, deltaV, accel, name, payload = 0):
        self.deltaV = deltaV
        self.accel = accel
        self.name = name
        self.payload = payload

munreturnburns = [
    burn(842, 1, "TMI"),
    burn(1000, 3, "mun land", payload = 0.5), # lander legs, ladder, etc
    burn(1000, 3, "mun depart"),
    burn(842, 1, "TKI", payload = 4.5) # capsule, 2 parachutes, RTG
]

# This is for a return from Laythe starting from Laythe orbit.
#startAltitude = None
#burns = [
#    burn(1500, 7.85*2, "Laythe liftoff"),
#    burn(1000, 1, "Laythe escape"),
#    burn(1200, 1, "Trans-Kerbin"),
#    burn(250, 1, "Kerbin reentry", payload = 12.63)
#]
# Answer: about 46 tonnes.

startAltitude = 0
laythereturnburns = [
    burn(1915, 1, "Trans-Jool"),
    # Payload here is left on Laythe: 4x hab modules, legs, chutes, instruments.
    burn(1500, 1, "Laythe intersect and land", 14.75),
    burn(2500, 7.85*2, "Laythe liftoff"),
    burn(1000, 1, "Laythe escape"),
    burn(1200, 1, "Trans-Kerbin"),
    # Payload here is what we return to Kerbin with: capsule, 2x modules,
    # panels, legs, chutes, instruments, ...
    burn(250, 1, "Kerbin reentry", payload = 10.25)
]
burns = munreturnburns

# We must specify the kind of symmetry we want our spacecraft to have.
# 2 is almost always optimal, but it may end up with a spacecraft that is too
# fat to be built in the VAB.  This can be a list; we'll try them all.
symmetry = 2

# Do the engines come below the final payload, or must they be side-mounted?
# 1 means to put the engines below; if side-mounted, specify how many sides.
# It's perfectly valid to have payload, 2 side-mounted engines, and thenceforth
# 4-fold symmetry, etc.
numBaseTowers = 2

# The automatic stage splitter will try splitting stages that achieve more than
# this much deltaV.  Smaller values give better rockets but longer
# optimization time.
minStageDeltaV = 750

# The search will ignore solutions that can only improve the current best by 1%
# or less (i.e. the mass of a new solution must be less than 0.99 times the old
# solution or else we don't really care).  You can change that; set to 1 to
# find the strict optimum given how the search works (which isn't actually
# optimal).  Set lower to allow a progressively less optimal solution, but
# search faster.
improvementRatio = 0.99



###########################################################################
# Standard gravity.
g0 = 9.80665


# Atmospheric pressure at a body is k e^{altitude/scale}
# Kerbin is 1*e^{altitude/5000}, Duna is 0.2e^{altitude/3000}, etc.
class planet(object):
    # F_{drag} = P v^2 m gamma, where P is atmospheric pressure, v is
    # the speed, m the mass (as a standin for cross-section).  gamma
    # combines a bunch of coefficients in one (including the 1/2 that would
    # normally be there).  It might not actually be a
    # constant, but it varies little over reported terminal velocities.
    gamma = 0.001

    # Following pretty much an optimal climb slope to a 100km orbit at Kerbin,
    # here's a correspondence between deltaV and altitude.  I haven't figured
    # out how to calculate this correctly for other planets... TODO
    climbSlope = [
        # deltaV, altitude
        (287, 923),
        (509, 2117),
        (731, 3461),
        (957, 5044),
        (1188, 6900),
        (1437, 9169),
        (1723, 12138),
        (2102, 16248),
        (2729, 22722),
        (3069, 32500)
    ]

    def __init__(self, name, gravity, radius, base, scale):
        self.name = name
        self.base = base
        self.scale = scale
        self.radius = radius # in km
        self._gravity = gravity

    def gravity(self, altitude = 0):
        """
        Return the gravitational acceleration at a given altitude above the
        surface.
        """
        if altitude == 0:
            return self._gravity
        else:
            baseradius = self.radius
            altradius = baseradius + (altitude / 1000)
            fraction = baseradius / altradius
            return self._gravity * fraction * fraction

    def terminalVelocity(self, altitude):
        """
        Return the terminal velocity at a given altitude.
        This is the speed at which drag is as strong as gravity.
        """

        # v = sqrt(g/ (alpha * gamma * e^{-altitude/beta}))
        # Then pull the e term out, and remember that gravity changes
        # (slightly) with altitude.
        return math.exp(0.5 * altitude / beta) * math.sqrt(
                    self.gravity(altitude) / (gamma * self.base) )

    def optimalDeltaV(self, altitude):
        """
        Return the deltaV needed to get to the given altitude on an optimal
        glide slop.  This means we have achieved terminal velocity.

        If the altitude is in orbit, return 4700 (the 100km orbit
        altitude).
        """
        if altitude is None: return 4700
        for (i, dataPoint) in enumerate(self.climbSlope):
            if dataPoint[1] == altitude:
                return dataPoint[0]
            elif dataPoint[1] > altitude:
                # linearly (!) interpolate
                if i == 0:
                    alpha = altitude / dataPoint[1]
                    return alpha * dataPoint[0]
                else:
                    prevPoint = self.climbSlope[i-1]
                    alpha = ((altitude - prevPoint[1]) /
                            (dataPoint[1] - prevPoint[1]))
                    return (alpha * dataPoint[0] +
                        (1-alpha) * prevPoint[0])
        return None

    def optimalAltitude(self, deltaV):
        """
        Return the altitude after burning deltaV on an optimal climb slope.
        This means we have achieved terminal velocity, and have burned an equal
        amount on gravity and drag.

        Return None if deltaV puts us outside the atmosphere.
        """
        # TODO: I don't understand the math yet, so just look it up.
        for (i, dataPoint) in enumerate(self.climbSlope):
            if dataPoint[0] == deltaV: return dataPoint[1]
            elif dataPoint[0] > deltaV:
                # linearly (!) interpolate
                if i == 0:
                    alpha = deltaV / dataPoint[0]
                    return alpha * dataPoint[1]
                else:
                    prevPoint = self.climbSlope[i-1]
                    alpha = ((deltaV - prevPoint[0]) /
                            (dataPoint[0] - prevPoint[0]))
                    return (alpha * dataPoint[1] +
                        (1-alpha) * prevPoint[1])
        return None

    def pressure(self, altitude):
        """
        Return the atmospheric pressure in Atm, which is a factor for
        interpolating the Isp of engines.
        """
        if altitude is None: p = 0
        else: p = self.base * math.exp(-altitude / self.scale)
        #print ("pressure at %s: %g" % (altitude, p))
        return p

kerbin = planet("kerbin", g0, 600, 1, 5000)

##############################
## Available engines
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

    def Isp(self, altitude):
        # Assumption: Isp is in a linear correspondence with pressure.
        # False for the turbojet.
        pressure = kerbin.pressure(altitude)
        Isp = pressure * self.IspAtm + (1.0 - pressure) * self.IspVac
        return Isp

    def __str__(self):
        return self.name

# We have a none engine for fuel-only stages in asparagus staging.  It has no
# mass, no thrust, Isp doesn't matter.  Make it radial, so that we can have
# radial engines above and we needn't add towers.
g_noEngine = engine("none", 0, 0, 0, 0, radial=True)

g_engines = [
    # One of the choices of engines is to have none...
    g_noEngine,

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
]

# To help the heuristics, choose the best possible Isp at a given altitude.
def maxIsp(altitude):
    ispMaxEngine = max(g_engines, key = lambda x: x.Isp(altitude))
    return ispMaxEngine.Isp(altitude)

# To help the heuristics, choose the best possible mass to achieve a given
# thrust.
maxThrustPerMassEngine = max(g_engines, key =
        lambda x: 0 if x.thrust == 0 else x.thrust / x.mass)
def lightestEngineForThrust(thrust):
    num = thrust / maxThrustPerMassEngine.thrust
    return (maxThrustPerMassEngine, num)

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

def combineIsp(engines, altitude):
    """
    Given a dictionary mapping engine -> count, compute the
    Isp of the system at the given altitude.

    This is the weighted sum of the impulses of each type, with weights
    based off the relative mass flow rate of the engines.
    """
    # the mass flow rate of an engine is thrust / (Isp * g0)
    # (unless there's no thrust, then it's 0 even if Isp is 0).
    # so the mass is sum(num_i * thrust_i / (Isp_i * g0))
    # the relative mass of engine type j is
    #   [num_j*thrust_j / (Isp_j*g0)] / sum(...)
    # and we cancel out the g0.
    def numerator(engcount):
        (engine, count) = engcount
        totalThrust = engine.thrust * count
        if totalThrust == 0: return 0
        return count * engine.thrust / engine.Isp(altitude)
    flows = [ numerator(x) for x in engines.iteritems() ]
    total = sum(flows)
    if total == 0: return 0
    alpha = [ f / total for f in flows ]
    weightedIsp = [
        engine.Isp(altitude) * a
            for (engine, a) in zip(engines.iterkeys(), alpha)
    ]
    return sum(weightedIsp)

# combine 3x nuke and sail
# print combineIsp( [(engines[1], 1), (engines[2], 3)] )


##############################
#
class WeakEngineException(Exception):
    def __init__(self, Isp): self.Isp = Isp
class MoarBoosters(Exception): pass

class stage(object):
    """
    """
    def __init__(self, deltaV, payload, engineType, nEngines, laterStages,
                numTowers, altitude, propMassOverride = None):

        # Calculate the effective Isp including later stages.
        allEngines = { engineType : nEngines }
        thrust = engineType.thrust * nEngines
        if engineType.vectoring:
            vectoringThrust = engineType.thrust * nEngines
        else:
            vectoringThrust = 0

        if not laterStages:
            Isp = engineType.Isp(altitude)
        else:
            self.collectUsableEngines(laterStages, allEngines)
            Isp = combineIsp(allEngines, altitude)
            thrust += laterStages.head.thrust
            vectoringThrust += laterStages.head.vectoringThrust

        # Calculate masses.
        # Engine mass is specified already.  We only count the engines
        # being dumped in this stage.
        # We add 0.05 for decouplers, struts, fuel lines.
        # Dry mass is payload, engines, and empty tanks.
        engineMass = engineType.mass * nEngines
        decouplerMass = self._decouplerConstant * numTowers
        dryMassNoTanks = payload + engineMass + decouplerMass

        # Get the propellant mass, and distribute it over the towers.
        #
        # TODO: handle SRBs.
        if propMassOverride:
            (propMass, tankMass) = propMassOverride
        else:
            (propMass, tankMass) = burnMass(deltaV, Isp, dryMassNoTanks)

        # Round up to fit an integer number of tanks on each tower.
        #
        # By the rocket equation, bringing extra unburned mass is OK, we'll
        # just finish our burn without finishing our fuel.  We can dump the
        # excess, or use it for the next stage burn.  So we round up to an
        # integer fuel mass per tower.
        propMass = math.ceil(propMass / numTowers) * numTowers
        tankMass = propMass / beta

        dryMass = dryMassNoTanks + tankMass
        fullMass = dryMass + propMass

        # Store a bunch of data (do we really need it all?)
        self.targetDeltaV = deltaV
        self.payload = payload
        self.engineType = engineType
        self.numEngines = nEngines
        self.numTowers = numTowers
        self.asparagus = bool(laterStages)
        self.engineMass = engineMass
        self.decouplerMass = decouplerMass
        self.propellantMass = propMass
        self.tankMass = tankMass
        self.dryMass = dryMass
        self.fullMass = fullMass
        self.Isp = Isp
        self.burnTime = burnTime(deltaV, Isp, thrust, dryMass)
        self.thrust = thrust
        self.vectoringThrust = vectoringThrust
        self.altitude = altitude

    def achievedDeltaV(self):
        return self.Isp * g0 * math.log(self.fullMass / self.dryMass)

    def acceleration(self):
        return self.thrust / self.fullMass

    def collectEngines(self, dict):
        if self.engineType in dict:
            dict[self.engineType] += self.numEngines
        else:
            dict[self.engineType] = self.numEngines

    @staticmethod
    def collectUsableEngines(stages, engineDict = None):
        """
        Collect all the engines usable in the remaining stages,
        listed from bottom to top.
        """
        if engineDict is None:
            engineDict = dict()
        for s in stages:
            s.collectEngines(engineDict)
            if not s.asparagus: break
        return engineDict

    # decoupler mass is assumed to be constant, no matter the circumstances, at
    # a mass of 0.05 (per tower that needs a decoupler)
    _decouplerConstant = 0.05

    def __str__(self):
        if self.engineType.name == "none":
            description = ("%dx %d T fuel" %
                (self.numTowers, self.propellantMass / self.numTowers))
        elif self.numEngines // self.numTowers == 1:
            description = ("%dx %s and %d T fuel" %
                (self.numTowers, self.engineType.name,
                 self.propellantMass // self.numTowers))
        else:
            description = ("%dx %d %s and %d T fuel" %
                (self.numTowers, self.numEngines // self.numTowers,
                 self.engineType.name, self.propellantMass // self.numTowers))
        return (
            "%g T:%s %s, %.2fs burn at %g kN (%.2f m/s^2), Isp %d"
            %   (self.fullMass,
                 " asparagus" if self.asparagus else "",
                 description,
                 self.burnTime,
                 self.thrust,
                 self.thrust / self.fullMass,
                 self.Isp)
        )

class burnRequirement(object):
    """
    Defines requirements for a burn.
    """
    def __init__(self, burnTime = None, acceleration = None):
        self.burnTime = burnTime
        self.acceleration = acceleration

    def satisfiedBy(self, stage):
        """
        Given a proposed stage, does it satisfy the burn requirements?
        """
        if self.burnTime and stage.burnTime > self.burnTime: return False
        if self.acceleration:
            # thrust is in kN, mass in tonnes, so the ratio is m/s^2
            accel = stage.thrust / stage.fullMass
            if accel < self.acceleration: return False
        # Nothing failed => success.
        return True

    def __str__(self):
        if self.burnTime is None and self.acceleration is None:
            return "unlimited burn"
        elif self.burnTime is None:
                return ("min %g m/s^2" % self.acceleration)
        elif self.acceleration is None:
            return ("max %gs burn" % self.burnTime)
        else:
            return ("min %g m/s^2, max %gs burn"
                        % (self.acceleration, self.burnTime))

def suggestEngineNumbers(symmetry, numBaseTowers, deltaV, payload, engineType, altitude,
              burnRequirement = None, laterStages = []):
    """
    Given an engine type, a maximum burn time,
    deltaV, and payload, return the stage description.
    This is where we choose the number of engines (the type is given).
    We maintain the desired symmetry.
    If we are asparagus staging, specify the list of later stages we can
    use.

    Returns a list of stages.  Each item is a choice for number of engines, and
    the required fuel, subject to symmetry requirements.

    If the engine is too wimpy, return an empty list.

    We are required to have at least 25% of our thrust be vectoring at all
    times.  If we're violating that, return an empty list.
    """
    # If we aren't asparagus staging (or this is the first stage), we need
    # a vectoring engine.
    if not laterStages and not engineType.vectoring:
        return []

    # How many towers do we have?
    # If we only have no engines, or only radial engines above, we can fit
    # on the base number of towers.  Otherwise we have a number of towers
    # according to symmetry, attached to the side of the previous stage.  Then,
    # further we can add any number of towers to this stage, according to
    # symmetry.  I'm allowing up to 2 more steps.
    for e in stage.collectUsableEngines(laterStages).iterkeys():
        if not e.radial:
            numBaseTowers = symmetry
            break
    numTowerChoices = [ numBaseTowers + symmetry * i for i in range(3) ]

    def tryNumTowers(numTowers):
        """
        Given a specified number of towers, try to make a stage with as few
        engines as possible.  Return None if we can't fit enough engines or
        the engine type is too weak, etc.
        """
        # Count up how many engines we might be able to use.
        # TODO: I'm not allowing making a tower with side towers.  That limits
        # a lot of things (particularly large loads).
        # TODO: I'm not allowing mixing engine types.  In particular, you could
        # have a stage with both standard and radial engines; or if you have
        # more than one tower, you could use a bicoupler or tricoupler with 2
        # types of engines.
        if engineType.large:
            # We can only fit one engine on each tower.
            numEngines = [ numTowers ]
            extraMass = [ 0 ]
        elif engineType.radial:
            # We can fit up to 8 Mark-55s, or 16 24-77s.
            # But I don't *want* to fit that many 24-77s, so set the max at 8.
            maxRadials = 8

            if numTowers == 1:
                # If we have one tower, we need to maintain symmetry.
                maxEngines = int(math.ceil(maxRadials/symmetry))
                numEngines = [ x * symmetry for x in range(1, maxEngines + 1) ]
            else:
                # Otherwise we can use any number, even a prime number.
                numEngines = [ x * numTowers for x in range(1,maxRadials + 1) ]
            extraMass = [ 0 for _ in range(len(numEngines)) ]
        else:
            # We can use 1, 2 (on a bicoupler), 3 (on a tricoupler), or 4 (on
            # chained bicouplers).  I suppose we could do crazy things too,
            # let's ignore that.
            numEngines = [ numTowers, 2*numTowers, 3*numTowers, 4*numTowers ]
            extraMass  = [ 0,
                0.1 * numTowers,  # bicoupler
                0.15 * numTowers, # tricoupler
                0.3 * numTowers ] # bicoupler with 2 bicouplers under it

        for (n, xmass) in zip(numEngines, extraMass):
            try:
                s = stage(deltaV, payload + xmass, engineType, n, laterStages,
                          numBaseTowers, altitude)
            except WeakEngineException, e:
                # Our Isp is too low to get anywhere.  Check if adding more
                # engines will improve the Isp.
                if e.Isp < engineType.Isp(altitude):
                    continue # Adding more may help.
                else:
                    return None # Adding more will get us nowhere.

            if s.vectoringThrust < 0.25 * s.thrust:
                if engineType.vectoring: continue
                else:
                    # We need more engines, but we also need to vector; try
                    # another type.
                    return None

            if not burnRequirement.satisfiedBy(s):
                # Not enough thrust, add more thrust.
                continue

            # If we get here, we like the stage.
            return s

        # If we get here, we can't fit enough engines on the specified number
        # of towers.
        return None

    # For each choice of number of towers, try that number of towers.  Return
    # only the first one we find.  We are thereby dropping some possibilities,
    # but otherwise the branching factor is just too ridiculous.
    for n in numTowerChoices:
        s = tryNumTowers(n)
        if s is not None:
            return [s]
    return []


def designStage(symmetry, numBaseTowers, deltaV, payload, altitude,
                laterStages = [],
                burnRequirement = None):
    """
    Given a payload mass (i.e. the mass of the next stage), calculate the
    type of engine, the number of engines needed if we have a limited burn
    time, and the propellant mass.

    We are required to have at least 25% of our thrust be vectoring at all
    times.

    Limitation: we only think of adding one type of engine, no mixing.

    deltaV in m/s
    payload in tonnes
    burnRequirement optional, an instance of a 'burn'
    laterStages: list of later stages whose thrust we can use in asparagus
        staging.

    Returns a list of possibilities in arbitrary order.
    """
    choices = []
    for engine in g_engines:
        choices.extend(
            suggestEngineNumbers(symmetry, numBaseTowers, deltaV, payload,
                    engine, altitude, burnRequirement, laterStages)
        )
    return choices

##############################
#
# We have n stages.
# We need a given total deltaV (in m/s), or a list of deltaV, one per
# stage.
# We have a given payload (in tonnes).
# We never fall below a given thrust:mass ratio.
# We can only use the engines in the map above.
# Tanks have arbitrary size (falsehood!), carry 4:1 prop:tank mass.
# Each burn must occur with the specified burn time, if specified.
#   None means we have unlimited burn time.
#   A number means to use that burn time for each stage.
#   A list specifies the burn time per stage (each one is a number or
#   None).
# We can use asparagus staging.  A false value means we don't, a true
#   value means we do.  If it's true, then at each stage we use all the
#   engines above.  Limitation: each stage only uses one type of engine,
#   though different stages can use different engines.
#
# Output the number and type of engines, and the
# number of tons of propellant at each stage.

class burnProfile(object):
    """
    Represent a set of burns to be done:
    * The number of burns.
    * The deltaV of each burn.
    * The altitude we'll be at for each burn.
    * The requirements (time or acceleration) for each part of the burns.
    * Payload each stage must carry (in addition to carrying the next stage).
    * To help the search, the best Isp available at each altitude.

    The burns come in order from top stage down.

    Use burnRequirements, deltaV, and altitude.
    """
    def __init__(self, burns, startAltitude = 0):
        # Important: store everything as tuples, not lists, because we need to
        # hash.
        burnRequirements = tuple(
            burnRequirement(acceleration = b.accel) for b in burns
        )
        deltaV  = tuple( b.deltaV for b in burns )
        payload = tuple( b.payload for b in burns )

        # Figure out the altitude -- but from the bottom up, not the top down.
        # We therefore need to reverse twice to get things back the way they
        # need to be.
        altitude = []
        lastDV = kerbin.optimalDeltaV(startAltitude)
        for dV in reversed(deltaV):
            altitude.append(kerbin.optimalAltitude(lastDV))
            lastDV = lastDV + dV
        altitude = tuple(reversed(altitude))

        # For each altitude, figure out what is the best Isp we can possible
        # achieve.
        maxIspEngine = [
            max(g_engines, key = lambda x: x.Isp(alt)) for alt in altitude
        ]
        maxIsp = tuple( x.Isp(alt) for (x, alt) in zip(maxIspEngine, altitude) )

        self.startAltitude = startAltitude
        self.deltaV           = deltaV
        self.altitude         = altitude
        self.burnRequirements = burnRequirements
        self.payload          = payload
        self.maxIsp           = maxIsp
        print ("Created new profile %s" % self)
        assert len(deltaV) > 0
        assert len(deltaV) == len(altitude) == len(burnRequirements) == len(payload) == len(maxIsp)

    def __hash__(self):
        return hash( (self.startAltitude, self.deltaV, self.burnRequirements, self.payload) )

    def __eq__(self, other):
        return (self.deltaV == other.deltaV
            and self.burnRequirements == other.burnRequirements
            and self.payloads == other.payloads
            and self.startAltitude == other.startAltitude)

    def __str__(self):
        # Not really intended for human consumption.
        prologue = ("%d stages" % len(self.deltaV))
        stages = [
            ("deltaV %d, altitude %s, require %s, %g payload, Isp is at best %d" % stage)
            for (stage) in zip(self.deltaV, self.altitude,
                    self.burnRequirements, self.payload, self.maxIsp) ]
        stages.insert(0, prologue)
        return "\n\t".join(stages)


class partialSolution(object):
    """
    Set up a partial solution with the given upper stages already
    selected.  Keep track of the required symmetry.
    If this is the first stage, set the payload.
    """
    def __init__(self, profile, stages, symmetry = None, numBaseTowers = None):
        self.profile  = profile
        self.stages   = stages
        self.symmetry = symmetry
        self.numBaseTowers = numBaseTowers
        self.complete = profile is None or (len(stages) == len(profile.deltaV))
        if symmetry is None or numBaseTowers is None:
            assert self.complete
            assert stages

        self.currentMass = stages.head.fullMass if stages else 0
        if self.complete:
            self.bestMass = self.currentMass
        else:
            bestMass = self.currentMass
            decouplers = stage._decouplerConstant
            numTowers = stages.head.numTowers if stages else numBaseTowers
            allEngines = stage.collectUsableEngines(stages)

            # For all lower stages, lower-bound the mass they will need.
            for i in xrange(len(stages), len(profile.deltaV)):
                (bestMass, allEngines) = self._lowerBound(
                    decouplers, numTowers, allEngines, bestMass, i)
            self.bestMass = bestMass

    def _lowerBound(self, decouplers, numTowers, allEngines, mass, i):
        """
        Lower bound the mass we'll need at stage i (where 0 is the top
        stage), assuming we need the given amount of mass at stage i-1.

        Improving the heuristic has a huge effect on runtime and memory
        use.
        """
        deltaV = self.profile.deltaV[i]
        Isp    = self.profile.maxIsp[i]
        req    = self.profile.burnRequirements[i]
        payload= self.profile.payload[i]
        allEngines = dict(allEngines)

        # We will need decouplers and other struts.
        mass += decouplers

        # We will need to carry the payload.
        mass += payload

        for _ in xrange(4):
            # Iterate:
            # 1. Check how much propellant we need, assuming the best
            #    possible Isp.
            # 2. Check how much thrust we need to push that propellant.
            # 3. Use the lightest engines to achieve it, if we don't
            #    already.
            # 4.

            # Assuming we achieve the best possible Isp for the given altitude,
            # compute the propellant use.  Don't add it yet!
            (bestProp, bestTank) = burnMass(deltaV, Isp, mass)

            # We are stuck using integer tonnage.  We might actually need more,
            # if we asparagus-stage.  But we might not, if we don't.
            # TODO: try both cases, take the better one.
            bestProp = math.ceil(bestProp)
            bestTank = bestProp / beta

            # Check how much thrust we need to push that mass.
            dryMass = mass + bestTank
            wetMass = mass + bestTank + bestProp
            if not req.burnTime: minTimeThrust = 0
            else:
                minTimeThrust = minThrustForBurnTime(deltaV, Isp,
                    dryMass, req.burnTime)
            if not req.acceleration: minAccelThrust = 0
            else:
                # F = ma
                minAccelThrust = wetMass * req.acceleration
            minThrust = max(minTimeThrust, minAccelThrust)
            # How much additional thrust do we need?
            curThrust = sum(e.thrust * n for (e, n) in allEngines.iteritems())
            if curThrust >= minThrust:
                # We have all the thrust we need to push the
                # propellant, we're done.
                break
            else:
                # Use the lightest possible engine.  It has worse Isp than
                # assumed; no matter, we're lower bounding.
                needThrust = minThrust - curThrust
                (engine, n) = lightestEngineForThrust(needThrust)
                if engine in allEngines:
                    allEngines[engine] += n
                else:
                    allEngines[engine] = n
                # Add in the mass, and iterate -- we'll need more
                # propellant now.
                mass += n * engine.mass

        # We got the engine mass all set up.  Now add in the final
        # propellant mass.
        mass += bestProp
        mass += bestTank

        return (mass, allEngines)



    def __lt__(self, other):
        if other is None: return True
        elif isinstance(other, Number): return self.bestMass < other
        else: return self.bestMass < other.bestMass

    def extend(self):
        """
        Return a set of options to add the next stage down on this
        rocket.  The options are sorted by mass.
        """
        assert not self.complete
        # i is the next burn to do
        i = len(self.stages)
        assert i < len(self.profile.deltaV)
        options = designStage(self.symmetry, self.numBaseTowers,
            self.profile.deltaV[i],
            self.currentMass + self.profile.payload[i],
            self.profile.altitude[i],
            burnRequirement = self.profile.burnRequirements[i],
            laterStages = self.stages) # asparagus staging
        options.extend( designStage(self.symmetry, self.numBaseTowers,
            self.profile.deltaV[i],
            self.currentMass + self.profile.payload[i],
            self.profile.altitude[i],
            burnRequirement = self.profile.burnRequirements[i],
            laterStages = nil) # standard staging
        )
        list.sort(options, key = lambda x: x.fullMass) # critical!
        solutions = []
        for stage in options:
            nextstages = cons(stage, self.stages)
            partial = partialSolution(self.profile, nextstages,
                        self.symmetry, self.numBaseTowers)
            solutions.append(partial)
        return solutions

    def __str__(self):
        return "\n".join( "\t" + str(s) for s in self.stages )


def designRocket(profiles, massToBeat = None,
        analyst = None, symmetries = 2, numBaseTowers = 1):
    if isinstance(symmetries, Number): symmetries = (symmetries,)
    if isinstance(numBaseTowers, Number): numBaseTowers = (numBaseTowers,)
    assert len(symmetries) == len(numBaseTowers)

    # bestKnown is None, a number in tonnes, or the best actual solution we've found
    bestKnown = massToBeat

    # trees is roughly speaking the list of profiles we're still pursuing
    trees = [ ]

    # Should we prune?  If the best known is an actual solution, we want pursue
    # a candidate if it might reduce the mass by at least 1%.  Less reduction,
    # we don't really care.

    def shouldKeep(candidate):
        if isinstance(bestKnown, partialSolution):
            if candidate.bestMass <= improvementRatio * bestKnown.bestMass:
                return True
        else:
            return candidate < bestKnown

    def generateSolutions(candidate):
        """
        Return a generator over all the complete and partial expansions of the
        candidate, pruning if the candidate isn't better than the best known
        solution.
        """
        if not shouldKeep(candidate):
            return

        # Generate all the children of this solution, i.e. all the
        # possibilities of adding an engine to the present solution.
        nextStageCandidates = candidate.extend()

        # There might not be any children (i.e. no engine can handle the next stage)
        if not nextStageCandidates:
            return

        if nextStageCandidates[0].complete:
            # Base case: return the first proposed solution -- they're sorted
            # by mass already, so either the first is an improvement or none of
            # them is.
            yield nextStageCandidates[0]
        else:
            # This is the general case: if the next stage is not bottom-most,
            # iterate over the children and generate their children in turn.
            for child in nextStageCandidates:
                # Recur first, then yield the child.
                # It's important to yield the child: otherwise, the iteration
                # blocks until we find a better solution.  This way, we can let
                # the caller stop searching in this tree.
                for x in generateSolutions(child):
                    yield x
                yield child

    try:
        def makeRoots(profiles):
            candidates = [ partialSolution(profile, nil, symmetry, numBase)
                for profile in profiles
                for (symmetry, numBase) in zip(symmetries, numBaseTowers) ]

            roots = [ generateSolutions(candidate) for candidate in candidates ]
            return roots

        def process(generator):
            """
            Move this generator along by one iterate, and return a pair
            (newBest, finished, profiles) where:
            * newBest is usually None, but otherwise is a new solution better
                than the best known
            * finished is usually False, otherwise True if the generator
                reached its last iterate.
            * profiles is a new list of profiles to try that might bring improvement
            """
            # try to push it a bit further (might not work: we might prune everything,
            # or the last one we got might have been the last solution to search)
            try:
                candidate = generator.next()
            except StopIteration:
                return (None, True, [])

            # Check if we improved the best known solution
            if (not candidate.complete) or (not (candidate < bestKnown)):
                return (None, False, [])

            # Improvement!
            newBest = candidate
            if not analyst:
                print ("Improved solution: %s" % newBest)
                return (newBest, False, [])
            else:
                # See if we have any good ideas for potentially improved stagings.
                # Optimality is violated right here: we only suggest for solutions
                # that are better than any prior.  That means we might miss a non-
                # optimal solution could have been locally improved to
                # the global optimum.
                analysis = analyst.analyze(newBest.stages)
                (improvement, profiles) = analyst.suggest(newBest.stages, analysis)
                profiles = list(profiles)

                while improvement:
                    asPartial = partialSolution(None, LinkedList(improvement))
                    assert (asPartial < newBest)
                    newBest = asPartial
                    analysis = analyst.analyze(newBest.stages)
                    (improvement, subprofiles) = analyst.suggest(newBest.stages, analysis)
                    profiles.extend(subprofiles)

                analyst.prettyPrint(newBest.stages, analysis)
                return (newBest, False, profiles)


        roots = makeRoots(profiles)
        nexpansions = 0
        nexpansionsLastPrinted = 0
        while(len(roots)):
            # Make a pass over all the roots, expanding them one step further.
            # Any suggested new profiles will be appended, and therefore
            # processed later in the loop.
            # Don't iterate; use the indices, so it's understandable what is
            # going on as we change the roots during the loop.
            j = 0
            n = len(roots)
            for i in xrange(len(roots)):
                (newBest, done, newprofiles) = process(roots[i])
                roots.extend(makeRoots(newprofiles))
                nexpansions += 1
                if newBest:
                    bestKnown = newBest
                if not done:
                    # We keep this one by copying it to the last unused
                    # location (which is, usually, i)
                    roots[j] = roots[i]
                    j += 1
                # If we added this profile this round, pick at it a bunch more:
                # once per root we had before starting.  Otherwise it takes
                # forever and a day to do any searching around that new
                # part of the solution space.
                if i > n:
                    for _ in xrange(n):
                        (newBest, done, newprofiles) = process(roots[i])
                        roots.extend(makeRoots(newprofiles))
                        nexpansions += 1
                        if newBest:
                            bestKnown = newBest
                        if done: break
            del roots[j:]

            # let the user know things are moving along
            if nexpansions >= nexpansionsLastPrinted + 1000:
                nexpansionsLastPrinted = nexpansions
                print ("%d choices considered, %d avenues remain"
                        % (nexpansions, len(roots)) )

        print ("Completed search after %d evaluations" % nexpansions)
    except KeyboardInterrupt:
        print ("Cancelled search after %d evaluations" % nexpansions)

    if isinstance(bestKnown, partialSolution):
        return bestKnown.stages
    else:
        return None


###########################################################################
# Convert the burns specified at the top of the file to the burnProfile
# interface, and split up the original burns into chunks.
#
# Any burn of more than 1km/s can be split into equal parts; generate a profile
# for each split possibility.  This is a combinatorial number of profiles, and
# then each one splits into a combinatorial number of possibilities!
#
# In practice, splitting in 500m/s seems to just take too long on large
# missions and not actually improve much on small missions.
#
def splitBurns(burns):
    profiles = []
    # subburns are burn instances, but stored in a cons list, so they're in
    # top-first order rather than bottom-first.  This is conveniently how all
    # the internal code wants it.
    def recur(subburns, burnIdx):
        if burnIdx == len(burns):
            # Base case:
            # convert the subburns to a complete burn profile
            profile = burnProfile(subburns,
                startAltitude = startAltitude)
            profiles.append(profile)
        else:
            # Inductive case: for every split of burn i, recur.
            b = burns[burnIdx]
            maxstages = int(math.ceil(b.deltaV / minStageDeltaV)) # max number of splits
            for nstages in xrange(1, maxstages+1):
                # Recur with nstages stages.
                mysubburns = subburns
                dV = b.deltaV / nstages
                # Add n-1 of the burns, then the last sub-burn has the payload.
                if nstages > 1:
                    for i in xrange(nstages - 1):
                        subname = ("%s [%d]" % (b.name, i))
                        mysubburns = cons(burn(dV, b.accel, subname), mysubburns)
                i = nstages - 1
                subname = ("%s [%d]" % (b.name, i))
                mysubburns = cons(burn(dV, b.accel, subname, b.payload), mysubburns)
                recur(mysubburns, burnIdx + 1)

    recur(nil, 0)
    return profiles

# Local search trick, and pretty printing.
class analyst(object):
    def __init__(self, burns):
        """
        Set up, relative to the *original* burns that the user wants.
        Said burns are ordered bottom up.
        """
        self.burns = burns
        self.totalDeltaV = sum(b.deltaV for b in self.burns)

    class stageData(object):
        def __init__(self, burnIds, dumpDV, payload):
            self.burnIds = burnIds
            self.dumpDV = dumpDV    # delta-V capacity we can't use
            self.payload = payload

    def _deltaVMap(self, deltaVs):
        """
        For each deltaV, report:
        * list of burns that are completed by this deltaV
        * an incomplete burn that this deltaV helps with, or None
        Ignore acceleration constraints.

        DeltaVs come in the same order as the burns: bottom up.
        """
        burnIdx = 0
        b = self.burns[burnIdx]
        burndV = b.deltaV
        totaldV = 0
        data = []
        for dV in deltaVs:
            completedBurns = []
            partialBurn = None
            while burnIdx < len(self.burns) and dV > 0:
                # We have dV left over, and we can use it.
                if dV < burndV - 1e-3:
                    # We can make progress on a burn.
                    burndV -= dV
                    partialBurn = b
                    break
                else:
                    # We can complete a burn (perhaps missing by 1mm/s).
                    completedBurns.append(b)
                    dV -= burndV
                    burnIdx += 1
                    if burnIdx == len(self.burns): break
                    b = self.burns[burnIdx]
                    burndV = b.deltaV
            data.append( (completedBurns, partialBurn) )
        return data

    def analyze(self, stages):
        """
        Returns a list of burns per stage.
        If the list is empty, the stage is wasted (for fuel; its thrust may
        still be important)!
        The list will be passed to suggest and to prettyPrint later.

        Stages come in order from bottom up.
        """
        burnIdx = 0
        b = self.burns[burnIdx]
        burndV = b.deltaV
        totaldV = 0
        data = []
        for (stageIdx, s) in enumerate(stages):
            dV = s.achievedDeltaV()
            startdV = totaldV
            totaldV += dV

            # Keep track of all burns affected by this stage.
            ids = [ ]
            mustDump = False
            payload = 0
            while burnIdx < len(self.burns) and dV > 0:
                if s.acceleration() < b.accel:
                    # We have dV left over, but we don't accelerate enough.
                    # TODO: actually the acceleration improves over time,
                    # so it might be enough now even though it wasn't at the
                    # start of the stage.
                    mustDump = dV
                    break
                else:
                    # We have dV left over, and we can use it.
                    ids.append(burnIdx)
                    if dV < burndV:
                        # We can make progress on a burn.
                        burndV -= dV
                        break
                    else:
                        # We can complete a burn.
                        dV -= burndV
                        payload += self.burns[burnIdx].payload
                        burnIdx += 1
                        if burnIdx == len(self.burns):
                            break
                        b = self.burns[burnIdx]
                        burndV = b.deltaV
            if burnIdx == len(self.burns) and dV > 0:
                # Wasted excess at the end of the flight, which we usually
                # want to dump instead of landing on!
                mustDump = dV

            data.append(self.stageData(ids, mustDump, payload))
        return data


    def suggest(self, stages, analysis):
        """
        Suggest new profiles based on the solution we have.
        I'm still struggling to get a non-buggy local improvement.

        * From the bottom up, look at each stage.  Produce a new profile
          that asks for the deltaV the stages provide, up to the total
          deltaV, and the max acceleration of any burn affected by this
          burn.  In other words, if we have wasted stages at the end of our
          trip, we will eliminate them.
        * From the top down, same thing.  In other words, any wasted deltaV at
          the end is pushed back, and hopefully we can squeeze some of it all
          the way back.

        Any other bright ideas?  Isp-weighted deltaV spreading?

        The stages are listed bottom-up.
        """
        # Local improvement.  Assign the best we find to this variable.
        improvedStages = None

        # The stages listed top-down, if that helps.
        revStages = tuple(reversed(tuple(stages)))
        revAnalysis = tuple(reversed(analysis))

        def suggestDeltaVs(stages):
            """
            Suggest deltaVs that don't overshoot, following the given stages
            (which are either bottom-up or top-down).
            """
            dVremaining = self.totalDeltaV
            deltaVs = []
            for s in stages:
                dV = s.achievedDeltaV()
                if dV == 0: continue
                if dV < dVremaining:
                    if dVremaining - dV < minStageDeltaV:
                        deltaVs.append(dVremaining)
                        break # after this, all stages are wasted
                    else:
                        deltaVs.append(dV)
                        dVremaining -= dV
                else:
                    deltaVs.append(dVremaining)
                    break # after this, all stages are wasted
            return deltaVs

        def makeProfile(deltaVs):
            """
            Return a burn profile that corresponds to achieving the given
            deltaVs, and matching up with acceleration requirements.
            """
            def makeBurn(deltaV, (fullBurns, partialBurn)):
                if partialBurn:
                    allBurns = [ b for b in fullBurns ]
                    allBurns.append(partialBurn)
                else:
                    allBurns = fullBurns
                assert len(allBurns) > 0

                # TODO: this is pessimistic.  Our acceleration will increase
                # over the stage.  However, we have no way to convey such
                # fine-grained information in the profile being optimized.
                accel = max( b.accel for b in allBurns )

                if len(allBurns) == 1:
                    name = allBurns[0].name
                else:
                    name = ",".join(b.name for b in allBurns)

                # Payloads that we have to lift, not counting later stages.
                # We don't care about the partial burn, which a later stage
                # will handle.
                payload = sum( b.payload for b in fullBurns ) if fullBurns else 0

                return burn(deltaV, accel, name, payload)

            # The deltaVs and burnsByStage are bottom-up, we need the burns top-down.
            burnsByStage = self._deltaVMap(deltaVs)
            newburns = [  makeBurn(*x) for x in reversed(zip(deltaVs, burnsByStage)) ]
            assert len(newburns) > 0
            return burnProfile(newburns, startAltitude = startAltitude)

        # suggest deltaV lists in both directions
        deltaVsUp = suggestDeltaVs(stages)
        deltaVsDown = tuple(reversed(suggestDeltaVs(revStages)))

        # Rebuild the stages, trying to reduce capacity to drop wasted deltaV.
        (lastStage, lastData) = (revStages[0], revAnalysis[0])
        newStages = nil
        for (s, d) in zip(revStages, revAnalysis):
            payload = (newStages.head.fullMass if newStages else 0) + d.payload
            if not d.burnIds:
                # This is a completely wasted stage.  Rebuild it with no fuel.
                # TODO: merge with the stage below if it has the same engine,
                # or no engine, or if we have no engine.
                newS = stage(0, payload, s.engineType, s.numEngines,
                    newStages if s.asparagus else None, s.numTowers, s.altitude)
            else:
                # This stage is useful.  However, its mass may be reduced.
                # How much dV we actually provide versus how much we can use
                # (neither of which are how much we asked for), and try to
                # trim our propellant use.
                targetDeltaV = s.achievedDeltaV() - d.dumpDV
                try:
                    newS = stage(targetDeltaV, payload, s.engineType, s.numEngines,
                        newStages if s.asparagus else None, s.numTowers, s.altitude)
                except WeakEngineException:
                    # This implies roundoff error, since we didn't change the engines.
                    newS = None
                if newS.propellantMass > s.propellantMass:
                    # This implies roundoff, since we reduced the mass above.
                    newS = None
                if not newS:
                    # I was hoping we'd get the same deltaV with less
                    # propellant, but roundoff bit me; just use the same
                    # propellant as before.
                    newS = stage(targetDeltaV, payload, s.engineType, s.numEngines,
                        newStages if s.asparagus else None, s.numTowers, s.altitude,
                        propMassOverride = (s.propellantMass, s.tankMass))

            newStages = cons(newS, newStages)
        if newStages.head.fullMass < stages.head.fullMass:
            print ("Locally improved from %g T to %g T" %
                    (stages.head.fullMass, newStages.head.fullMass))
            improvedStages = newStages

        # Return the improved set of stages, and also return the potentially
        # better profiles to search.
        return ( improvedStages, ( makeProfile(deltaVsUp), makeProfile(deltaVsDown) ) )


    def prettyPrint(self, stages, data):
        m0 = self.burns[-1].payload
        m1 = stages.head.fullMass
        Isp = self.totalDeltaV / (g0 * math.log(m1 / m0)) # log is ln
        print ("Solution in %d stages:" % (len(stages),))
        print ("  Final payload: %g T" % m0)
        print ("  Launch mass:   %g T" % m1)
        print ("  Ratio:         %g" % (m1/m0))
        print ("  Implied Isp:   %g" % Isp)
        # TODO: we're not counting waste when we're forced to dump fuel.
        print ("  Wasted deltaV: %g m/s" %
                (sum(s.achievedDeltaV() for s in stages) - self.totalDeltaV))
        if startAltitude == 0:
            print ("  Start on Kerbin surface")
        elif startAltitude is None:
            print ("  Start in the vacuum")
        else:
            print ("  Start at %dm altitude" % startAltitude)

        for (sIdx, (s, d)) in enumerate(zip(stages, data)):
            names = [ self.burns[burnId].name for burnId in d.burnIds ]
            if len(names) == 0:
                print ("  Stage %d: %g m/s, wasted stage!" % (sIdx, s.achievedDeltaV()))
            elif len(names) == 1:
                print ("  Stage %d: %g m/s, used for %s" % (sIdx, s.achievedDeltaV(), names[0]))
            else:
                print ("  Stage %d: %g m/s, used for %d burns: %s" %
                    (sIdx, s.achievedDeltaV(), len(names), ", ".join(names)))
            print ("    %s" % (s,))
            if d.dumpDV and sIdx != len(data) - 1:
                print ("\t* must dump %g m/s before next burn" % d.dumpDV)
            if d.payload:
                print ("\t* includes %g T payload" % d.payload)


startDeltaV = kerbin.optimalDeltaV(startAltitude)
atmosphericDeltaV = 4700 - startDeltaV
if atmosphericDeltaV > 0:
    burns.insert( 0, burn(atmosphericDeltaV, 2.2*g0, "Kerbin launch") )

profiles = splitBurns(burns)

# Design the rocket!
shrink = analyst(burns)
soln = designRocket(profiles,
            massToBeat = None, analyst = shrink, symmetries = symmetry,
            numBaseTowers = numBaseTowers)
if soln:
    data = shrink.analyze(soln)
    shrink.prettyPrint(soln, data)
else:
    print "You will not be going to space today."
