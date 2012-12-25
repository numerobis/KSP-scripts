from __future__ import division # / means float div always

import math
from numbers import Number
import heapq
from LinkedList import LinkedList, cons, nil

import ascent
import engine
import physics

##############################
#

class stage(object):
    """
    """
    def __init__(self, deltaV, payload, engineType, nEngines, laterEngines,
                numTowers, planet, altitude, propMassOverride = None):

        # Calculate the effective Isp including later stages.
        allEngines = dict(laterEngines.iteritems() if laterEngines else [])
        if engineType in allEngines:
            allEngines[engineType] += nEngines
        else:
            allEngines[engineType] = nEngines
        thrust = sum(e.thrust * count for (e, count) in allEngines.iteritems())
        vectoringThrust = sum(e.thrust * count for (e, count) in
            allEngines.iteritems() if e.vectoring)
        Isp = engine.combineIsp(allEngines, planet, altitude)

        # Calculate masses.
        # Engine mass is specified already.  We only count the engines
        # being dumped in this stage.
        # We add 0.05 per decoupler.
        # Dry mass is payload, engines, and empty tanks.
        engineMass = engineType.mass * nEngines
        decouplerMass = self._decouplerConstant * numTowers
        dryMassNoTanks = payload + engineMass + decouplerMass

        # Get the propellant mass, and distribute it over the towers.
        #
        # Round up to fit an integer number of tanks on each tower.
        #
        # By the rocket equation, bringing extra unburned mass is OK, we'll
        # just finish our burn without finishing our fuel.  We can dump the
        # excess, or use it for the next stage burn.  So we round up to an
        # integer fuel mass per tower.
        #
        # TODO: handle SRBs.
        #
        if propMassOverride:
            (propMass, tankMass) = propMassOverride
        else:
            (propMass, tankMass) = engine.burnMass(deltaV, Isp, dryMassNoTanks)
            propMass = math.ceil(propMass / numTowers) * numTowers
            tankMass = propMass / engine.beta

        dryMass = dryMassNoTanks + tankMass
        fullMass = dryMass + propMass

        # Store a bunch of data (do we really need it all?)
        self.targetDeltaV = deltaV
        self.payload = payload
        self.engineType = engineType
        self.numEngines = nEngines
        self.numTowers = numTowers
        self.asparagus = bool(laterEngines)
        self.engineMass = engineMass
        self.decouplerMass = decouplerMass
        self.propellantMass = propMass
        self.tankMass = tankMass
        self.dryMass = dryMass
        self.fullMass = fullMass
        self.Isp = Isp
        self.burnTime = engine.burnTime(deltaV, Isp, thrust, dryMass)
        self.thrust = thrust
        self.vectoringThrust = vectoringThrust
        self.altitude = altitude
        self.planet = planet

    def achievedDeltaV(self):
        return self.Isp * physics.g0 * math.log(self.fullMass / self.dryMass)

    def acceleration(self):
        return self.thrust / self.fullMass

    def collectEngines(self, d = None):
        if d is None: d = dict()
        if self.engineType in d:
            d[self.engineType] += self.numEngines
        else:
            d[self.engineType] = self.numEngines
        return d

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

def suggestEngineNumbers(symmetry, numBaseTowers, deltaV, payload, engineType,
        altitude = None, planet = None,
        laterEngines = dict(), acceleration = None):
    """
    Given an engine type, return the stage description.
    This is where we choose the number of engines (the type is given, but we
    have to achieve a certain acceleration).
    We maintain the desired symmetry.
    If we are asparagus staging, specify the list of later engines we can use.

    Returns a list of stages.  Each item is a choice for number of engines, and
    the required fuel, subject to symmetry requirements.

    If the engine is too wimpy, return an empty list.

    We are required to have at least 25% of our thrust be vectoring at all
    times.  If we're violating that, return an empty list.
    """
    # If we aren't asparagus staging (or this is the first stage), we need
    # a vectoring engine.
    if not laterEngines and not engineType.vectoring:
        return []

    # How many towers do we have?
    # If we only have no engines, or only radial engines above, we can fit
    # on the base number of towers.  Otherwise we have a number of towers
    # according to symmetry, attached to the side of the previous stage.  Then,
    # further we can add any number of towers to this stage, according to
    # symmetry.  I'm allowing up to 2 more steps.
    for e in laterEngines:
        if not e.radial:
            numBaseTowers = symmetry
            break
    if engineType == engine.noEngine:
        # don't try adding more of these to get more thrust!
        numTowerChoices = [ numBaseTowers ]
    else:
        numTowerChoices = [ numBaseTowers + symmetry * i for i in range(3) ]

    def tryNumTowers(numTowers, nEnginesAttempted):
        """
        Given a specified number of towers, try to make a stage with as few
        engines as possible.  Return None if we can't fit enough engines or
        the engine type is too weak, etc.
        """
        # Count up how many engines we might be able to use.
        # TODO: I'm not allowing mixing engine types.  In particular, you could
        # have a stage with both standard and radial engines; or if you have
        # more than one tower, you could use a bicoupler or tricoupler with 2
        # types of engines.
        if engineType == engine.noEngine:
            # The "none" engine is just an asparagus fuel stage.
            numEngines = [ 1 ]
            extraMass = [ 0 ]
        elif engineType.large:
            # We can only fit one engine on each tower.
            numEngines = [ numTowers ]
            extraMass = [ 0 ]
        elif engineType.radial:
            # We can fit up to 8 Mark-55s, or 16 24-77s.
            # But the part count gets ridiculous fast, so keep it much lower.
            maxRadials = 4

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
            if n in nEnginesAttempted: continue
            else: nEnginesAttempted.add(n)
            try:
                s = stage(deltaV, payload + xmass, engineType, n, laterEngines,
                          numBaseTowers, planet, altitude)
            except engine.WeakEngineException, e:
                # Our Isp is too low to get anywhere.  Check if adding more
                # engines will improve the Isp.
                if e.Isp < engineType.Isp(planet, altitude):
                    continue # Adding more may help.
                else:
                    return None # Adding more will get us nowhere.

            if s.vectoringThrust < 0.25 * s.thrust:
                if engineType.vectoring: continue
                else:
                    # We need more engines, but we also need to vector; try
                    # another type.
                    return None

            if s.thrust < acceleration * s.fullMass:
                # Not enough thrust, add more thrust.
                continue

            # If we get here, we like the stage.
            return s

        # If we get here, we can't fit enough engines on the specified number
        # of towers.
        return None

    # For each choice of number of towers, try that number of towers.
    # That's too ridiculous for radials, so instead, return just the first
    # number of radials that gives enough thrust.
    choices = []
    nEnginesAttempted = set()
    for n in numTowerChoices:
        s = tryNumTowers(n, nEnginesAttempted)
        if s is None: continue
        if engineType.radial:
            return [s]
        else:
            choices.append(s)

    return choices


def designStage(symmetry, numBaseTowers, deltaV, payload,
                altitude, planet,
                laterEngines = [],
                acceleration = None):
    """
    Given a payload mass (i.e. the mass of the next stage), calculate the
    type of engine, the number of engines needed if we have a limited burn
    time, and the propellant mass.

    We are required to have at least 25% of our thrust be vectoring at all
    times.

    Limitation: we only think of adding one type of engine, no mixing.

    deltaV in m/s
    payload in tonnes

    laterEngines: list of engines on upper stages that we can use in asparagus
        staging.

    Returns a list of possibilities in arbitrary order.
    """
    choices = []
    for eType in engine.types:
        choices.extend(
            suggestEngineNumbers(symmetry, numBaseTowers, deltaV, payload,
                    eType, altitude=altitude, planet=planet,
                    acceleration = acceleration, laterEngines = laterEngines)
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

class rawBurn(object):
    def __init__(self, name, deltaV, accel, payload, altitude = None, planet = None):
        self.name = name
        self.deltaV = deltaV
        self.acceleration = accel
        self.payload = payload
        if planet is None or altitude is None:
            self.altitude = None
            self.planet = None
        else:
            self.altitude = altitude
            self.planet = planet

    def __str__(self):
        return ("raw burn %s: %g m/s at %g m/s^2, payload %g%s"
            % (self.name, self.deltaV, self.acceleration, self.payload,
               ("%g m above %s" % (self.altitude, self.planet))
                    if self.planet else " in space"))

class deepSpaceBurn(object):
    def __init__(self, name, deltaV, accel = 1, payload = 0):
        self.name = name
        self.deltaV = deltaV
        self.acceleration = accel # m/s
        self.payload = payload

    def convert(self, n):
        """
        Return a list of up to n raw burns corresponding to this burn.
        """
        burns = []
        if n > 1:
            for i in range(n-1):
                burns.append(
                    rawBurn(self.name, self.deltaV / n, self.acceleration, 0)
                )
        burns.append(rawBurn(self.name, self.deltaV / n, self.acceleration,
                    self.payload))
        return burns

class liftoffBurn(object):
    def __init__(self, name, planet, orbit,
            acceleration = None, initialAltitude = 0, initialVelocity = None,
            payload = 0):
        """
        Define a sequence of burns for a given planet, one per altitude.
        orbit (m): the height of the orbit to achieve.  The orbit is
            circular, equatorial, and in the same direction as sidereal
            rotation.
        planet: by default Kerbin, other planets can be specified.
        acceleration: if specified, we ask the ascent module to assume this is
            the acceleration that will be used.  Otherwise, we let the ascent
            module use its default.
        initialAltitude: the ascent will start at this altitude.
        initialVelocity: this is an initial velocity vector that will magically
            happen; it applies to the velocity at the first altitude.  If a
            number, it is assumed to be vertical velocity.  If a 2-vector, it
            is relative to the surface tangent.
        payload (tonnes): The payload is extra payload that will stay in the
            given orbit, beyond whatever further burns are required.
        """
        self.name = name
        self.payload = payload
        self.planet = planet
        self.slope = ascent.climbSlope(planet, orbit,
            initialAltitude = initialAltitude,
            initialVelocity = initialVelocity,
            acceleration = acceleration)
        self.deltaV = self.slope.deltaV()
        if acceleration is not None:
            self.acceleration = acceleration
        else:
            # TODO: really we need acceleration to be by a given deltaV,
            # since the required acceleration changes over the climb (more
            # early and when the atmosphere starts to thin, less in between).
            self.acceleration = max(x.thrust for x in self.slope._climbSlope)


    def convert(self, n):
        """
        Return this burn split into n stages, bottom first.
        """
        # TODO: take account of Isp varying over pressure: instead of
        # splitting by equal deltaV we should split by equal Isp.
        burns = []
        accel = max(x.thrust for x in self.slope._climbSlope)
        if n > 1:
            for i in range(n-1):
                altitude = self.slope.altitudeAtDeltaV(self.deltaV * i / n, None)
                burns.append(
                    rawBurn(self.name, self.deltaV / n, accel, 0,
                        planet = self.planet, altitude = altitude)
                )
        lastAltitude = self.slope.altitudeAtDeltaV(self.deltaV * (n-1) / n, None)
        burns.append(rawBurn(self.name, self.deltaV / n, accel,
                    planet = self.planet, altitude = lastAltitude,
                    payload = self.payload))
        return burns





class burnProfile(object):
    """
    Represent a set of burns to be done:
    * The number of burns.
    * The deltaV of each burn.
    * If the burn is atmospheric, the altitude and the planet.
    * The acceleration requirement for each part of the burns.
    * Payload each stage must carry (in addition to carrying the next stage).
    * To help the search, the best Isp available at each altitude.

    The burns come in order from top stage down, and are items with a
    convert(n) function that converts to a rawBurn.
    """
    def __init__(self, rawburns):
        self.rawburns = tuple(rawburns)
        self.maxIsp = tuple( engine.maxIsp(b.planet, b.altitude if b.planet else None)
                for b in rawburns )

    def __str__(self):
        # Intended just for debugging...
        prologue = ("%d burns" % len(self.rawburns))
        burnstrs = [
            ("deltaV %d, %s, require %s m/s, %g payload, Isp is at best %d" %
                (b.deltaV,
                ("%gm altitude on %s" % (b.altitude, b.planet)
                    if b.planet else "in space"),
                b.acceleration, b.payload, isp))
            for (b, isp) in zip(self.rawburns, self.maxIsp)
        ]
        burnstrs.insert(0, prologue)
        return "\n\t".join(burnstrs)


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
        self.complete = profile is None or (len(stages) == len(profile.rawburns))
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
            try:
                for i in xrange(len(stages), len(profile.rawburns)):
                    (bestMass, allEngines) = self._lowerBound(
                        decouplers, numTowers, allEngines, bestMass, i)
                self.bestMass = bestMass
            except engine.WeakEngineException:
                # Totally impossible even with the best engines.  Set this to
                # basically infinite lower bound.
                self.bestMass = 1e30

    def __lt__(self, other):
        if other is None:
            return True
        elif isinstance(other, Number):
            return self.bestMass < other
        else:
            return self.bestMass < other.bestMass

    def remainingStages(self):
        if self.complete:
            # Locally-improved solutions have no profile, so the general
            # code would fail.
            return 0
        else:
            return len(self.profile.rawburns) - len(self.stages)

    def _lowerBound(self, decouplers, numTowers, allEngines, mass, i):
        """
        Lower bound the mass we'll need at stage i (where 0 is the top
        stage), assuming we need the given amount of mass at stage i-1.

        Improving the heuristic has a huge effect on runtime and memory
        use.
        """
        b          = self.profile.rawburns[i]
        deltaV     = b.deltaV
        Isp        = self.profile.maxIsp[i]
        accel      = b.acceleration
        payload    = b.payload
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
            # 3. Use the lightest engines to achieve it, which adds mass.
            # 4. goto 1.

            # Assuming we achieve the best possible Isp for the given altitude,
            # compute the propellant use.  Don't add it yet!
            (bestProp, bestTank) = engine.burnMass(deltaV, Isp, mass)

            # We are stuck using integer tonnage.  We might actually need more,
            # if we asparagus-stage.  But we might not, if we don't.
            # TODO: try both cases, take the better one.
            bestProp = math.ceil(bestProp)
            bestTank = bestProp / engine.beta

            # Check how much thrust we need to push that mass.
            dryMass = mass + bestTank
            wetMass = mass + bestTank + bestProp
            minThrust = wetMass * accel
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
                (e, n) = engine.lightestEngineForThrust(needThrust)
                if e in allEngines:
                    allEngines[e] += n
                else:
                    allEngines[e] = n
                # Add in the mass, and iterate -- we'll need more
                # propellant now.
                mass += n * e.mass

        # We got the engine mass all set up.  Now add in the final
        # propellant mass.
        mass += bestProp
        mass += bestTank

        return (mass, allEngines)

    def extend(self):
        """
        Return a set of options to add the next stage down on this
        rocket.  The options are sorted by mass.
        """
        assert not self.complete
        # i is the next burn to do
        i = len(self.stages)
        assert i < len(self.profile.rawburns)
        b = self.profile.rawburns[i]
        engines = stage.collectUsableEngines(self.stages)
        options = designStage(self.symmetry, self.numBaseTowers,
            b.deltaV, self.currentMass + b.payload, b.altitude, b.planet,
            acceleration = b.acceleration,
            laterEngines = engines) # asparagus staging
        options.extend( designStage(self.symmetry, self.numBaseTowers,
            b.deltaV, self.currentMass + b.payload, b.altitude, b.planet,
            acceleration = b.acceleration,
            laterEngines = dict()) )# standard staging: no later engines for use
        list.sort(options, key = lambda x: x.fullMass) # critical!
        solutions = []
        for s in options:
            nextstages = cons(s, self.stages)
            partial = partialSolution(self.profile, nextstages,
                        self.symmetry, self.numBaseTowers)
            solutions.append(partial)
        return solutions

    def __str__(self):
        strs = []
        if self.complete:
            strs.append("Complete candidate with mass %g, %d stages" % (self.bestMass, len(self.stages)))
        else:
            strs.append("Partial candidate with mass %g, %d stages" % (self.bestMass, len(self.stages)))
        strs.extend( "\t" + str(s) for s in self.stages )
        return "\n".join(strs)


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
    improvementRatio = 1
    def shouldKeep(candidate):
        if isinstance(bestKnown, partialSolution):
            if candidate.bestMass <= improvementRatio * bestKnown.bestMass:
                return True
        else:
            return candidate < bestKnown

    def greedySolution(partial):
        """
        Greedily extend the partial solution to completion.
        Return None if the greedy solution got pruned, or no solution exists.
        """
        while not partial.complete:
            if not shouldKeep(partial):
                return
            children = partial.extend()
            if not children: return
            partial = children[0]
        return partial

    def semiGreedySolutions(candidate, depth):
        """
        Returns a generator that allows iterating over greedy completions of
        all possibilities of engine choice for the top 'depth' stages.
        Depth 0 is purely greedy, whereas in an n-stage rocket depth n is
        completely optimal for the burn profile.
        This is support for implementing an iterative deepening search.

        The concept is that engine choice in low stages matters less to the
        overall mass than does engine mass in upper stages: e.g. it takes about
        30T to lift 1T to LKO.  Example: Using an engine one tonne too heavy
        when leaving LKO will thus cost us 30T, whereas an equally wrong engine
        at the first stage costs us 1T.  Depth-first does an exhaustive search
        for the best low-stage engine for a fixed upper-stage engine, which is
        precisely wrong.  Breadth-first with greedy extension does an
        exhaustive search of upper-stage engines and doesn't worry much about
        the low-stage engines.
        """
        if not shouldKeep(candidate):
            yield None
            return

        if depth == 0:
            yield greedySolution(candidate)
        else:
            for child in candidate.extend():
                for x in semiGreedySolutions(child, depth - 1):
                    yield x

    def generateSolutions(candidate):
        """
        Create a generator that allows iterating in over all the
        completions of the candidate.

        Yields None periodically if the search is being unfruitful, to allow
        trying other search trees in parallel.

        The search order is iterative deepening: all depth-0 solutions, then
        all depth-1 solutions, ...
        """
        # TODO: The algorithm repeats solutions: at depth 1, the first greedy
        # solution is exactly the depth 0 greedy solution.  At depth 2, the
        # first is again the depth 0 greedy solution, but also for each child
        # of the candidate, the first greedy solution is also the corresponding
        # depth-1 greedy solution.  So the first greedy solution of a leaf is
        # repeated up to n times on an n-stage rocket -- though pruning can
        # reduce that.  Might be nice to fix.

        # Note: If we have 1 stage remaining, the greedy solution is optimal,
        # so the upper bound on the range is correct.
        for depth in xrange(candidate.remainingStages()):
            for soln in semiGreedySolutions(candidate, depth):
                yield soln

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
        assert (not candidate) or (candidate.complete)
        if (not candidate) or (not (candidate < bestKnown)):
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

            # Do only one round of improvement; assumption is the analyst
            # already looped.
            if improvement and improvement.head.fullMass < newBest.currentMass:
                asPartial = partialSolution(None, LinkedList(improvement))
                assert (asPartial < newBest)
                newBest = asPartial
                analysis = analyst.analyze(newBest.stages)

            analyst.prettyPrint(newBest.stages, analysis)
            return (newBest, False, profiles)


    try:
        roots = makeRoots(profiles)
        nexpansions = 0
        nexpansionsLastPrinted = 0
        print ("Starting to search %d configurations" % (len(roots)))
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
def splitBurns(burns, minStageDeltaV):
    profiles = []
    # subburns are rawburn instances, but stored in a cons list, so they're in
    # top-first order rather than bottom-first.  This is conveniently how all
    # the internal code wants it.
    def recur(subburns, burnIdx):
        if burnIdx == len(burns):
            # Base case:
            # convert the subburns to a complete burn profile
            profile = burnProfile(subburns)
            profiles.append(profile)
        else:
            # Inductive case: for every split of burn i, recur.
            b = burns[burnIdx]
            maxstages = int(math.ceil(b.deltaV / minStageDeltaV)) # max number of splits
            for nstages in xrange(1, maxstages+1):
                mysubburns = subburns
                for subb in b.convert(nstages):
                    mysubburns = cons(subb, mysubburns)
                recur(mysubburns, burnIdx + 1)

    recur(nil, 0)
    return profiles

# Local search trick, and pretty printing.
class analyst(object):
    def __init__(self, burns, minStageDeltaV):
        """
        Set up, relative to the *original* burns that the user wants.
        Said burns are ordered bottom up.
        """
        self.burns = burns
        self.totalDeltaV = sum(b.deltaV for b in self.burns)
        self.minStageDeltaV = minStageDeltaV

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
                if s.acceleration() < b.acceleration:
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


    def _locallyImprove(self, stages, analysis):
        """
        Rebuild the stages, trying to reduce capacity to drop wasted
        deltaV.

        Returns a new set of stages if we succeed at reducing the mass,
        None if we fail.
        """
        # TODO: currently we only look at excess off the end of a series of
        # stages.  We could look at excess off the start.
        revStages = tuple(reversed(tuple(stages)))
        revAnalysis = tuple(reversed(analysis))
        (lastStage, lastData) = (revStages[0], revAnalysis[0])
        newStages = nil


        def rebuildStage(deltaV, payload, oldStage, laterEngines):
            """
            Copy the old stage, but change its target deltaV and payload.
            """
            if not oldStage.asparagus:
                laterEngines = None

            try:
                newS = stage(deltaV, payload,
                    oldStage.engineType, oldStage.numEngines, laterEngines,
                    oldStage.numTowers, oldStage.planet, oldStage.altitude)
            except engine.WeakEngineException:
                # This implies roundoff error, since we didn't change the engines.
                newS = None
            if newS.propellantMass > oldStage.propellantMass:
                # This implies roundoff, since we reduced the mass above.
                newS = None
            if not newS:
                # I was hoping we'd get the same deltaV with less
                # propellant, but roundoff bit me; just use the same
                # propellant as before.
                newS = stage(deltaV, payload,
                    oldStage.engineType, oldStage.numEngines, laterEngines,
                    oldStage.numTowers, oldStage.planet, oldStage.altitude,
                    propMassOverride = (s.propellantMass, s.tankMass))
            return newS

        # Rebuilding all the stages:
        # We start at the top and go down.
        for (revStageIdx, (s, d)) in enumerate(zip(revStages, revAnalysis)):
            laterEngines = stage.collectUsableEngines(newStages)
            payload = (newStages.head.fullMass if newStages else 0)

            if d.burnIds:
                # This stage is useful.  However, its mass may be reduced by
                # taking less propellant along.  Cruise stages just before
                # a descent or ascent stage typically have some slop.  Even if
                # there's no slop, stages above might now be lighter.
                targetDeltaV = s.achievedDeltaV() - d.dumpDV
                payload = d.payload + (newStages.head.fullMass if newStages else 0)
                laterEngines = stage.collectUsableEngines(newStages)

                newS = rebuildStage(targetDeltaV, payload, s, laterEngines)
                newStages = cons(newS, newStages)
            else:
                # Completely wasted stage.  Do we actually need the thrust for
                # lower stages?
                # For each lower stage, if it's asparagus with us, it might
                # need our thrust.
                # Check the acceleration without the mass or thrust of our
                # engines and fuel.  If that's not sufficient for the burns it
                # participates in, then we need our thrust (but not the fuel).
                # Things are slightly better in reality: we reduce our mass, so
                # lower stages can reduce theirs too!  So we could do better
                # than this estimate.
                thrust = s.engineType.thrust * s.numEngines
                vthrust = s.engineType.thrust * s.numEngines if s.engineType.vectoring else 0
                stageMassDiff = s.fullMass - s.payload
                def canIgnoreStage():
                    for (lowerStage, lowerData) in zip(revStages[revStageIdx+1:], revAnalysis[revStageIdx:]):
                        if not lowerStage.asparagus: break

                        massNoStage = lowerStage.fullMass - stageMassDiff
                        thrustNoStage = lowerStage.thrust - thrust
                        vthrustNoStage = lowerStage.vectoringThrust - vthrust
                        if vthrustNoStage < 0.25 * thrustNoStage:
                            return False # not enough vectoring thrust without the stage
                        for bId in lowerData.burnIds:
                            accel = self.burns[bId].acceleration
                            if thrustNoStage < massNoStage * accel:
                                return False # not enough thrust:mass ratio without the stage
                    return True # no reason it can't work => it works!
                if canIgnoreStage():
                    # Yay, we ditched the fuel *and* the engines *and* the decouplers!
                    continue
                else:
                    # Make a 0-deltaV stage with just the engines instead.
                    # It's too much trouble for too little improvement to merge
                    # engines with the next stage below.
                    zeroStage = rebuildStage(0, payload, s, laterEngines)
                    newStages = cons(zeroStage, newStages)


        if newStages.head.fullMass < stages.head.fullMass:
            print ("Locally improved from %g T to %g T" %
                    (stages.head.fullMass, newStages.head.fullMass))
            return newStages

    def _suggestDeltaVs(self, stages):
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

    def suggest(self, stages, analysis):
        """
        * First, try reducing the propellant in each stage to eliminate
          wasted deltaV (which is typically at the end, or just before a
          descent/ascent stage).

        TODO: reinstate these!
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
        # Local improvement.  Keep improving as long as we can!
        improved = self._locallyImprove(stages, analysis)
        while improved:
            stages = improved
            analysis = self.analyze(stages)
            improved = self._locallyImprove(stages, analysis)


#        def makeProfile(deltaVs):
#            """
#            Return a burn profile that corresponds to achieving the given
#            deltaVs, and matching up with acceleration requirements.
#            """
#            def makeBurn(deltaV, (fullBurns, partialBurn)):
#                if partialBurn:
#                    allBurns = [ b for b in fullBurns ]
#                    allBurns.append(partialBurn)
#                else:
#                    allBurns = fullBurns
#                assert len(allBurns) > 0
#
#                # TODO: this is pessimistic.  Our acceleration will increase
#                # over the stage.  However, we have no way to convey such
#                # fine-grained information in the profile being optimized.
#                accel = max( b.accel for b in allBurns )
#
#                if len(allBurns) == 1:
#                    name = allBurns[0].name
#                else:
#                    name = ",".join(b.name for b in allBurns)
#
#                # Payloads that we have to lift, not counting later stages.
#                # We don't care about the partial burn, which a later stage
#                # will handle.
#                payload = sum( b.payload for b in fullBurns ) if fullBurns else 0
#
#                return burn(deltaV, accel, name, payload)
#
#            # The deltaVs and burnsByStage are bottom-up, we need the burns top-down.
#            burnsByStage = self._deltaVMap(deltaVs)
#            newburns = [  makeBurn(*x) for x in reversed(zip(deltaVs, burnsByStage)) ]
#            assert len(newburns) > 0
#            return burnProfile(newburns, startAltitude = startAltitude)
#
#        # reverse the stages (painful)
#        revStages = tuple(reversed(tuple(stages)))
#
#        # suggest deltaV lists in both directions
#        deltaVsUp = self._suggestDeltaVs(stages)
#        deltaVsDown = tuple(reversed(self._suggestDeltaVs(revStages)))
#
        # Return the improved set of stages (or the original set), and also
        # return the potentially better profiles to search.
        return ( stages, [] ) #( makeProfile(deltaVsUp), makeProfile(deltaVsDown) ) )


    def prettyPrint(self, stages, data):
        m0 = self.burns[-1].payload
        m1 = stages.head.fullMass
        Isp = self.totalDeltaV / (physics.g0 * math.log(m1 / m0)) # log is ln
        print ("Solution in %d stages:" % (len(stages),))
        print ("  Final payload: %g T" % m0)
        print ("  Launch mass:   %g T" % m1)
        print ("  Ratio:         %g" % (m1/m0))
        print ("  Implied Isp:   %g" % Isp)
        # TODO: we're not counting waste when we're forced to dump fuel.

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
            if d.dumpDV:
                if sIdx != len(data) - 1:
                    print ("\t* must dump %g m/s before next burn" % d.dumpDV)
                else:
                    print ("\t* end mission with %g m/s to spare" % d.dumpDV)
            if d.payload:
                print ("\t* includes %g T payload" % d.payload)


def design(burns, minStageDeltaV = 750, symmetry = 2, numBaseTowers = 1):
    """
    Design a rocket to perform the given burns, instances of liftoffBurn and
    deepSpaceBurn.  Prints to stdout.

    symmetries details the symmetry required.  By default, we try 2-way symmetry.

    numBaseTowers says how many towers we can start with.  By default, we try 1 tower.
        Increase this if you want to stick the rocket around the sides of your payload,
        rather than under it.

    minStageDeltaV gives a limit to how much splitting to perform: we will split
        any stage with more than 750 m/s dV into a variety of numbers of
        smaller stages.
        It is probably a bad algorithm...
    """
    print ("Designing for burns totalling %g m/s, payload total %g T" %
        (sum(b.deltaV for b in burns),
         sum(b.payload for b in burns)))
    profiles = splitBurns(burns, minStageDeltaV)
    print ("%g Profiles:" % (len(profiles)))
    for p in profiles:
        print str(p)
    shrink = analyst(burns, minStageDeltaV)

    # This code takes a long, long time...  Hit C-c at the shell to stop it and
    # return the best yet.
    soln = designRocket(profiles,
                massToBeat = None,
                analyst = shrink,
                symmetries = symmetry,
                numBaseTowers = numBaseTowers)

    if soln:
        data = shrink.analyze(soln)
        shrink.prettyPrint(soln, data)
    else:
        print "You will not be going to space today."
