from __future__ import division

import math

import engine
import jets
import lift
import planet
import physics

class part(object):
    def __init__(self, parttype, n = 1, AoAdegrees = 0, extraMass = 0, Cd = None):
        self._parttype = parttype
        self._n = n
        self._AoA = AoAdegrees
        self._extraMass = extraMass
        self._staticCd = Cd

    def mass(self):
        if self._parttype is None:
            return self._n * self._extraMass
        else:
            return self._n * (self._parttype.mass + self._extraMass)

    def zeroThrustForces(self, flightPitchDegrees, AoADegrees, v, altitude, planet):
        """
        Return the forces other than thrust that this part produces.
        Returns a tuple (drag, lift, apparentGravity) with quantities in kN.
        Lift and drag are vectors; gravity is a number (just the y coordinate, positive means up).
        """
        Cd = self._staticCd
        Clift = 0
        if isinstance(self._parttype, lift.wing):
            Cd = self._parttype.dragCoeff(self._AoA + AoADegrees)
            Clift = self._parttype.liftCoeff(self._AoA + AoADegrees)
        elif isinstance(self._parttype, jets.intake):
            Cd = self._parttype.dragCoeff(self._AoA + AoADegrees)
        elif isinstance(self._parttype, jets.jetengine):
            Cd = 0.2
        elif isinstance(self._parttype, engine.engine):
            Cd = 0.2

        dragMagnitude = planet.dragForce(altitude, v, self.mass(), Cd)
        liftMagnitude = v * planet.pressure(altitude) * Clift

        flightPitchRad = math.radians(flightPitchDegrees)
        cosPitch = math.cos(flightPitchRad)
        sinPitch = math.sin(flightPitchRad)

        horizontalSpeedSrf = v * cosPitch
        horizontalSpeedOrb = horizontalSpeedSrf + planet.siderealSpeed(altitude)
        vOrbit = planet.orbitalVelocity(altitude)
        centrifuge = horizontalSpeedOrb * horizontalSpeedOrb / (vOrbit * vOrbit)
        gravity = planet.gravity(altitude)
        downForceMagnitude = self.mass() * (centrifuge - gravity)

        dragVector = physics.vector(-cosPitch * dragMagnitude, -sinPitch * dragMagnitude)
        liftVector = physics.vector(-sinPitch * liftMagnitude, cosPitch * liftMagnitude)

        return (dragVector, liftVector, downForceMagnitude)

    def thrustForceAtMaxThrottle(self, v, flightPitchDegrees, AoAdegrees):
        if isinstance(self._parttype, jets.jetengine):
            t = self._parttype.thrust(v)
        elif isinstance(self._parttype, engine.engine):
            t = self._parttype.thrust
        else:
            t = 0

        enginePitchRad = math.radians(flightPitchDegrees + AoAdegrees + self._AoA)
        return physics.vector(math.cos(enginePitchRad) * t, math.sin(enginePitchRad) * t)

    def airProvided(self, AoADegrees, v, altitude, options):
        if isinstance(self._parttype, jets.intake):
            return self._n * self._parttype.intakeair(altitude, v, AoADegrees + self._AoA, options)
        else:
            return 0

    def airRequired(self, v, altitude, options):
        if isinstance(self._parttype, jets.jetengine):
            return self._n * self._parttype.airRequired(altitude, options = options)
        else:
            return 0

    def fuelRequired(self, altitude, options):
        if isinstance(self._parttype, jets.jetengine):
            return self._n * self._parttype.fuelRequired(altitude, options = options)
        else: return 0

class plane(object):
    def __init__(self, *parts):
        self._parts = parts

    def forces(self, flightPitch, AoA, v, altitude = 0, jetoptions = jets.kerbonormative, throttle = 1, verbose = False):
        """
        Given:
        the flight path pitch (the angle of prograde with the surface),
        the angle of attack (the angle of pitch with prograde),
        the speed,
        the altitude and the planet,
        the throttle (default max we can achieve at this altitude)

        Prints a report if 'verbose' is set.

        Returns a tuple:
            (net, throttle, lift, drag, thrust, apparentGravity)
        net, lift, drag, and thrust are vectors with X being forward and Y being up.
        throttle is, from 0 to 1, how much the jet engines will be throttled up assuming proper intake layout
        apparentGravity is the force of gravity minus centrifugal effects.  Normally negative.
        """
        drag = physics.vector(0,0)
        lift = physics.vector(0,0)
        down = 0
        maxthrust = physics.vector(0,0)
        airProvided = 0
        airRequired = 0
        for part in self._parts:
            (d,l,g) = part.zeroThrustForces(flightPitch, AoA, v, altitude, jetoptions.planet)
            drag = drag.add(d)
            lift = lift.add(l)
            down += g
            maxthrust = maxthrust.add(part.thrustForceAtMaxThrottle(v, flightPitch, AoA))
            airProvided += part.airProvided(AoA, v, altitude, jetoptions)
            airRequired += part.airRequired(v, altitude, jetoptions)

        # account for throttle setting due to air
        if throttle * airRequired > airProvided:
            throttle = airProvided / airRequired
        thrust = maxthrust.scale(throttle)

        net = thrust.add(lift.add(drag.add(physics.vector(0, down))))

        if verbose:
            print ("throttle: %.1f%%" % (throttle * 100))
            print ("lift:    (%7.2f, %7.2f) kN" % (lift[0], lift[1]))
            print ("drag:    (%7.2f, %7.2f) kN" % (drag[0], drag[1]))
            print ("thrust:  (%7.2f, %7.2f) kN" % (thrust[0], thrust[1]))
            print ("gravity: (       , %7.2f) kN" % down)
            print ("net:     (%7.2f, %7.2f) kN" % (net[0], net[1]))

        return (net, throttle, lift, drag, thrust, down)

    def equilibriumHorizontalSpeed(self, AoA, altitude = 0, jetoptions = jets.kerbonormative, throttle = 1):
        """
        Given:
        angle of attack
        optional altitude
        optional throttle
        optional planet and physics deltaT

        Find the speed at which drag and thrust balance perfectly.
        """
        def willSpeedUp(v):
            # we expect zero flight path angle
            netForce = self.forces(0, AoA, v, altitude, jetoptions, throttle)[0]
            if netForce[0] > 0: return True
            else: return False

        # start at 100 m/s, keep doubling until we find a speed where we don't speed up.
        vmin = 0
        vmax = 100
        while willSpeedUp(vmax):
            vmin = vmax
            vmax = vmax * 2

        # Now we know the equilibrium lies between vmin and vmax; binary search.
        while vmax - vmin > 1e-3:
            vhalf = 0.5 * (vmax + vmin)
            if willSpeedUp(vhalf):
                vmin = vhalf
            else:
                vmax = vhalf

        return vhalf

    def levelThrottle(self, AoA, altitude = 0, jetoptions = jets.kerbonormative):
        """
        Given:
        angle of attack
        optional altitude
        optional planet and physics deltaT

        Find the throttle setting in [0,1] that has us be at equilibrium:
        thrust and lift offset apparent gravity when at the equilibrium horizontal speed.

        Return None if that's not possible, i.e. you need more than throttle of 1.

        Assumption: engines aren't pointed down.
        """
        def isFalling(t):
            v = self.equilibriumHorizontalSpeed(AoA, altitude, jetoptions, t)
            netForce = self.forces(0, AoA, v, altitude, jetoptions, t)[0]
            if netForce[1] < 0:
                return True
            else:
                return False

        if isFalling(1):
            return None
        if not isFalling(0):
            return 0

        minThrottle = 0
        maxThrottle = 1
        while maxThrottle - minThrottle > 1e-4:
            half = 0.5 * (maxThrottle + minThrottle)
            if isFalling(half):
                minThrottle = half
            else:
                maxThrottle = half
        return half

    def fuelConsumption(self, altitude = 0, throttle = 1, jetoptions = jets.kerbonormative):
        """
        How much fuel does this plane consume at the given throttle setting and altitude?
        Returns a value in kg/s.
        """
        fuel = 0
        for part in self._parts:
            fuel += part.fuelRequired(altitude, options = jetoptions)
        return throttle * fuel

    def optimizeFuelConsumptionAtAltitude(self, altitude = 0, jetoptions = jets.kerbonormative, verbose = False):
        """
        Given that we're flying at a certain altitude, achieve
        level flight for the minimum fuel consumed per distance.
        Returns a tuple:
            (pitch [0,90] degrees, throttle [0,1], speed, fuel consumption in kg/m)
        Returns None if you can't get level flight at that altitude no matter what.
        """
        bestPitch = -1
        bestThrottle = 1
        bestSpeed = 0
        bestFuel = 1e308

        for pitch in range(90):
            if verbose: print "Trying pitch %d, best so far is %d" % (pitch, bestPitch)
            throttle = self.levelThrottle(pitch, altitude = altitude, jetoptions = jetoptions)
            if throttle is None:
                continue
            speed = self.equilibriumHorizontalSpeed(pitch, altitude = altitude, throttle = throttle, jetoptions = jetoptions)
            if speed < 1e-6:
                # no stationary hovering allowed
                continue
            fuelPerTime = self.fuelConsumption(altitude, throttle, jetoptions)
            fuelPerDist = fuelPerTime / speed
            if fuelPerDist < bestFuel:
                bestPitch = pitch
                bestThrottle = throttle
                bestSpeed = speed
                bestFuel = fuelPerDist

        if bestPitch == -1:
            return None
        else:
            return (bestPitch, bestThrottle, bestSpeed, bestFuel)

    def optimizeFuelConsumption(self, minAltitude = 0, maxAltitude = None, altitudeStep = 500, jetoptions = jets.kerbonormative, verbose = False):
        """
        Optimize the fuel consumption over all altitudes.

        Optionally specify the min altitude (so you don't fly through terrain),
        the max altitude, and the step between altitudes that we consider.
        """
        if maxAltitude is None:
            maxAltitude = jetoptions.planet.topOfAtmosphere()

        bestAltitude = -1
        bestResult = None

        altitude = minAltitude
        while altitude < maxAltitude:
            if verbose: print ("Testing altitude %d, best so far is %d with %s" % (altitude, bestAltitude, bestResult))
            result = self.optimizeFuelConsumptionAtAltitude(altitude = altitude, jetoptions = jetoptions)
            if result is None:
                break
            if bestResult is None or result[3] < bestResult[3]:
                bestAltitude = altitude
                bestResult = result
            altitude += altitudeStep

        if bestResult is None:
            if verbose: print ("Level flight is impossible with this plane")
            return None

        (pitch, throttle, speed, fuelKgPerM) = bestResult

        if verbose:
            fuelUnitsPerKm = fuelKgPerM * 1000 / 5
            print ("Fly at %d m altitude, pitch %d degrees, throttle %.3f, to achieve %.2f m/s using %g U/km" %
                    (bestAltitude, pitch, throttle, speed, fuelUnitsPerKm))

        return (bestAltitude, pitch, throttle, speed, fuelKgPerM)


###########################################################################
### Below here it's all just tests to run each routine.

def forcesTest():
    """
    The plane has a basic jet and two radial intakes.
    It's got 2 AV-R8 for its elevator, two deltaDeluxe for ailerons,
    and two deltaDeluxe for wings, which are tilted 10 degrees.
    There's also two fuel tanks, full.
    We're climbing at a 5 degree angle, with pitch 7.5 degrees (so angle of attack 2.5 degrees).
    We're going Mach 0.85, at 10km altitude.

    This plane is not at equilibrium: we've got horizontal speed to gain,
    and vertical speed is going down.  Specifically, we should see:
    net:     (  25.72,  -15.09) kN
    """
    import planeDesigner, jets, lift
    p = plane(planeDesigner.part(jets.basicjet),
             planeDesigner.part(jets.radialIntake, 2),
             planeDesigner.part(lift.avr8, 2),
             planeDesigner.part(lift.deltaDeluxe, 2),
             planeDesigner.part(lift.deltaDeluxe, 2, AoAdegrees = 10),
             planeDesigner.part(None, n = 2, Cd = 0.2, extraMass = 0.15 + 150 * 0.005), # two fuel tanks
            )
    p.forces(5, 2.5, 0.85 * 342, altitude = 10000, throttle = 0.5)

def speedTest():
    import planeDesigner, jets, lift
    p = plane(planeDesigner.part(jets.basicjet),
             planeDesigner.part(jets.radialIntake, 2),
             planeDesigner.part(lift.avr8, 2),
             planeDesigner.part(lift.deltaDeluxe, 2),
             planeDesigner.part(lift.deltaDeluxe, 2, AoAdegrees = 10),
             planeDesigner.part(None, n = 2, Cd = 0.2, extraMass = 0.15 + 150 * 0.005), # two fuel tanks
            )
    speed = p.equilibriumHorizontalSpeed(2.5, altitude = 10000, throttle = 0.5)
    print ("equilibrium speed: %g m/s" % speed)
    p.forces(0, 2.5, speed, altitude = 10000, throttle = 0.5)

def throttleTest():
    import planeDesigner, jets, lift
    p = plane(planeDesigner.part(jets.basicjet),
             planeDesigner.part(jets.radialIntake, 2),
             planeDesigner.part(lift.avr8, 2),
             planeDesigner.part(lift.deltaDeluxe, 2),
             planeDesigner.part(lift.deltaDeluxe, 2, AoAdegrees = 10),
             planeDesigner.part(None, n = 2, Cd = 0.2, extraMass = 0.15 + 150 * 0.005), # two fuel tanks
            )
    throttle = p.levelThrottle(10, altitude = 10000)
    if throttle is None:
        print ("level flight is impossible")
    else:
        speed = p.equilibriumHorizontalSpeed(2.5, altitude = 10000, throttle = throttle)
        print ("equilibrium speed: %g m/s" % speed)
        p.forces(0, 2.5, speed, altitude = 10000, throttle = throttle)


def optimizePitchTest():
    import planeDesigner, jets, lift
    p = plane(planeDesigner.part(jets.basicjet),
             planeDesigner.part(jets.radialIntake, 2),
             planeDesigner.part(lift.avr8, 2),
             planeDesigner.part(lift.deltaDeluxe, 2),
             planeDesigner.part(lift.deltaDeluxe, 2, AoAdegrees = 10),
             planeDesigner.part(None, n = 2, Cd = 0.2, extraMass = 0.15 + 150 * 0.005), # two fuel tanks
            )
    opt = p.optimizeFuelConsumptionAtAltitude(altitude = 10000, verbose = True)
    if opt is None:
        print ("Level flight is impossible with any pitch and throttle at this altitude")
    else:
        (pitch, throttle, speed, fuelKgPerM) = opt
        fuelUnitsPerKm = fuelKgPerM * 1000 / 5
        print ("Fly at %d degrees, throttle %.3f, to achieve %.2f m/s using %g U/km" % (pitch, throttle, speed, fuelUnitsPerKm))

def optimizeFullyTest():
    import planeDesigner, jets, lift
    # use the same plane as above but with a turbojet because that shows more interesting behaviour
    p = plane(planeDesigner.part(jets.turbojet),
             planeDesigner.part(jets.radialIntake, 2),
             planeDesigner.part(lift.avr8, 2),
             planeDesigner.part(lift.deltaDeluxe, 2),
             planeDesigner.part(lift.deltaDeluxe, 2, AoAdegrees = 10),
             planeDesigner.part(None, n = 2, Cd = 0.2, extraMass = 0.15 + 150 * 0.005), # two fuel tanks
            )
    p.optimizeFuelConsumption(verbose = True, altitudeStep = 1000)
