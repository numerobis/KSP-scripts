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

def forces(parts, flightPitch, AoA, v, altitude = 0, planet = planet.kerbin, deltaT = 0.04, throttle = 1, quiet = False):
    """
    Given a description of the plane, 
    the flight path pitch (the angle of prograde with the surface),
    the angle of attack (the angle of pitch with prograde),
    the speed,
    the altitude and the planet, 
    the throttle (default max we can achieve at this altitude)

    Prints a report unless 'quiet' is set.

    Returns a tuple:
        (net, throttle, lift, drag, thrust, apparentGravity)
    'net', 'lift', 'drag', and 'thrust' are vectors with X being forward and Y being up.
    'throttle' is, from 0 to 1, how much the jet engines will be throttled up assuming proper intake layout
    'apparentGravity' is the force of gravity minus centrifugal effects.  Normally negative.
    """
    pitch = flightPitch + AoA
    jetOptions = jets.standardoptions(planet, deltaT)

    drag = physics.vector(0,0)
    lift = physics.vector(0,0)
    down = 0
    maxthrust = physics.vector(0,0)
    airProvided = 0
    airRequired = 0
    for part in parts:
        (d,l,g) = part.zeroThrustForces(flightPitch, AoA, v, altitude, planet)
        drag = drag.add(d)
        lift = lift.add(l)
        down += g
        maxthrust = maxthrust.add(part.thrustForceAtMaxThrottle(v, flightPitch, AoA))
        airProvided += part.airProvided(AoA, v, altitude, jetOptions)
        airRequired += part.airRequired(v, altitude, jetOptions)

    # account for throttle setting due to air
    if throttle * airRequired > airProvided:
        throttle = airProvided / airRequired
    thrust = maxthrust.scale(throttle)

    net = thrust.add(lift.add(drag.add(physics.vector(0, down))))

    if not quiet:
        print ("throttle: %.1f%%" % (throttle * 100))
        print ("lift:    (%7.2f, %7.2f) kN" % (lift[0], lift[1]))
        print ("drag:    (%7.2f, %7.2f) kN" % (drag[0], drag[1]))
        print ("thrust:  (%7.2f, %7.2f) kN" % (thrust[0], thrust[1]))
        print ("gravity: (       , %7.2f) kN" % down)
        print ("net:     (%7.2f, %7.2f) kN" % (net[0], net[1]))

    return (net, throttle, lift, drag, thrust, down)

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
    plane = (planeDesigner.part(jets.basicjet),
             planeDesigner.part(jets.radialIntake, 1),
             planeDesigner.part(lift.avr8, 2),
             planeDesigner.part(lift.deltaDeluxe, 2),
             planeDesigner.part(lift.deltaDeluxe, 2, AoAdegrees = 10),
             planeDesigner.part(None, n = 2, Cd = 0.2, extraMass = 0.15 + 150 * 0.005), # two fuel tanks
            )
    planeDesigner.forces(plane, 5, 2.5, 0.85 * 342, altitude = 10000, throttle = 0.5)
