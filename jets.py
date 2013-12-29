# KSP Jet engine information.
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

#
# This file is set up in two parts: first, information about the jet engines
# and intakes; second, functions that help you design a plane.
#

from __future__ import division # / means float div always
import math
from physics import g0, PiecewiseLinearCurve
import planet

class standardoptions(object):
    """
    Represent some standard options that don't depend on which
    function you're calling, but which affect most of them.
    """
    def __init__(self, planet, deltaT):
        self.planet = planet
        self.deltaT = deltaT

kerbonormative = standardoptions(planet.kerbin, 0.03)
laythenormative = standardoptions(planet.laythe, 0.03)

class jetengine(object):
    def __init__(self, 
            name,
            ispCurve,   # key-value pairs
            thrustCurve,# key-value pairs
            maxthrust,
            mass): # kN
        self.name        = name
        self.ispCurve    = ispCurve
        self.thrustCurve = thrustCurve
        self.maxthrust   = maxthrust
        self.mass        = mass

    def Isp(self, altitude, options = kerbonormative):
        pressure = options.planet.pressure(altitude)
        return self.ispCurve.evaluate(pressure)

    def thrust(self, v):
        """
        Return the thrust (kN) provided by the jet at the given speed.
        """
        return self.maxthrust * self.thrustCurve.evaluate(v)

    def airRequired(self, altitude, throttle = 1, options = kerbonormative):
        """
        Intake air (kg/s) required at the given 
        altitude (m) and throttle ([0-1], default 1).
        """
        # Total mass used is:
        #       mdot = throttle * maxthrust / (Isp g0)
        # Air is 15/16 = 0.9375 of that.
        # maxthrust is kN so multiply by 1000 to get kg/s
        # fold all the constants together: 1000 * 0.9375 / 9.82 = 95.4...
        Isp = self.Isp(altitude, options = options)
        return 95.46843177189409 * throttle * self.maxthrust / Isp

    def fuelRequired(self, altitude, throttle = 1, options = kerbonormative):
        """
        LiquidFuel (kg/s) required at the given
        altitude (m) and throttle ([0-1], default 1).
        """
        # Fuel:air is 1:15, so we need 1/15 of the amount of air we need.
        return self.airRequired(altitude, throttle, options) / 15.0

turbojet = jetengine("TurboJet Engine",
    PiecewiseLinearCurve.fromfile("data/turbojet-isp.txt"),
    PiecewiseLinearCurve.fromfile("data/turbojet-thrust.txt"),
    225, # max thrust
    1.2 # mass
)

basicjet = jetengine("Basic Jet Engine",
    PiecewiseLinearCurve.fromfile("data/jet-isp.txt"),
    PiecewiseLinearCurve.fromfile("data/jet-thrust.txt"),
    150, # max thrust
    1 # mass
)

rapier = jetengine("RAPIER atmospheric phase",
    PiecewiseLinearCurve.fromfile("data/turbojet-isp.txt"),
    PiecewiseLinearCurve.fromfile("data/rapierjet-thrust.txt"),
    190,
    1.75
)

class intake(object):
    def __init__(self, name, drymass, area, capacity):
        self.name = name
        self.area = area         # m^2
        self.capacity = capacity # units
        self.drymass = drymass      # tonnes

    def intakeair(self, altitude, v, AoA = 0, options = kerbonormative):
        """
        Return the intake air (kg/s) provided at the given
        altitude (m) and airspeed (m/s), assuming the given angle 
        of attack (radians, default 0) and physics deltaT (s, default 0.03)

        The angle of attack must be less than pi/2, otherwise a different
        formula would hold.
        """
        # For low angle of attack, the formula is:
        #       I = area * airdensity * cos(AoA) * u^2 (6 intakespeed + v)
        # in t/s ; multiply by 1000 to get kg/s
        # u is 0.2 for all intakes; named unitscalar in the code.  No idea what
        # it's about.

        D = options.planet.airdensity(altitude)

        # After folding the u = 0.2, intakespeed = 100, and the *1000, 
        # this is what we get:
        air = self.area * D * math.cos(AoA) * 40 * (600 + v)

        # The intake can hold no more than its capacity at any time step.
        # This means your plane flies better with a smaller deltaT!
        maxair = self.capacity * 5 / options.deltaT

        return min(air, maxair)

    def drag(self, altitude, v, AoA = 0, options = kerbonormative):
        """
        Return the drag (kN) caused by this intake at the given altitude and
        velocity.
        """
        # The total mass is the dry mass of the intake, plus the mass of the
        # air in the intake!  This actually matters at low altitude.

        # intakeair is kg/s, but the amount we currently have in the intake
        # is the amount we gather in one timestep.  And we need it in tonnes,
        # not kg.
        airkgpers = self.intakeair(altitude, v, AoA, options=options)
        airtonnes = 0.001 * options.deltaT * airkgpers
        mass = self.drymass + airtonnes
        # TODO: the coefficient is not always 2, but I don't know the formula.
        dragCoeff = 2
        return options.planet.dragForce(altitude, v, mass, 0.3 + dragCoeff)

ramAirIntake = intake("Ram Air Intake", 0.01, 0.01, 0.2)
radialIntake = intake("XM-G50 Radial Air Intake", 0.01, 0.004, 1.0)
circularIntake = intake("Circular Intake", 0.01, 0.008, 0.2)
nacelleIntake = intake("Engine Nacelle", 0.3, 0.002, 0.2)


###########################################################################
# Jet design functions.

class planeInfo(object):
    """
    Compute a bunch of information about a plane from its list of parts
    and its mass (mass should include engines and intakes).
    """
    def __init__(self, altitude, v, mass, parts, AoA = 0, options = kerbonormative):
        provided = 0
        required = 0
        maxthrust = 0
        intakedrag = 0
        intakemass = 0
        for part in parts:
            try:
                # assume it's a list-like thing
                typ = part[0]
                n = part[1]
            except TypeError:
                # not a list-like thing, must be just a single part
                typ = part
                n = 1
            if isinstance(typ, jetengine):
                maxthrust += n * typ.thrust(v)
                required += n * typ.airRequired(altitude, options = options)
            else:
                provided += n * typ.intakeair(altitude, v, AoA = AoA, 
                                        options = options)
                intakedrag += n * typ.drag(altitude, v, AoA = AoA, 
                                        options = options)
                intakemass += n * typ.drymass

        self.maxthrust = maxthrust      # kN
        self.required = required        # kg/s of air
        self.provided = provided        # kg/s of air

        if required <= provided:
            self.maxthrottle = 1        # [0,1]
            self.thrust = maxthrust     # kN
        else:
            maxthrottle = provided / required
            self.maxthrottle = maxthrottle
            self.thrust = maxthrust * maxthrottle
        
        self.intakedrag = intakedrag    # kN
        self.intakemass = intakemass    # tonnes, dry mass of intakes
        self.stdmass = mass - intakemass# tonnes, not including intakes
        self.totalMass = mass           # tonnes, not including mass of the air

        # kN, total drag on the airframe
        self.dragForce = (intakedrag + 
            options.planet.dragForce(altitude, v, self.stdmass, 0.2))

        # kN, net force on the airframe, >0 means accelerate, <0 means
        # decelerate.
        self.netForce = self.thrust - self.dragForce

        # m/s^2 at the top achievable throttle
        if abs(self.netForce) < 1e-6:
            self.acceleration = 0
        else:
            self.acceleration = self.totalMass / self.netForce



def maxSpeed(parts, mass, altitude, tolerance = 1e-3, options = kerbonormative):
    """
    Return a speed close to the max speed.
    The speed is guaranteed to have positive net force, but to be within
    tolerance (in m/s) of a speed with negative net force.

    jets.maxSpeed( (jets.turbojet, (jets.ramAirIntake, 2)), 5, 25250 )
    1858.9853515625
    """
    maxVWithPosForce = 0
    minVWithNegForce = None
    v = 1

    # binary search to define the upper bound of the binary search
    # could do better via physics, but this will work
    while minVWithNegForce is None:
        v *= 2
        plane = planeInfo(altitude, v, mass, parts)
        if plane.netForce >= 0:
            maxVWithPosForce = v
        else:
            minVWithNegForce = v

    # binary search to find the exact max speed, to within tolerance
    while minVWithNegForce - maxVWithPosForce > tolerance:
        v = 0.5*(minVWithNegForce + maxVWithPosForce)
        plane = planeInfo(altitude, v, mass, parts)
        if plane.netForce >= 0:
            maxVWithPosForce = v
        else:
            minVWithNegForce = v

    return maxVWithPosForce

def maxSpeedAltitude(parts, mass, tolerance = 1e-3, options = kerbonormative):
    """
    Return the max speed and the altitude that can provide it.

    jets.maxSpeedAltitude( (jets.turbojet, (jets.ramAirIntake, 2)), 5)
    (1860.71484375, 24631.498671722413)
    """
    # top speed using turbojets is at low as you can while running out of air
    def speed(altitude):
        return maxSpeed(parts, mass, altitude, tolerance, options)
    def data(altitude):
        v = speed(altitude)
        plane = planeInfo(altitude, v, mass, parts, options=options)
        plane.v = v
        plane.altitude = altitude
        return plane

    left = 0
    right = options.planet.topOfAtmosphere() - 1

    while right - left > 2:
        mid = 0.5 * (right + left)
        midPlane = data(mid)
        if midPlane.required > midPlane.provided:
            right = mid
        else:
            left = mid

    return midPlane.v, mid


def intakesForSpeed(parts, mass, targetSpeed = None, tolerance = 1e-3, options = kerbonormative):
    """
    Given an initial configuration, return the min number of intakes to add in
    order to achieve the target speed.
    If target speed is not stated, it's the speed of a circular orbit at the
    top of the atmosphere.
    """
    if targetSpeed is None:
        targetSpeed = options.planet.orbitalVelocity(options.planet.topOfAtmosphere()) - options.planet.siderealRotationSpeed

    (speed, _) = maxSpeedAltitude(parts, mass, tolerance, options)
    if speed > targetSpeed:
        return 0

    # todo: binary search, but I am not sure the intake benefit grows monotonically
    for i in range(1000):
        iparts = [(ramAirIntake, i)] + [x for x in parts]
        (speed, _) = maxSpeedAltitude(iparts, mass, tolerance, options)
        if speed > targetSpeed:
            return i
    raise Exception("impossible to reach " + targetSpeed + " even with 1000 intakes")


def machingbird(parts, mass, tolerance = 1e-3, options = kerbonormative):
    speedAlt = maxSpeedAltitude(parts, mass, tolerance, options)

    orbitSpeed = options.planet.orbitalVelocity(speedAlt[1])
    siderealSpeed = options.planet.siderealRotationSpeed
    
    # speed to gain to reach orbit at the optimal altitude
    subOrbitSpeed = orbitSpeed - speedAlt[0]
    if subOrbitSpeed > siderealSpeed:
        inclination = 0
    elif subOrbitSpeed < -siderealSpeed:
        inclination = 180
    else:
        inclination = math.degrees(math.acos(subOrbitSpeed / siderealSpeed))

    return speedAlt[0], speedAlt[1], inclination
