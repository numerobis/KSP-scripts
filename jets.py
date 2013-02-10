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
            maxthrust): # kN
        self.name        = name
        self.ispCurve    = PiecewiseLinearCurve(ispCurve)
        self.thrustCurve = PiecewiseLinearCurve(thrustCurve)
        self.maxthrust   = maxthrust

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
        # fold all the constants together: 1000 * 0.9375 / 9.81 = 95.5...
        Isp = self.Isp(altitude, options = options)
        return 95.565749235474 * throttle * self.maxthrust / Isp


turbojet = jetengine("TurboJet Engine",
    ( (1, 800), (0.3, 2500), (0, 1200) ), # isp curve
    ( (0, 0.5), (1000, 1), (2000, 0.5), (2400, 0) ), # thrust curve
    225 # max thrust
)

basicjet = jetengine("Basic Jet Engine",
    ( (1, 2000), (0.3, 1800), (0, 1000) ), # isp
    ( (0, 1), (850, 0.2), (1000, 0) ), # thrust
    150 # max thrust
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
        return options.planet.dragForce(altitude, v, mass, 2.0)

ramAirIntake = intake("Ram Air Intake", 0.01, 0.01, 0.2)
radialIntake = intake("XM-G50 Radial Air Intake", 0.1, 0.004, 1.0)
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

        if required >= provided:
            self.maxthrottle = 1        # [0,1]
            self.thrust = maxthrust     # kN
        else:
            maxthrottle = required / provided
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
