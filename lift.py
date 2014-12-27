# KSP aerodynamic lift information.
# Copyright 2013 Benoit Hudson
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
import planet

# Given the angle of attack, return the angle of the velocity vector with the
# direction a normal person would call up (but Unity calls it forward).
def upDegrees(AoAdegrees):
    if AoAdegrees > 90:
        return 270 - AoAdegrees
    else:
        return 90 + AoAdegrees

# Wings are things that don't have control authority on their own.
# They are part = winglet in their part.cfg (but only one of the winglets
# matches it, so I renamed to wing for clarity).
class wing(object):
    def __init__(self, name, mass, lift, drag):
        self.name = name
        self.mass = mass
        self.lift = lift
        self.drag = drag

    @classmethod
    def deflectionLift(cls, AoAdegrees):
        cosAoA = math.cos(math.radians(AoAdegrees))
        cosUp = math.cos(math.radians(upDegrees(AoAdegrees)))
        return -cosUp * (1.0 - abs(cosUp)) * cosAoA

    def dragCoeff(self, AoAdegrees):
        return math.sin(math.radians(AoAdegrees)) * self.drag

    def liftCoeff(self, AoAdegrees):
        return self.deflectionLift(AoAdegrees) * self.lift

    def liftForce(self, AoAdegrees, v, altitude = 0, planet = planet.kerbin):
        return v * planet.pressure(altitude) * self.liftCoeff(AoAdegrees)

    def dragForce(self, AoAdegrees, v, altitude = 0, planet = planet.kerbin):
        Cd = self.dragCoeff(AoAdegrees)
        return planet.dragForce(altitude, v, self.mass, Cd)

    def forceVector(self, flightPitch, AoA, v, altitude = 0, planet = planet.kerbin, includeGravity = False):
        """
        Returns (x, y) forces due to the airfoil.
        Includes drag, lift, and the parasitic drag due to lift.
        If includeGravity is set, also includes apparent gravity (taking into
        account centrifugal effects).

        flightPitch: angle of velocity vector relative to the surface (degrees)
        AoA: angle of the airfoil relative to the velocity vector (degrees)
        """
        lift = self.liftForce(AoA, v, altitude, planet)
        drag = self.dragForce(AoA, v, altitude, planet)

        liftAngle = math.radians(flightPitch + 90)
        dragAngle = math.radians(flightPitch + 180)

        (liftX, liftY) = (lift * math.cos(liftAngle), lift * math.sin(liftAngle))
        (dragX, dragY) = (drag * math.cos(dragAngle), drag * math.sin(dragAngle))
        if not includeGravity:
            return (liftX + dragX, liftY + dragY)

        # now compute apparent gravity
        vOrbital = planet.orbitalVelocity(altitude)
        vSrfHorizontal = math.cos(math.radians(flightPitch)) * v
        vOrbHorizontal = vSrfHorizontal + planet.siderealSpeed(altitude)
        gravityAccel = planet.gravity(altitude)
        centrifugeAccel = vOrbHorizontal * vOrbHorizontal / (vOrbital * vOrbital)
        apparentGravityForce = (centrifugeAccel - gravityAccel) * self.mass
        return (liftX + dragX, liftY + dragY + apparentGravityForce)

    def liftDragRatio(self, AoAdegrees, v):
        lift = self.liftForce(AoAdegrees, v)
        drag = self.dragForce(AoAdegrees, v)
        return lift / drag

# Control surfaces include most of the winglets, the control surfaces, and
# the canards.  They have ModuleControlSurface.
# Most of the code is shared with the Winglet class, so I derive from its
# counterpart here.
class controlSurface(wing):
    def __init__(self, name, mass, lift, drag):
        wing.__init__(self, name, mass, lift, drag)

    @classmethod
    def deflectionLift(cls, AoAdegrees):
        cosUp = math.cos(math.radians(upDegrees(AoAdegrees)))
        return -cosUp

# Lifting Surfaces include many of the fuel tanks.
class liftingSurface(wing):
    def __init__(self, name, mass, lift, dragMax, dragMin):
        wing.__init__(self, name, mass, lift, 1)
        self.dragMin = dragMin
        self.dragMax = dragMax

    def dragCoeff(self, AoAdegrees):
        sinAoA = math.sin(math.radians(AoAdegrees))
        return sinAoA * self.dragMax + (1 - sinAoA) * self.dragMin


# For design, compare two lifting surfaces
def compare(srf1, srf2, AoAdegrees=30, quiet=False):
    """
    compare(srf1, srf2)

    Scale srf2 to have equal mass, and compare the lift and drag.
    Returns a tuple
            (#srf2 to make the mass of one srf1,
             ratio of scaled srf2:srf1 lift,
             ratio of scaled srf2:srf1 drag)

    Set AoAdegrees to the angle you want, 30 is the default.
    Set quiet=False to turn off printing.
    """
    def div(a, b):
        if a == 0 and b == 0:
            return 1
        if b == 0:
            return float('inf')
        return a/b

    num2 = srf1.mass / srf2.mass
    lift1 = srf1.liftForce(AoAdegrees, 1)
    lift2 = num2 * srf2.liftForce(AoAdegrees, 1)
    drag1 = srf1.dragForce(AoAdegrees, 1)
    drag2 = num2 * srf2.dragForce(AoAdegrees, 1)

    retval = (num2, div(lift2, lift1), div(drag2, drag1))
    if not quiet:
        print (("Use %g " + srf1.name + " for each " + srf2.name + ".\n"
               + "Lift is %gx ; drag is %gx.\n") % retval)
    return retval

def speedForLift(F, liftSurfaces = [], AoAdegrees = 0, altitude = 0, planet = planet.kerbin):
    """
    Return the speed we need to go in order to get a given lift force.
    Useful to figure out the liftoff speed of a plane.
    F is in kN
    liftSurfaces is a list of triples (n, liftType, AoA) with AoA being relative to your craft.
        - you can elide n to mean 1, or elide AoA to mean 0
    """
    coeff = 0
    for item in liftSurfaces:
        if item is wing:
            coeff += item.liftCoeff(AoAdegrees)
        else:
            try:
                (n, x, AoA) = item
                coeff += n * x.liftCoeff(AoA + AoAdegrees)
            except:
                try:
                    (n, x) = item
                    coeff += n * x.liftCoeff(AoAdegrees)
                except:
                    # last-ditch; this may throw
                    (x, AoA) = item
                    coeff += x.liftCoeff(AoA + AoAdegrees)
    # coeff is the lift coefficient of the whole craft
    return F / (coeff * planet.pressure(altitude))

# stock wings
sweptWing = wing("Swept Wings", 0.05, 1.6, 0.6)
tailFin = wing("Tail Fin", 0.02, 0.3, 0.5)
avt1 = wing("AV-T1 Winglet", 0.05, 0.3, 0.5)
wingConnector = wing("Wing Connector A", 0.1, 2, .5)
wingConnector2 = wing("Wing Connector B", 0.1, 2, .5)
wingConnector3 = wing("Wing Connector C", .05, 1, .25)
wingConnector4 = wing("Wing Connector D", .025, .5, .12)
wingConnector5 = wing("Wing Connector E", .025, .5, .12)
deltaWing = wing("Delta Wing", .1, 2, .6)
delta_small = wing("Small Delta Wing", .025, .5, .1)

wings = [
    sweptWing,
    tailFin,
    avt1,
    wingConnector,
    wingConnector2,
    wingConnector3,
    wingConnector4,
    wingConnector5,
    deltaWing,
    delta_small
]

# stock control surfaces
canard = controlSurface("Canard", 0.04, 0.7, 0.5)
avr8 = controlSurface("AV-R8", 0.02, 0.4, 0.5)
deltaDeluxe = controlSurface("Delta-Deluxe", 0.02, 0.7, 0.6)
elevon1 = controlSurface("Elevon 1", .02, .25, .25)
elevon2 = controlSurface("Elevon 2", .03, .33, .33)
elevon3 = controlSurface("Elevon 3", .04, .5, .5)
elevon4 = controlSurface("Elevon 4", .02, .25, .25)
elevon5 = controlSurface("Elevon 5", .05, .6, .6)

controlSurfaces = [
    canard,
    avr8,
    deltaDeluxe,
    elevon1,
    elevon2,
    elevon3,
    elevon4,
    elevon5,
]

# stock lifting surfaces
mk2_1m_Bicoupler = liftingSurface("Mk2 Bicoupler", .2, .2, .2, .1)
mk2_1m_AdapterLong = liftingSurface("Mk2 to 1.25m Adapter Long", .4, .6, .3, .1)
mk2SpacePlaneAdapter = liftingSurface("mk2 to 1.25m Adapter", .2, .3, .2, .1)
mk2CargoBayLarge = liftingSurface("Mk2 cargo (large)", .5, .6, .4, .1)
mk2CargoBaySmall = liftingSurface("Mk2 cargo (small)", .25, .3, .2, .1)
mk2Cockpit_Inline = liftingSurface("Mk2 Inline Cockpit", 2, .4, .4, .1)
mk2Cockpit_Standard = liftingSurface("Mk2 Cockpit", 2, .6, .1, .03)
mk2CrewCabin = liftingSurface("Mk2 Crew Cabin", 2, .3, .15, .1)
mk2DockingPort = liftingSurface("Mk2 Clamp-O-Tron", .3, .15, .3, .1)
mk2FuselageLong = liftingSurface("Mk2 Fuselage (long)", .5, .6, .3, .1)
mk2FuselageShort = liftingSurface("Mk2 Fuselage (short)", .25, .3, .15, .1)
wingStrake = liftingSurface("Wing Strake", .025, .75, .2, 0)
structuralWingA = liftingSurface("Structural Wing A", 0.05, 1, .3, 0)
structuralWingB = liftingSurface("Structural Wing B", 0.05, 1, .3, 0)
structuralWingC = liftingSurface("Structural Wing C", 0.025, .5, .15, 0)
structuralWingD = liftingSurface("Structural Wing D", 0.012, .25, .08, 0)
sweptWing12 = liftingSurface("Swept Wing Type A/B", 0.05, 1, .3, 0)

liftingSurfaces = [
    mk2_1m_Bicoupler,
    mk2_1m_AdapterLong,
    mk2SpacePlaneAdapter,
    mk2CargoBayLarge,
    mk2CargoBaySmall,
    mk2Cockpit_Inline,
    mk2Cockpit_Standard,
    mk2CrewCabin,
    mk2DockingPort,
    mk2FuselageLong,
    mk2FuselageShort,
    wingStrake,
    structuralWingA,
    structuralWingB,
    structuralWingC,
    structuralWingD,
    sweptWing12,
]
