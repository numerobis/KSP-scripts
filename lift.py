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
    def __init__(self, mass, lift, drag):
        self.mass = mass
        self.lift = lift
        self.drag = drag

    @classmethod
    def deflectionLift(cls, AoAdegrees):
        cosAoA = math.cos(math.radians(AoAdegrees))
        cosUp = math.cos(math.radians(upDegrees(AoAdegrees)))
        return -cosUp * (1.0 - abs(cosUp)) * cosAoA

    def dragCoeff(cls, AoAdegrees):
        return math.sin(math.radians(AoAdegrees)) * self.drag

    def _liftFactor(self, AoAdegrees, altitude = 0, planet = planet.kerbin):
        return self.deflectionLift(AoAdegrees) * planet.pressure(altitude) * self.lift

    def liftForce(self, AoAdegrees, v, altitude = 0, planet = planet.kerbin):
        return v * self._liftFactor(AoAdegrees, altitude, planet)

    def speedForLift(self, AoAdegrees, F, altitude = 0, planet = planet.kerbin):
        return F / self._liftFactor(AoAdegrees, altitude, planet)

    def dragForce(self, AoAdegrees, v, altitude = 0, planet = planet.kerbin):
        Cd = self.dragCoeff(AoAdegrees)
        return planet.dragForce(altitude, v, self.mass, Cd)

    def forceVector(self, flightPitch, AoA, v, altitude = 0, planet = planet.kerbin):
        """
        Returns (x, y) forces due to the airfoil.
        Includes drag, lift, and the parasitic drag due to lift.

        flightPitch: angle of velocity vector relative to the surface (degrees)
        AoA: angle of the airfoil relative to the velocity vector (degrees)
        """
        lift = self.liftForce(AoA, v, altitude, planet)
        drag = self.dragForce(AoA, v, altitude, planet)

        liftAngle = math.radians(flightPitch + 90)
        dragAngle = math.radians(flightPitch + 180)

        (liftX, liftY) = (lift * math.cos(liftAngle), lift * math.sin(liftAngle))
        (dragX, dragY) = (drag * math.cos(dragAngle), drag * math.sin(dragAngle))
        return (liftX + dragX, liftY + dragY)


# Control surfaces include most of the winglets, the control surfaces, and
# the canards.  They have ModuleControlSurface.
# Most of the code is shared with the Winglet class, so I derive from its
# counterpart here.
class controlSurface(wing):
    def __init__(self, mass, lift, drag):
        wing.__init__(self, mass, lift, drag)

    @classmethod
    def deflectionLift(cls, AoAdegrees):
        cosUp = math.cos(math.radians(upDegrees(AoAdegrees)))
        return -cosUp

# Lifting Surfaces include many of the fuel tanks.
class liftingSurface(wing):
    def __init__(self, mass, lift, dragMax, dragMin):
        wing.__init__(self, mass, lift, 1)
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
        print (("Use %g of the latter for each of the former.\n"
               + "Lift is %gx ; drag is %gx.\n") % retval)
    return retval

# stock wings
sweptWing = wing(0.05, 1.6, 0.6)
tailFin = wing(0.02, 0.3, 0.5)
avt1 = wing(0.05, 0.3, 0.5)
wingConnector = wing(0.1, 2, .5)
wingConnector2 = wing(0.1, 2, .5)
wingConnector3 = wing(.05, 1, .25)
wingConnector4 = wing(.025, .5, .12)
wingConnector5 = wing(.025, .5, .12)
deltaWing = wing(.1, 2, .6)
delta_small = wing(.025, .5, .1)

# stock control surfaces
smallControlSurface = controlSurface(0.01, 0.5, 0.5)
largeControlSurface = controlSurface(0.04, 0.7, 0.5)
canard = controlSurface(0.04, 0.7, 0.5)
avr8 = controlSurface(0.02, 0.4, 0.5)
deltaDeluxe = controlSurface(0.02, 0.7, 0.6)
elevon1 = controlSurface(.02, .25, .25)
elevon2 = controlSurface(.03, .33, .33)
elevon3 = controlSurface(.04, .5, .5)
elevon4 = controlSurface(.02, .25, .25)
elevon5 = controlSurface(.05, .6, .6)

# stock lifting surfaces
mk2_1m_Bicoupler = liftingSurface(.2, .2, .2, .1)
mk2_1m_AdapterLong = liftingSurface(.4, .6, .3, .1)
mk2SpacePlaneAdapter = liftingSurface(.2, .3, .2, .1)
mk2CargoBayLarge = liftingSurface(.5, .6, .4, .1)
mk2CargoBaySmall = liftingSurface(.25, .3, .2, .1)
mk2Cockpit_Inline = liftingSurface(2, .4, .4, .1)
mk2Cockpit_Standard = liftingSurface(2, .6, .1, .03)
mk2CrewCabin = liftingSurface(2, .3, .15, .1)
mk2DockingPort = liftingSurface(.3, .15, .3, .1)
mk2FuselageLong = liftingSurface(.5, .6, .3, .1)
mk2FuselageShort = liftingSurface(.25, .3, .15, .1)
wingStrake = liftingSurface(.025, .75, .2, 0)
structuralWingA = liftingSurface(0.05, 1, .3, 0)
structuralWingB = liftingSurface(0.05, 1, .3, 0)
structuralWingC = liftingSurface(0.025, .5, .15, 0)
structuralWingD = liftingSurface(0.012, .25, .08, 0)
sweptWing1 = liftingSurface(0.05, 1, .3, 0)
sweptWing2 = liftingSurface(0.05, 1, .3, 0)
