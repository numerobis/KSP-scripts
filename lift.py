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

    @classmethod
    def deflectionDrag(cls, AoAdegrees):
        return math.sin(math.radians(AoAdegrees))

    def _liftFactor(self, AoAdegrees, altitude = 0, planet = planet.kerbin):
        return self.deflectionLift(AoAdegrees) * planet.pressure(altitude) * self.lift

    def liftForce(self, AoAdegrees, v, altitude = 0, planet = planet.kerbin):
        return v * self._liftFactor(AoAdegrees, altitude, planet)

    def speedForLift(self, AoAdegrees, F, altitude = 0, planet = planet.kerbin):
        return F / self._liftFactor(AoAdegrees, altitude, planet)

    def dragForce(self, AoAdegrees, v, altitude = 0, planet = planet.kerbin):
        Cd = self.deflectionDrag(AoAdegrees) * self.drag
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

# stock wings
deltaWing = wing(0.07, 1.9, 0.6)
sweptWing = wing(0.05, 1.6, 0.6)
tailFin = wing(0.02, 0.3, 0.5)
wingConnector = wing(0.05, 1, 0.4)
structuralWing = wing(0.05, 1, 0.4)
avt1 = wing(0.05, 0.3, 0.5)

# stock control surfaces
smallControlSurface = controlSurface(0.01, 0.5, 0.5)
largeControlSurface = controlSurface(0.04, 0.7, 0.5)
canard = controlSurface(0.04, 0.7, 0.5)
avr8 = controlSurface(0.02, 0.4, 0.5)
deltaDeluxe = controlSurface(0.02, 0.7, 0.6)
