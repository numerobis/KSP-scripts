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

class winglet(object):
    def __init__(self, mass, lift):
        self._mass = mass
        self._lift = lift

    def liftForce(self, AoAdegrees, v, altitude = 0, planet = planet.kerbin):
        cosAoA = math.cos(math.radians(AoAdegrees))
        cosUp = math.cos(math.radians(AoAdegrees + 90))
        liftMag = -v * cosUp * (1.0 - abs(cosUp)) * cosAoA * self._lift
        return liftMag * planet.pressure(altitude)

smallControlSurface = winglet(0.01, 0.5)
