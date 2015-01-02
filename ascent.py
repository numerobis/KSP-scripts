#! /usr/bin/python

from __future__ import division
import math
from physics import vector, g0
import planet

# Ascent simulation for an airless world:
# algorithm: thrust up just enough to cancel out gravity, taking into account
# the current horizontal speed, ignoring the atmosphere.

class simulation(object):
    def __init__(self, planet, thrust, Isp, mass, v = None, altitude = 0, deltaT = 0.03):
        self.t = 0
        self.p = vector(0, planet.radius + altitude)
        if v is None:
            self.v = vector(planet.siderealSpeed(altitude), 0)
        else:
            self.v = vector(*v)
        self.deltaT = deltaT
        self.mass = mass
        self.initialMass = mass
        self.thrust = thrust
        self.Ve = Isp * g0
        self.planet = planet

    @staticmethod
    def _polarToCartesian(r, theta):
        return vector(r * math.cos(theta), r * math.sin(theta))

    def _updatePosition(self, thrustAngle):
        """
        Return the new position, velocity, and acceleration
        given the thrust vector.

        Angle is in radians.
        """
        thrustAccel = self._polarToCartesian(self.thrust / self.mass, thrustAngle)

        a = thrustAccel.add(self.gravityVector())

        v1 = self.v.add(a.scale(self.deltaT))
        p1 = self.p.add(v1.scale(self.deltaT))

        return (p1, v1, a)

    def _isValidThrustAngle(self, angle):
        """
        Is the thrust angle valid given the altitude curve?

        Angle is in radians.
        """
        (p, v, a) = self._updatePosition(angle)
        r2New = p.sqrMagnitude()
        r2Old = self.p.sqrMagnitude()
        return r2New >= r2Old

    def optimizeThrustAngle(self, tolerance = 1e-6):
        """
        Find which way to point given current conditions such that
        our altitude next frame is the same as the altitude this
        frame.

        This is a numeric solve, hence the tolerance.  We return
        a solution that has altitude no lower than before.

        Returns the angle in radians.
        """
        # just binary search in the first quadrant for a solution that gets us
        # the same or larger altitude.
        verticalAngle = math.atan2(self.p[1], self.p[0])
        horizontalAngle = verticalAngle - math.pi / 2

        minAngle = horizontalAngle
        maxAngle = verticalAngle

        # Try horizontal and vertical first.
        if not self._isValidThrustAngle(maxAngle):
            raise "Insufficient thrust"
        if self._isValidThrustAngle(minAngle):
            return minAngle

        # Now, binary search.
        while maxAngle - minAngle > 1e-6:
            angle = 0.5 * (maxAngle + minAngle)
            if self._isValidThrustAngle(angle):
                maxAngle = angle
            else:
                minAngle = angle
        # err on the side of too steep
        return maxAngle

    def altitude(self):
        return self.p.L2() - self.planet.radius

    def horizonAngle(self):
        verticalAngle = math.atan2(self.p[1], self.p[0])
        horizonAngle = verticalAngle - math.pi / 2
        return horizonAngle

    def horizonVector(self):
        angle = self.horizonAngle()
        return vector(math.cos(angle), math.sin(angle))

    def localVelocity(self):
        """
        Return the velocity as (horizontal, vertical) in the local frame
        of reference.
        """
        horizon = self.horizonVector()
        up = vector(-horizon[1], horizon[0])
        return vector(self.v.dot(horizon), self.v.dot(up))

    def gravityVector(self):
        """
        Return the gravity vector in global coordinates.
        """
        horizon = self.horizonVector()
        down = vector(horizon[1], -horizon[0])
        g = self.planet.gravity(self.altitude())
        return down.scale(g)

    def achievedOrbit(self):
        vOrbit = math.sqrt(self.planet.mu / self.p.L2())
        return self.localVelocity()[0] >= vOrbit

    def deltaVBurned(self):
        return self.Ve * math.log(self.initialMass / self.mass)

    def step(self):
        """
        Simulate one timestep.

        Return false if we have burned all our mass or achieved orbit.
        """
        if self.mass <= 0:
            return False

        # TODO: first achieve Ap, then raise Pe
        if self.achievedOrbit():
            return False

        # If we run out of mass this frame, shorten the deltaT.
        # That affects all future computations, so we have to do this
        # first.
        mdot = self.thrust / self.Ve
        finalMass = self.mass - mdot * self.deltaT
        if finalMass < 0:
            self.deltaT = self.mass / mdot
            finalMass = 0

        # figure out which way to thrust
        angle = self.optimizeThrustAngle()

        (p1, v1, a) = self._updatePosition(angle)

        self.mass = finalMass
        self.p = p1
        self.v = v1
        self.t += self.deltaT

        altitude = self.p.L2() - self.planet.radius
        vLocal = self.localVelocity()
        angleDegrees = math.degrees(angle - self.horizonAngle())

        print ("t=%.1fs altitude %gm, speed (%g, %g) m/s, thrust angle %g degrees, mass %gt"
                % (self.t, altitude, vLocal[0], vLocal[1], angleDegrees, self.mass))
        return True

    def run(self):
        while self.step():
            pass


def simpleSim(planet, engine, m1, altitude = 0, nEngines = 1):
    """
    Run a simulation with nothing complicated going on:
    - given engine
    - given initial mass
    - assume we can burn all the mass
    - ignore atmosphere and terrain
    - optional altitude
    """
    sim = simulation(planet, engine.thrust * nEngines, engine.IspVac, m1, altitude=altitude)
    sim.run()
    if sim.achievedOrbit():
        print ("achieved orbit in %.1fs, burned %gt fuel (aka %g m/s deltaV)" 
                % (sim.t, m1 - sim.mass , sim.deltaVBurned()))
    else:
        print ("failed to achieve orbit, need more than %gt fuel"
                % m1)
