# KSP Atmospheric ascent planning module
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


from __future__ import division

import math
from math import cos, sin, sqrt

import physics
from physics import L2


def angles(p, v):
    """
    Given a position and velocity in a cartesian coordinate frame centered on a
    planet, return the following angles:
    1. The angle p makes with the horizontal.
    2. The angle v makes with the horizontal.
    3. The angle v makes with the surface.
    """
    # psi is the angle of our position relative to the core, in the
    # initial coordinate system (90 is straight up from the core).
    psi = math.atan2(p[1], p[0])

    # theta is the angle we're going in the initial coordinate
    # system.  0 = flat, 90 = vertical.  At launch with speed 0,
    # match psi.
    if sum(x*x for x in v) == 0:
        theta = psi
    else:
        theta = math.atan2(v[1], v[0])

    # thetaSurf is the angle relative to the surface
    thetaSurf = theta + math.pi/2 - psi

    return (psi, theta, thetaSurf)



class BadFlightPlanException(Exception): pass

class climbSlope(object):
    def __init__(self,
        planet,
        orbitAltitude = None,
        gravityTurnStart = None,
        gravityTurnEnd = None,
        targetVelocity = None,
        initialAltitude = 0,
        initialVelocity = None,
        launchInclination = 0,
        acceleration = None,
        timestep = 1):
        """
        Compute a climb slope for exiting the atmosphere and achieving orbit.

        The result has a series of timesteps that can be interpolated with functions below.

        By default the orbital altitude is set to 1km above the top of the
        atmosphere.  This can be changed.

        By default, we start on the surface at the equator and we're going due
        East, so we take advantage of the sidereal rotation.

        This is not presently optimal: when should we start and end the gravity
        turn, and what should the curve look like?  Currently the assumption is
        that we climb through 85% of the atmosphere then start turning until TOA,
        with an angle that varies linearly with altitude.  The start and end points
        of the turn can be changed, the curve cannot.

        There's no model of variable acceleration.  It would be easy
        to add in since we are using numerical integration; you'd want to
        take in a function from (deltaV, altitude, velocity) ->
        acceleration to model jets.  By default, acceleration is assumed to be 2.2g,
        which is good low down (except the first ~100m) but too little at altitude.

        Specify the timestep to change the accuracy; by default we
        use a timestep of 1s.  It is a fixed timestep.  Smaller timesteps
        improve simulation accuracy (until we start to build up numerical error),
        but take longer.  This uses Euler integration because I'm too lazy
        to code up Runge-Kutta.
        """
        # The math:
        # At each timestep, we compute:
        # - the current gravity
        # - the current drag
        # - the current terminal velocity
        # - the current velocity (theta) and thrust (phi) angles
        # We control the thrust:
        # - maximum thrust subject to:
        # - don't exceed terminal velocity
        # - don't exceed apoapsis
        # Then we update:
        # - altitude goes up by the current speed's y component.
        # - velocity increases according to thrust
        # - deltaV goes up by the thrust applied
        # - time goes up by the timestep
        # - drag loss goes up
        class ClimbPoint(object):
            def __init__(self, alt, v, thrust, t, dV, dragLoss, thrustLimited):
                self.altitude = alt
                self.velocity = v
                self.thrust   = thrust
                self.time     = t
                self.deltaV   = dV
                self.dragLoss = dragLoss
                self.thrustLimited = thrustLimited

            def __str__(self):
                theta = math.atan2(self.velocity[1], self.velocity[0])
                return (
                    "%g s: %gm altitude, %g m/s at pitch %g; %gm/s deltaV, %g drag, thrust %g m/s^2%s"
                    % (self.time, self.altitude,
                        L2(self.velocity), math.degrees(theta),
                        self.deltaV,
                        self.dragLoss,
                        self.thrust,
                        "" if self.thrustLimited == 0 else
                                (" needs %g m/s^2 more thrust" % self.thrustLimited)
                    ))

        self.planet = planet
        self.launchInclination = launchInclination # used to estimate circularization

        if initialVelocity is None:
            v = [0,0]
        else:
            v = list(initialVelocity)

        # p and v are the current position and velocity; we update these.
        # We use a coordinate system centered at the core of the planet,
        # aligned so that the initial altitude forms the y axis.
        # Angles theta and phi refer to this coordinate system.
        p = [0, initialAltitude + planet.radius] # (m,m)
        h = L2(p)                                # m, update when p is updated
        dV = 0                                   # m/s
        t = 0                                    # s
        dragLoss = 0                             # m/s
        if acceleration is None:
            a_thrust = 2.2 * planet.gravity(initialAltitude) # m/s^2
        else:
            a_thrust = acceleration

        TOA = planet.topOfAtmosphere() # m
        if orbitAltitude is None:
            orbitAltitude = TOA + 1000
        self.orbitAltitude = orbitAltitude
        v_orbit = planet.orbitalVelocity(orbitAltitude)

        if gravityTurnStart is None:
            # No idea what's optimal here...
            gravityTurnStart = planet.altitude(planet.datumPressure / 8)

        if gravityTurnEnd is None:
            # No idea what's optimal here...
            gravityTurnEnd = min(TOA + 1000, orbitAltitude)
            assert gravityTurnEnd > gravityTurnStart

        #print ("Aiming for orbit %d m (%d m/s), start turn at %d, end at %d" %
        #    (orbitAltitude, v_orbit, gravityTurnStart, gravityTurnEnd))

        # Numerical integration time: every timestep until we reach apoapsis,
        # decide how much to burn and update our position.
        targetApoapsis = orbitAltitude + planet.radius
        climbSlope = []             # ClimbPoint list
        while h < targetApoapsis and t < 1000:
            if h < planet.radius:
                # We crashed!  Bring more thrust next time.
                raise BadFlightPlanException

            alt = h - planet.radius

            # get:
            # the angle of our position relative to coordinate horizontal
            # the angle of velocity relative to coordinate horizontal
            # the angle of velocity relative to surface horizontal
            (psi, theta, thetaSurf) = angles(p, v)

            # phiSurf is the angle of thrust relative to the surface.
            # Straight up before the gravity turn, straight sideways after,
            # and interpolate linearly during the turn.
            if alt <= gravityTurnStart:
                phiSurf = 0
            elif alt >= gravityTurnEnd:
                phiSurf = - (math.pi / 2)
            else:
                # What fraction of the turn have we done?
                ratio = (alt - gravityTurnStart) / (gravityTurnEnd - gravityTurnStart)
                # The more we did, the closer we want to be to pi/2 from
                # the straight-up orientation.
                phiSurf = - (ratio * math.pi / 2)
            # phi is the angle of thrust in the original coordinate frame.
            phi = phiSurf + psi

            # Compute the gravity and drag.
            a_grav = planet.gravity(alt)
            a_drag = planet.drag(alt, L2(v))

            # Terminal velocity for our altitude:
            #print ("drag %g m/s^2, grav %g m/s^2, term %g m/s" % (a_drag, a_grav, v_term))

            # Compute the velocity after the given amount of thrust.
            # First, compute the loss on both axes.
            loss = (
                - (a_drag * cos(theta) + a_grav * cos(psi)),
                - (a_drag * sin(theta) + a_grav * sin(psi))
            )
            v_nothrust = [ v[i] + loss[i] * timestep for i in range(2) ]
            def thrustResult(thrust):
                return (
                    v_nothrust[0] + thrust * cos(phi) * timestep,
                    v_nothrust[1] + thrust * sin(phi) * timestep,
                )

            # Compute the acceleration that gets us to terminal velocity in the
            # direction of thrust.
            # Let X be (cos phi, sin phi) => transform to the direction of thrust.
            # Obviously X^T X === X^2 = 1
            # The velocity after alpha thrust is:
            #       v_no + alpha X dt
            # How much of that is in the direction of thrust?  Dot with X:
            #       (v_no + alpha X dt)^T X
            # We want the dot product to be at most v_term in absolute value
            # (we don't care which direction it's going).  In fact, we want it to be equal,
            # since we want to go as fast as possible:
            #       | (v_no + alpha X dt)^T X | = v_term
            # Equivalently, just square both sides:
            #       [(v_no + alpha X dt)^T X]^2 = v_term^2
            # Now solve for alpha.
            #       (v_no^T X)^2 + 2 alpha dt v_no^T X + alpha^2 dt^2 X^2 = v_term^2
            #       dt^2 alpha^2 + 2 dt (v_no^T X) alpha + (v_no^T X)^2 - v_term^2 = 0
            # Makes a quadratic equation.
            #
            # Return the bigger solution, or zero if both solutions are negative.
            def findTerminalThrust():
                # above the top of the atmosphere, there is no limit!
                if alt >= TOA: return 1e30
                v_term = planet.terminalVelocity(alt)
                a = timestep * timestep
                v_noTX = v_nothrust[0] * cos(phi) + v_nothrust[1] * sin(phi)
                b = 2 * timestep * v_noTX
                c = v_noTX - v_term * v_term
                soln = physics.quadratic(a,b,c)
                return max(soln)

            # Binary search the thrust to achieve the apoapsis.
            # I'm sure this could be worked out analytically.
            # Use the terminal velocity thrust to reduce our search space.
            def findApoapsisThrust(thrust_term):

                # The most we could want to speed up is enough to immediately
                # get our momentum to match what it will when we have a
                # circular orbit:
                # Generally,
                #       P = r.v
                # now: r.v = (|v| cos theta, |v| sin theta) . p
                #       |v| = P/p.(cos theta, sin theta)
                # with P unknown.
                # at the circular orbit we are targetting, P = apoapsis * v_orbit.
                # So the most speed we want now is:
                #   |v| = v_orbit * targetApoapsis / (p[0] cos theta + p[1] sin theta)
                # And we want to reach that speed in one timestep.
                # This is an upper bound: if we thrust that much now, we'll
                # send our apoapsis much higher than the target, since we're
                # generally not currently at the present apoapsis.
                v_targetMax = (v_orbit * targetApoapsis /
                               (p[0]*cos(thetaSurf) + p[1]*sin(thetaSurf)))
                thrust_targetMax = (v_targetMax - L2(v_nothrust)) / timestep

                thrust_lo = 0
                thrust_hi = min(thrust_targetMax, thrust_term)
                # while we are off by more than 1mm/s, binary search
                # we waste a bit of time searching for really big solutions.
                while thrust_hi - thrust_lo > 1e-3:
                    thrust = (thrust_hi + thrust_lo) / 2
                    v_next = thrustResult(thrust)
                    (apoapsis, _) = planet.determineOrbit2(p, v_next)

                    # if apoapsis is too high, or negative (i.e. hyperbolic), reduce thrust
                    if apoapsis > targetApoapsis or apoapsis < 0:
                        thrust_hi = thrust
                    else:
                        thrust_lo = thrust
                return thrust_lo # or hi, it's only 1mm/s difference

            thrust_term = findTerminalThrust()
            thrust_apo = findApoapsisThrust(thrust_term)

            thrust = min(thrust_term, thrust_apo)
            #print ("thrust: %g to reach term, %g to reach apo" % (thrust_term, thrust_apo))
            if thrust < a_thrust:
                thrustLimit = 0
                if thrust < 0:
                    thrust = 0
            else:
                thrustLimit = thrust - a_thrust
                thrust = a_thrust

            # Update everything!
            for i in (0,1): p[i] += v[i] * timestep
            h = L2(p)
            v = thrustResult(thrust)
            dV += thrust * timestep
            dragLoss += a_drag * timestep
            t += timestep
            climbSlope.append(ClimbPoint(alt, v, thrust, t, dV, dragLoss, thrustLimit))

        if t == 1000:
            # Timed out...
            raise BadFlightPlanException

        self._climbSlope = climbSlope
        self.planet = planet
        self.orbit = orbitAltitude

    def _bsearch(self, attrname, query):
        """
        Return the index that has the biggest key that is no bigger
        than the query key.  Must be something that grows monotonically
        (deltaV, altitude, losses, or time).
        Raises a KeyError if the query is outside the range we simulated.
        """
        lo = 0
        data = self._climbSlope
        hi = len(data)
        while lo < hi - 1:
            mid = (lo + hi) // 2
            val = getattr(data[mid], attrname)
            if val == query:
                return mid
            elif val < query:
                lo = mid
            else:
                hi = mid
        if lo == len(data) - 1:
            # Not found.
            raise KeyError
        else:
            return lo

    def _interpolate(self, keyname, query, valuename):
        """
        Search the climb slope for the query, which will generally lie
        between two data points, and linearly interpolate between them.
        The key must be sorted (this means deltaV, altitude, losses, or
        time).
        Raises a KeyError if the query is outside the range we simulated.
        """
        index = self._bsearch(keyname, query)
        before = self._climbSlope[index]
        beforekey = getattr(before, keyname)
        if beforekey == query:
            return getattr(before, valuename)
        if query < beforekey:
            raise KeyError

        after  = self._climbSlope[index + 1]
        afterkey = getattr(after, keyname)
        assert afterkey > query
        assert afterkey > beforekey

        beforevalue = getattr(before, valuename)
        aftervalue = getattr(after, valuename)

        ratio = (query - beforekey) / (afterkey - beforekey)
        return ratio * aftervalue + (1-ratio) * beforevalue


    def deltaVToAltitude(self, altitude):
        """
        deltaV we spent up and including the given altitude
        At 0 meters, we include the cost of the first second of burn.

        Raises a KeyError if the query is outside the range we simulated.
        """
        return self._interpolate("altitude", altitude, "deltaV")

    def deltaV(self):
        """
        Tally up how much total deltaV required to get to orbit.
        - deltaV we spent up to apoapsis
        - deltaV we need in order to circularize
        - take account of sidereal rotation depending on the launch angle
        """
        v_orbit = self.planet.orbitalVelocity(self.orbitAltitude)
        v_last = L2(self._climbSlope[-1].velocity)
        v_last += cos(self.launchInclination) * self.planet.siderealRotationSpeed
        dV_circ = v_orbit - v_last
        return self._climbSlope[-1].deltaV + dV_circ

    def deltaVBetween(self, altitude0, altitude1):
        """
        Return the deltaV to transition between altitude 0 and altitude 1
        along the climb slope.
        """
        return self.deltaVToAltitude(altitude1) - self.deltaVToAltitude(altitude0)

    def dragLosses(self):
        """
        Tally up how much we lost fighting aerodynamic drag to get up to
        apoapsis.
        """
        return self._climbSlope[-1].dragLoss

    def altitudeAtDeltaV(self, deltaV, default = KeyError):
        """
        Return the altitude after a given amount of deltaV.

        Raises a KeyError if the query is outside the range we simulated
        (i.e. it corresponds to deltaV that exceeds what we spent to get to
        apoapsis).
        """
        if deltaV <= self._climbSlope[0].deltaV:
            return self._climbSlope[0].altitude
        try:
            return self._interpolate("deltaV", deltaV, "altitude")
        except KeyError:
            if default == KeyError:
                raise
            else:
                return default

    def __str__(self):
        return "\n".join( (str(p) for p in self._climbSlope) )

