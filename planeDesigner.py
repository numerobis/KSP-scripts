from __future__ import division

import math

import engine
import jets
import lift
import planet
from physics import vectoradd, vectorscale

def forces(parts, mass, flightPitch, AoA, v, altitude = 0, planet = planet.kerbin, deltaT = 0.04):
    """
    Given a description of the plane, its total mass, 
    the flight path pitch (the angle of prograde with the surface),
    the angle of attack (the angle of pitch with prograde),
    the speed,
    the altitude and the planet, 

    Return the forces assuming we're at the max thrust we can muster at
    this altitude and controls are flat.
    """
    vlift = [0, 0]
    vdrag = [0, 0]
    maxthrust = [0, 0]
    airUse = 0
    airGet = 0
    accountedMass = 0

    pitch = flightPitch + AoA
    jetOptions = jets.standardoptions(planet, deltaT)

    flightPitchRad = math.radians(flightPitch)
    (flightPitchCos, flightPitchSin) = (math.cos(flightPitchRad), math.sin(flightPitchRad))

    def standardDrag(mass):
        f = planet.dragForce(altitude, v, mass, 0.2)
        return (-f * flightPitchCos, -f * flightPitchSin)

    def liftVector(typ, angle):
        if isinstance(typ, lift.wing):
            f = typ.liftForce(AoA + angle, v, altitude, planet)
            # lift is 90 degrees from prograde
            return (-f * flightPitchSin, f * flightPitchCos)
        else:
            return (0,0)

    def dragVector(typ, angle):
        if isinstance(typ, lift.wing):
            f = typ.dragForce(AoA + angle, v, altitude, planet)
            return (-f * flightPitchCos, -f * flightPitchSin)
        elif isinstance(typ, jets.intake):
            f = typ.drag(altitude, v, AoA + angle, options = jetOptions)
            return (-f * flightPitchCos, -f * flightPitchSin)
        else:
            return standardDrag(typ.mass)

    def thrustVector(typ, angle):
        if isinstance(typ, jets.jetengine):
            t = typ.thrust(v)
        elif isinstance(typ, engine.engine):
            t = typ.thrust
        else:
            t = 0
        enginePitch = math.radians(angle + pitch)
        return (math.cos(enginePitch) * t, math.sin(enginePitch) * t)

    def airRequired(typ):
        if isinstance(typ, jets.jetengine):
            return typ.airRequired(altitude, options = jetOptions)
        else:
            return 0

    def airProvided(typ, angle):
        if isinstance(typ, jets.intake):
            return typ.intakeair(altitude, v, AoA + angle, options = jetOptions)
        else:
            return 0

    for part in parts:
        try:
            (typ, n, angle) = part
        except:
            try:
                (typ, n) = part
                angle = 0
            except:
                typ = part
                n = 1
                angle = 0
        vlift = vectoradd(vlift, vectorscale(n, liftVector(typ, angle)))
        vdrag = vectoradd(vdrag, vectorscale(n, dragVector(typ, angle)))
        maxthrust = vectoradd(maxthrust, vectorscale(n, thrustVector(typ, angle)))
        airUse += n * airRequired(typ)
        airGet += n * airProvided(typ, angle)
        accountedMass += typ.mass

    # account for throttle setting due to air
    if airUse <= airGet:
        throttle = 1
    else:
        throttle = airGet / airUse
    thrust = vectorscale(throttle, maxthrust)

    # account for drag on unaccounted mass (at Cd = 0.2)
    unaccountedMass = mass - accountedMass
    if unaccountedMass > 0:
        vdrag = vectoradd(vdrag, standardDrag(unaccountedMass))

    gravity = -planet.gravity(altitude) * mass
    net = [0, gravity]
    net = vectoradd(net, vlift)
    net = vectoradd(net, vdrag)
    net = vectoradd(net, thrust)

    print ("throttle: %.1f%%" % (throttle * 100))
    print ("lift:    (%7.2f, %7.2f) kN" % (vlift[0], vlift[1]))
    print ("drag:    (%7.2f, %7.2f) kN" % (vdrag[0], vdrag[1]))
    print ("thrust:  (%7.2f, %7.2f) kN" % (thrust[0], thrust[1]))
    print ("gravity: (       , %7.2f) kN" % gravity)
    print ("net:     (%7.2f, %7.2f) kN" % (net[0], net[1]))

    return net
