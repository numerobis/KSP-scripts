# an SSTO to Laythe
import engine
import lift
import planet
import math
import jets

def stage(payload, dV, Isp, beta, nextStageFn):
    # The tanks we add here affect the lower stages.
    tanks = 0
    payloadToOrbit = payload
    for i in range(0,3):
        (prop, tanks) = engine.burnMass(dV, Isp, payloadToOrbit, beta)
        (payloadToOrbit, beta4tanks, beta8tanks) = nextStageFn(payload + tanks)
    if beta == 4:
         beta4tanks += tanks
    else:
         beta8tanks += tanks
    return (payloadToOrbit + prop, beta4tanks, beta8tanks)

def lastStage(payload):
    # Transfer from laythe to kerbin on year 2, day 264, costs
    # 1556 Jool-relative; that's 1050 from a low Laythe orbit.
    # Exit at 3am.  The window is large.
    #
    # There's potential for a much cheaper escape via an assist at Tylo
    # apparently.
    #
    (prop, beta8tanks) = engine.burnMass(1100, 800, payload)
    return (payload + beta8tanks + prop, 0, beta8tanks)

def laytheUpStage(payload):
    return stage(payload, 4000, 10000, 4, lastStage)

def lkoToJool(payload):
    # Transfer from kerbin 100km to jool on year 1, day 171 costs
    # 1996 m/s from Kerbin.  Exit at 9:36. Arrive year 2, day 208.
    # Potential for a much cheaper transfer via two hops at Eve.
    return stage(payload, 2050, 800, 8, laytheUpStage)

def toLKO(payload):
    return stage(payload, 6000, 10000, 4, lkoToJool)

#payload = sum( (
#            3.5,
#            3 * 2.5,
#            0.8,
#            4 * 1.2,
#            196 * 0.01,
#            2.25,
#            0.4,
#            77 * 0.07,
#            ))


nintakes = 0
ncontrol = 1
nturbo = 2
takeoffAoA = 20 # pull back to 20 degree pitch
takeoffSpeed = 60
runwayAltitude = 70
liftSurfaceAoA = 60

for i in range(1, 10):
    payload = sum( (
                0.08, # nose cone science
                0.05, # nose cone adapter
                3.5, # cockpit
                2 * 0.01, # cockpit canards
                3.5, # lab
                5 * 0.005, # antenna and sensors
                0.2, # mat bay
                0.15, # goo
                0.08, # RTG
                3 * 2.5, # habs
                2.25, # LV-N
                4 * 0.01, # tail
                6 * 0.05, # wing connectors
                nturbo * 1.2, # turbojets
                nintakes * 0.01, # intakes
                ncontrol * 0.01, # small control
                ))
    (totalMass, beta4, beta8) = toLKO(payload)

    needLiftForce = totalMass * planet.kerbin.gravity(runwayAltitude)
    perWinglet = lift.smallControlSurface.lift(10 + takeoffAoA, takeoffSpeed, runwayAltitude)
    newNcontrol = math.ceil(needLiftForce / perWinglet)

    newNintakes = nintakes + jets.intakesForSpeed(
            ( (jets.turbojet, nturbo),
              (jets.ramAirIntake, nintakes) ),
            totalMass,
            2175) # I want to achieve circular orbit at 40km

    converged = (newNcontrol == ncontrol and newNintakes == nintakes)
    ncontrol = newNcontrol
    nintakes = newNintakes

    print ("iteration %d" % i)
    print ("payload %g needs %gt bipropellent tanks and %gt jet tanks, total mass %g" % (payload, beta8, beta4, totalMass))
    print ("suggest %d small control surfaces" % ncontrol)
    print ("suggest %d intakes" % nintakes)
    if converged:
        print "converged"
        break
