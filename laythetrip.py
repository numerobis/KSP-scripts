# an SSTO to Laythe
import engine
import lift
import planet
import math
import jets

class Phase(object):
    def __init__(self, name, deltaV, motor):
        self._name = name
        self._deltaV = deltaV
        if isinstance(motor, engine.engine):
            self._Isp = motor.IspVac
            self._beta = 8
        else:
            # it's a jet!
            self._Isp = 10000
            self._beta = 4

    def evaluate(self, payload):
        (prop, tanks) = engine.burnMass(self._deltaV, self._Isp, 
                payload, self._beta)
        if self._beta == 4:
            return (prop, tanks, 0)
        else:
            return (prop, 0, tanks)

lvn = engine.getEngine('lv-n')
phases = (
        Phase("KSC-LKO", 6000, jets.turbojet),
        Phase("LKO circ", 30, lvn),
        Phase("LKO-Jool", 2100, lvn),
        Phase("Laythe-LLO", 3000, jets.turbojet),
        Phase("LLO-Kerbin", 1100, lvn)
)

class Stage(object):
    def __init__(self, payload, prop, jetTanks, rocketTanks):
        self.jetTanks = jetTanks
        self.rocketTanks = rocketTanks
        self.mass = sum((payload, prop, jetTanks, rocketTanks))

def evaluateStages(payload, index):
    (prop, jetTanks, rocketTanks) = phases[index].evaluate(payload)

    if index + 1 == len(phases):
        return [ Stage(payload, prop, jetTanks, rocketTanks) ]
    else:
        for _ in xrange(4):
            payloadAndTanks = payload + jetTanks + rocketTanks
            laterStages = evaluateStages(payloadAndTanks, index + 1)
            (prop, jetTanks, rocketTanks) = phases[index].evaluate(laterStages[-1].mass)
        laterStages.append(Stage(laterStages[-1].mass, prop, jetTanks, rocketTanks))
        return laterStages

def reportLiftoffSpeed(name, mass, planet, altitude):
    liftForceNeeded = planet.gravity(altitude) * mass
    takeoffLift = liftForceNeeded  / nlift
    liftoffSpeed = lift.smallControlSurface.speedForLift(90, takeoffLift, altitude, planet)
    print("%s speed %g m/s" % (name, liftoffSpeed))


nintakes = 0
nlift = 0
nturbo = 2
takeoffSpeed = 60
runwayAltitude = 70

fixedPayload = sum( (
            0.08, # nose cone science
            0.05, # nose cone adapter
            3.5, # cockpit
            3.5, # lab
            6 * 0.005, # antenna and sensors and ladder
            0.2, # mat bay
            0.15, # goo
            5 * 0.005, # redundant fixed solar panels (front/back/bottom/top-left/top-right)
            3 * 2.5, # habs
            lvn.mass,
            2 * 0.05, # rover wheels
            (2 + # canard
             2 + # ailerons
             4)  # elevator and rudder
            * lift.smallControlSurface.mass
            ))

for i in range(1, 10):
    payload = sum( (
                fixedPayload,
                nturbo * jets.turbojet.mass, # turbojets
                nintakes * jets.ramAirIntake.mass, # intakes
                nlift * lift.smallControlSurface.mass, # small control surfaces tilted 90 down
                ))
    stages = list(reversed(evaluateStages(payload, 0)))
    launchMass = stages[0].mass

    needLiftForce = launchMass * planet.kerbin.gravity(runwayAltitude)
    liftPerSCS = lift.smallControlSurface.liftForce(90, takeoffSpeed, runwayAltitude)
    newNcontrol = math.ceil(needLiftForce / liftPerSCS)


    parts = ((jets.turbojet, nturbo), (jets.ramAirIntake, nintakes))
    newNintakes = nintakes + jets.intakesForSpeed(
            parts,
            launchMass,
            2175) # I want to achieve circular orbit at 40km

    converged = (newNcontrol == nlift and newNintakes == nintakes)
    nlift = newNcontrol
    nintakes = newNintakes

    print ("iteration %d mass %g" % (i, launchMass))
    if converged:
        print "converged"
        break

print ("payload %g, launch mass %g" % (fixedPayload, launchMass))
print ("%d turbojets, %d intakes, %d small control surfaces pointed down"
        % (nturbo, nintakes, nlift))

jetTanks = sum( [ x.jetTanks for x in stages ] )
rocketTanks = sum( [ x.rocketTanks for x in stages ] )
print ("%gt rocket fuel and %gt jet fuel" % (rocketTanks * 8, jetTanks * 4))
print ("%g U LiquidFuel, %g U Oxidizer" % (
            (rocketTanks * 8 * 0.9 + jetTanks * 4) * 200,
            (rocketTanks * 8 * 1.1) * 200))
for (phase, stage) in zip(phases, stages):
    name = phase._name
    if phase._beta == 4:
        fuelName = 'jet'
        tankMass = stage.jetTanks
    else:
        fuelName = 'rocket'
        tankMass = stage.rocketTanks
    startMass = stage.mass
    fuelMass = phase._beta * tankMass
    fuelUnits = fuelMass * 200
    print ("%10s:  %7.3f t %s fuel (%g units)" %
            (name, fuelMass, fuelName, fuelUnits))


print ""
reportLiftoffSpeed("KSC liftoff", launchMass, planet.kerbin, runwayAltitude)
laytheLandingMass = stages[2].mass
reportLiftoffSpeed("Laythe landing & liftoff", stages[2].mass, planet.laythe, 200)
kerbinLandingMass = fixedPayload + jetTanks + rocketTanks
reportLiftoffSpeed("KSC landing", kerbinLandingMass, planet.kerbin, runwayAltitude)
