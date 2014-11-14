# an SSTO to Laythe
import engine
import lift
import planet
import math
import jets

# what kind of rocket, intake, and lift surface will you use?
rocket = engine.getEngine('lv-909')
intake = jets.radialIntake
liftSurface = lift.deltaDeluxe
jetBeta = 6 # fuel : tank ratio, 6 for the Mk2 tanks
rocketBeta = 8 # fuel : tank, 8 for the old stock tanks

class Phase(object):
    def __init__(self, name, deltaV, motor):
        self._name = name
        self._deltaV = deltaV
        if isinstance(motor, engine.engine):
            self._Isp = motor.IspVac
            self._isJet = False
        else:
            # it's a jet!
            self._Isp = 10000
            self._isJet = True

    def evaluate(self, payload):
        (prop, tanks) = engine.burnMass(self._deltaV, self._Isp, 
                payload, jetBeta if self._isJet else rocketBeta)
        if self._isJet:
            return (prop, tanks, 0)
        else:
            return (prop, 0, tanks)

    def isJet(self):
        return self._isJet

    def beta(self):
        if self._isJet:
            return jetBeta
        else:
            return rocketBeta

phases = (
        Phase("KSC-LKO", 6000, jets.turbojet),
        Phase("LKO circ", 30, rocket),
        Phase("LKO-Jool", 2100, rocket),
        Phase("Laythe-LLO", 3000, jets.turbojet),
        Phase("LLO-Kerbin", 1100, rocket)
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
    liftoffSpeed = liftSurface.speedForLift(liftSurfaceAoA, takeoffLift, altitude, planet)
    print("%s speed %g m/s" % (name, liftoffSpeed))


nintakes = 0
nlift = 0
nturbo = 1
takeoffSpeed = 60
runwayAltitude = 70
liftSurfaceAoA = 60

fixedPayload = sum( (
            7.960, # cockpits, lab, etc
            2 * rocket.mass,
            6 * lift.avr8.mass, # 2 each for roll/yaw/pitch
            ))

for i in range(1, 10):
    payload = sum( (
                fixedPayload,
                nturbo * jets.turbojet.mass, # turbojets
                nintakes * intake.mass, # intakes
                nlift * liftSurface.mass,
                ))
    stages = list(reversed(evaluateStages(payload, 0)))
    launchMass = stages[0].mass

    needLiftForce = launchMass * planet.kerbin.gravity(runwayAltitude)
    liftPerSCS = liftSurface.liftForce(liftSurfaceAoA, takeoffSpeed, runwayAltitude)
    newNcontrol = math.ceil(needLiftForce / liftPerSCS)

    parts = ((jets.turbojet, nturbo), (intake, nintakes))
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
print ("%d turbojets, %d intakes, %d lift surfaces pointed %g degrees down"
        % (nturbo, nintakes, nlift, liftSurfaceAoA))

jetTanks = sum( [ x.jetTanks for x in stages ] )
rocketTanks = sum( [ x.rocketTanks for x in stages ] )
print ("%gt rocket fuel and %gt jet fuel" % (rocketTanks * rocketBeta, jetTanks * jetBeta))
print ("%g U LiquidFuel, %g U Oxidizer" % (
            (rocketTanks * rocketBeta * 0.9 + jetTanks * jetBeta) * 200,
            (rocketTanks * rocketBeta * 1.1) * 200))
for (phase, stage) in zip(phases, stages):
    name = phase._name
    if phase.isJet():
        fuelName = 'jet'
        tankMass = stage.jetTanks
    else:
        fuelName = 'rocket'
        tankMass = stage.rocketTanks
    startMass = stage.mass
    fuelMass = phase.beta() * tankMass
    fuelUnits = fuelMass * 200
    print ("%10s:  %7.3f t %s fuel (%g units)" %
            (name, fuelMass, fuelName, fuelUnits))


print ""
reportLiftoffSpeed("KSC liftoff", launchMass, planet.kerbin, runwayAltitude)
laytheLandingMass = stages[2].mass
reportLiftoffSpeed("Laythe landing & liftoff", stages[2].mass, planet.laythe, 200)
kerbinLandingMass = fixedPayload + jetTanks + rocketTanks
reportLiftoffSpeed("KSC landing", kerbinLandingMass, planet.kerbin, runwayAltitude)
