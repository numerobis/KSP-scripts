import rockets
import planet

# Automated launcher for a 100-ton payload.

eveOrbit = 100000
EveUp = rockets.liftoffBurn("Depart Eve", planet.eve, orbit = eveOrbit)

v_orbit = planet.eve.orbitalVelocity(eveOrbit)
v_eject = 4572.47

# Tiny return stage; give 100m/s for corrections.
KerbinReturn = rockets.deepSpaceBurn("Eve->Kerbin", v_eject - v_orbit + 100, payload = 0.3)

# Build like a mofo.  This will take forever: Eve requires about 11.7km/s to
# orbit, which means we're splitting into more than 40 stages.  2-way symmetry
# really means asparagus spirals.
rockets.design( ( EveUp, KerbinReturn ), minStageDeltaV = 250 , symmetry = 2)
