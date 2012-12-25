import rockets
import planet

# Capsule-manned mission to launch a 2-ton package into LKO
kerbinorbit = 100000
KerbinLiftoff = rockets.liftoffBurn("Depart Kerbin", planet.kerbin, orbit = kerbinorbit, payload = 2)
KerbinArrival = rockets.deepSpaceBurn("De-orbit", 500, payload = 5)

rockets.design( ( KerbinLiftoff, KerbinArrival ), minStageDeltaV = 250 )
