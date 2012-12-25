import rockets
import planet

# Automated launcher for a 100-ton payload.
kerbinorbit = 100000
KerbinLiftoff = rockets.liftoffBurn("Depart Kerbin", planet.kerbin, orbit = kerbinorbit, 
                    payload = 100)
KerbinArrival = rockets.deepSpaceBurn("De-orbit", 500, payload = 0.5)

rockets.design( ( KerbinLiftoff, KerbinArrival ), minStageDeltaV = 250 , symmetry = 4)
