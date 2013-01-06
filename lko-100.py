import rockets
import planet

# Automated launcher for a 100-ton payload.  Deorbit the launcher (which weighs
# almost nothing).
kerbinorbit = 100000
KerbinLiftoff = rockets.liftoffBurn("Depart Kerbin", planet.kerbin, orbit = kerbinorbit, 
                    payload = 115)
KerbinArrival = rockets.deepSpaceBurn("De-orbit", 500, payload = 0.08)

rockets.design( ( KerbinLiftoff, KerbinArrival ), minStageDeltaV = 250 , symmetry = 4)
