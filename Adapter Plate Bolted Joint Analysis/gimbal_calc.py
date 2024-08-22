import csv
from utils import max_bolt_thread_shear, max_bolt_tensile, max_internal_thread_shear, max_bolt_shear, torque
from collections import OrderedDict
from pprint import pprint
import math

"""
Generates a table of max loads for a range of bolts (given the constants), then checks the input load against that
table, then calculates torques for the usable bolts.

Cross-checked against this proof load chart: https://www.almabolt.com/pages/catalog/bolts/proofloadtensile.htm
"""


# Constants:

# Torque Coefficient
K = 0.17  # (Stainless 404 w/ threadlocker)

# Stainless 404 / grade 8 steel
PROOF_STRENGTH = 130000 * .915  # Can't find online, so approximated.

LENGTH_ENGAGEMENT = 0.375  # nut
SAFETY_FACTOR = 5  


bolt_strengths = OrderedDict()

with open("ASME_B1_1_1989_Bolt_Specs.csv") as bolt_specs:
    reader = csv.reader(bolt_specs, delimiter=",")
    next(reader)  # get rid of the column names
    for row in reader:
        name = row[0]
        D = float(row[1])
        pitch = float(row[2])
        D2 = float(row[3])
        D1 = float(row[5])
        bolt_strengths[(name, D)] = [
            max_bolt_tensile(D, pitch, PROOF_STRENGTH, SAFETY_FACTOR),
            max_bolt_thread_shear(D, pitch, PROOF_STRENGTH, LENGTH_ENGAGEMENT, D1, D2, SAFETY_FACTOR),
            max_internal_thread_shear(D, pitch, PROOF_STRENGTH, LENGTH_ENGAGEMENT, D2, SAFETY_FACTOR),
            max_bolt_shear(D, pitch, PROOF_STRENGTH, SAFETY_FACTOR)
        ]

failure_modes = ("bolt tensile", "bolt thread shear", "internal thread shear", "bolt shear")

# pprint(bolt_strengths)


load = float(input("bolt load (lbf): "))
load_angle = float(input("load angle (measured from vertical) (deg): "))
num_bolts = int(input("# of bolts:"))
dia = float(input("diameter of bolt circle (in):"))
z_dist = float(input("z-distance between load and (mating) surface:"))

tensile_load = math.sin((90 - load_angle)*2*math.pi/360) * load
shear_load = math.cos((90 - load_angle)*2*math.pi/360) * load

# Following load calculations based on AISC Steel Construction Manual 7-12
# Assumes that the load is located along the line concentric to the circle center (more conservative anyways)
# The half of the bolts opposite the direction of the load have a "prying" tensile force applied

shear_load = shear_load / num_bolts

r = dia/2
dM = 2*((4*math.pi*r)/(3*math.pi))

shear_induced_tensile_load = (shear_load*z_dist)/((num_bolts/2)*dM)
print("shear-induced tensile load:", shear_induced_tensile_load)
tensile_load += shear_induced_tensile_load  # TODO: if the *into* is into the plate, is does the bolt still under tensile load?



print("\nAllowable Bolts:")
print(f"Tensile load: {tensile_load} | Shear load: {shear_load}")
for (name, dia), strengths in bolt_strengths.items():
    if tensile_load <= min(strengths[:-1]) and shear_load <= strengths[-1]:
        print(f"{name} | Max Load (tensile, shear): {round(min(strengths[:-1]))}, {round(strengths[-1])} lbf ({round(tensile_load/min(strengths[:-1])*100)}%, {round(shear_load/strengths[-1]*100)}%) | Failure mode: {failure_modes[strengths.index(min(strengths))]} | Torque: {round(torque(K, load, dia, SAFETY_FACTOR)/12, 3)} lb-ft")
