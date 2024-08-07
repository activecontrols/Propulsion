import csv
from utils import max_bolt_shear, max_bolt_tensile, max_internal_shear, torque
from collections import OrderedDict
from pprint import pprint

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
            max_bolt_shear(D, pitch, PROOF_STRENGTH, LENGTH_ENGAGEMENT, D1, D2, SAFETY_FACTOR),
            max_internal_shear(D, pitch, PROOF_STRENGTH, LENGTH_ENGAGEMENT, D2, SAFETY_FACTOR)
        ]

failure_modes = ("bolt tensile", "bolt thread shear", "internal thread shear")

# pprint(bolt_strengths)


load = float(input("bolt load (lbf): "))

print("Allowable Bolts:")
for (name, dia), strengths in bolt_strengths.items():
    if load <= min(strengths):
        print(f"{name} | Max Load: {round(min(strengths))} lbf ({round(load/min(strengths)*100)}%) | Failure mode: {failure_modes[strengths.index(min(strengths))]} | Torque: {round(torque(K, load, dia, SAFETY_FACTOR)/12, 3)} lb-ft")
