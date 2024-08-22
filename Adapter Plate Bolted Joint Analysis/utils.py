import math

"""
Derived from:
ASME B1.1-1989 (https://www.iranqc.com/wp-content/uploads/2019/06/ASME-B1.1-1989-EBVvQBF.pdf)
https://purdue-space-program.atlassian.net/wiki/spaces/PAC/pages/710771240/Fastener+Standards+and+Torque+Sizing

----------

Used Nomenclature/Notes: 
- D  Basic (nominal) bolt diameter

- d = External thread Major diameter
- d2 = External thread Pitch diameter

- D1 = Internal thread Minor diameter 
- D2 = Internal thread Pitch diameter

- K = Torque coefficient (given)

- D2(max) = d2(min)
- D = d   
"""


def max_bolt_tensile(D, pitch, proof_strength, safety_factor):
    area = 0.7854 * (D - (0.9743 / pitch)) ** 2
    return (proof_strength * area)/safety_factor


def max_bolt_thread_shear(D, pitch, proof_strength, length_engagement, D1, d2, safety_factor):
    area = ((length_engagement * D1 * math.pi)/pitch) * ((pitch/2) + math.tan(math.pi/6)*(d2 - D1))
    return (0.5 * proof_strength * area)/safety_factor


def max_internal_thread_shear(d, pitch, yield_strength, length_engagement, D2, safety_factor):
    area = ((length_engagement * d * math.pi) / pitch) * ((pitch/2) + math.tan(math.pi/6) * (d - D2))
    return (0.5 * yield_strength * area) / safety_factor

def max_bolt_shear(D, pitch, proof_strength, safety_factor):
    return 0.6 * max_bolt_tensile(D, pitch, proof_strength, safety_factor)  # common approximation !? 


def torque(K, load, D, safety_factor):
    return K*load*D*safety_factor