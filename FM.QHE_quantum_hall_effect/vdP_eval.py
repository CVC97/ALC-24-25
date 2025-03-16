import numpy as np
import matplotlib.pyplot as plt
import os

# needed data
d = 100 # sample thickness
I_hall = 0.002 # in A
B = np.array([200,400,600])
e = 1.6021766*10**(-19)

# measured values
Vr1234 = np.array([-0.322060, -0.322060, -0.322060])
Vr2134 = np.array([0.333062, 0.333062, 0.333062])
Vr2341 = np.array([-0.244881, -0.244881, -0.244881])
Vr3241 = np.array([0.237033, 0.237033, 0.237033])
Vr3412 = np.array([-0.330457, -0.330457, -0.330457])
Vr4312 = np.array([0.313789, 0.313789, 0.313789])
Vr4123 = np.array([-0.227492, -0.227492, -0.227492])
Vr1423 = np.array([0.247476, 0.247476, 0.247476])

vh1342_plus = np.array([0.088215, 0.098297, 0.108439])
vh3142_plus = np.array([-0.100658, -0.110628, -0.120678])
vh2413_plus = np.array([-0.073427, -0.063018, -0.052603])
vh4213_plus = np.array([0.081455, 0.071872, 0.062328])

vh1342_minus = np.array([0.068840, 0.058322, 0.048213])
vh3142_minus = np.array([-0.081760, -0.071103, -0.060891])
vh2413_minus = np.array([-0.093200, -0.104175, -0.114615])
vh4213_minus = np.array([0.099623, 0.109538, 0.119487])
# Resistance values (Hall voltages divided by current)
R3142_plus = vh3142_plus / I_hall
R1342_plus = vh1342_plus / I_hall
R2413_plus = vh2413_plus / I_hall
R4213_plus = vh4213_plus / I_hall

R3142_minus = vh3142_minus / I_hall
R1342_minus = vh1342_minus / I_hall
R2413_minus = vh2413_minus / I_hall
R4213_minus = vh4213_minus / I_hall

# Compute RH
R_H = (1 / 8) * (d / B) *(
    (R3142_plus - R1342_plus) + (R1342_minus - R3142_minus) +
    (R4213_plus - R2413_plus) + (R2413_minus - R4213_minus)
)

# Print the computed RH values
print("Hall Coefficient (RH):", R_H)

##########################################################
# factors
Q_A= (Vr2134-Vr1234) / (Vr3241-Vr2341)
Q_B = (Vr4312-Vr3412) / (Vr1423-Vr4123)
A_fac= ((Q_A-1)/(Q_A+1))**2
B_fac = ((Q_B-1)/(Q_B+1))**2
f_A = 1 - 0.34657*A_fac - 0.09236*A_fac**2 
f_B = 1 - 0.34657*B_fac - 0.09236*B_fac**2


rho = np.pi * d / (np.log(2)*8) * ( f_A *(Vr2134-Vr1234 + Vr3241- Vr2341) + f_B* (Vr4312 - Vr3412 + Vr1423 - Vr4123) )

print(rho)
mobility = R_H/rho
carrier_concentration_per_thickness = 1/(e*R_H/d)
print("mobilities", mobility)
print("carrier concentration per thickness", carrier_concentration_per_thickness)

# Hall resistivity data

# NOW: HALL MEASUREMENTS

folder_path = "/home/caroline/Documents/Physics_MSc/WiSe2024_2025/lab_course/FM_QHE/FMQHE-Daten/250312"  
for file in os.listdir(folder_path):
    if file.endswith(".txt") and file.startswith("Caro_Carlo_vdPauw"):  # Process only vdP file
        file_path = os.path.join(folder_path, file)

        temperature = None # temperature in K
        B = [] # magnetic field
        V_H = [[] for i in range (0,8)] # hall voltage
        #print(len(V_H))
        B_index = 2
        V_H_indices = range(12,20)
        
        # Read file line by line
        with open(file_path, "r") as f:
            for line in f:
                if "#" in line:
                    continue
                elif not line.startswith("#") and line.strip():  
                    values = line.split()
                    if len(values) == 20:
                        #print(values)
                        for i in V_H_indices:
                            #print(values[i])
                            #print(i)
                            V_H[i-V_H_indices[0]].append(float(values[i]))
                        
                        B.append(values[B_index])

V_H = np.array(V_H)
V_H_labels = np.array(["$V_{1342}(B+)$", 	"$V_{3142}(B+)$",	"$V_{2413}(B+)$", 	"$V_{4213}(B+)$", 	"$V_{1342}(B-)$", 	"$V_{3142}(B-)$", 	"$V_{2413}(B-)$", 	"$V_{4213}(B-)$"])
B = np.array(B)
#  check
if V_H.shape[1] != len(B): # dimension check
    raise ValueError(f"Mismatch: B has {len(B)} elements, but V_H has {V_H.shape[1]} columns.")

#plotting
fig, ax = plt.subplots(figsize = (6,3), layout="constrained")
for i in range(0, len(V_H)):
        ax.plot(B, V_H[i], marker='x', linestyle='-', label=f"{V_H_labels[i]}")  
ax.set_xlabel("$B$ in mT")
ax.set_ylabel("$V_H$ in V")
ax.set_yscale("linear")

ax.minorticks_on()
ax.grid(True, which = "both")
fig.legend(loc='outside right upper')
plt.savefig("V_H_over_B.pdf")
plt.show()