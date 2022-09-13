from tkinter import *
import numpy as np
from cmath import sqrt as csqrt, phase, pi
from scipy import special
from scipy.optimize import newton_krylov

# GUI root window
root = Tk()
root.title("renner-solve-qp4")

# Entries
per_label = Label(root, text="per")
per_label.grid(row=0, column=0)
per = Entry(root, width=50, bg="#303030", fg="white", borderwidth=1)
per.grid(row=0, column=1)

ri_label = Label(root, text="ri")
ri_label.grid(row=1, column=0)
ri = Entry(root, width=50, bg="#303030", fg="white", borderwidth=1)
ri.grid(row=1, column=1)

del_qp_label = Label(root, text="del_qp")
del_qp_label.grid(row=2, column=0)
del_qp = Entry(root, width=50, bg="#303030", fg="white", borderwidth=1)
del_qp.grid(row=2, column=1)

PH_qp1_label = Label(root, text="PH_qp1")
PH_qp1_label.grid(row=3, column=0)
PH_qp1 = Entry(root, width=50, bg="#303030", fg="white", borderwidth=1)
PH_qp1.grid(row=3, column=1)

output = Label(root, text="")
output.grid(row=0, column=2)

run_button = Button(root, text="Run", padx=50, pady=10, command=lambda: QPphase(float(per.get()), float(ri.get()), float(del_qp.get()), float(PH_qp1.get())))
run_button.grid(row=4, column=1, columnspan=2)

# e.insert(0, "Enter your name: ") # default text for input field

# function for button  
def QPphase(per, ri, del_qp, PH_qp1):
    #per = period in seconds
    #ri = well radius in meters
    #del_qp = P/Q in MPa(m^3/s)
    #PH_qp1 = phase lag in fraction of cycle
    omega = 2*pi/per          # radial freq 
    func = lambda D1: float(PH_qp1- np.angle(csqrt(omega*1j/D1)* \
        special.kv(1, np.sqrt(omega*1j/D1)*ri)/\
        special.kv(0, np.sqrt(omega*1j/D1)*ri))/2/pi) #phase function for optirmizer

    D_init = 1e-6 #Diffusivity initial guess
    D_sol = newton_krylov(func, D_init, method='lgmres', verbose=1)
    #Solve for diffusivity
    eta = np.sqrt(omega*1j/D_sol)
    K_0_ri = special.kv(0, eta*ri)
    K_1_ri = special.kv(1, eta*ri)
    T = 1/del_qp/(2*pi*ri/10000*abs(eta*K_1_ri/K_0_ri))/1e6
    #del_qp = P/Q MPa(m^3/s)
    # amplitude from q 1*pi*T*ri/10000*

    print(per, " | ", D_sol, " | ",  T, " | ", T/D_sol )
    output["text"]= f"{per} {D_sol} {T} {T/D_sol}"
    return func(0.1)

# exit button
exit_button = Button(root, text="Exit", padx=30, pady=10, command=root.quit)
exit_button.grid(row=5, column=2, columnspan=1)

root.mainloop()