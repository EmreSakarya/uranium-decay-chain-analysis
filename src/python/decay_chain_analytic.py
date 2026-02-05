import numpy as np
import math

# --- CONFIGURATION & CONSTANTS ---
try:
    # Use higher precision if available
    DT = np.longdouble
except AttributeError:
    DT = np.float64

EPS_REL = DT(1e-28)
EPS_ABS = DT(0.0)
EPS_ZERO_DIV = DT(1e-99)

NA = DT('6.02214076e23')   # Avogadro constant (mol^-1)
G_PER_KG = DT(1000)        # grams per kg
T_END = DT('4.503e9')      # Age of Earth (years)

# Isotope Data: (Name, Molar Mass [g/mol], Half-life [years])
ISOTOPES = [
    ("U-238",           DT(238.0), DT('4.5e9')),
    ("U-234",           DT(234.0), DT('2.25e5')),
    ("Th-230",          DT(230.0), DT('8.0e4')),
    ("Ra-226",          DT(226.0), DT('1.6e3')),
    ("Pb-210",          DT(210.0), DT('22.0')),
    ("Pb-206 (stable)", DT(206.0), np.inf),
]

NAMES = [n for n, _, _ in ISOTOPES]
A_MASS = np.array([a for _, a, _ in ISOTOPES], dtype=DT)
T12 = np.array([t for _, _, t in ISOTOPES], dtype=object)

# Decay Constants (lambda = ln(2) / t_1/2)
LAM = np.array([
    (DT(math.log(2.0)) / DT(t)) if np.isfinite(t) else DT(0) 
    for t in T12
], dtype=DT)

NUM_ISOTOPES = len(LAM)

# --- HELPER FUNCTIONS ---

def mass_to_atoms(mass_kg, molar_mass):
    """Converts mass (kg) to number of atoms."""
    moles = (mass_kg * G_PER_KG) / molar_mass
    return moles * NA

def atoms_to_mass(N, molar_mass):
    """Converts number of atoms to mass (kg)."""
    moles = N / NA
    return (moles * molar_mass) / G_PER_KG

# --- SOLVERS ---

def bateman_solution(t, lam, N0_parent):
    """
    Computes analytical solution using Bateman equations.
    Assumes only the first parent (N0) is present at t=0.
    """
    N = np.zeros(NUM_ISOTOPES, dtype=DT)
    
    if t == 0:
        N[0] = N0_parent
        return N

    # Parent Nucleus
    N[0] = N0_parent * np.exp(-lam[0] * t)

    # Daughter Nuclei
    for k in range(1, NUM_ISOTOPES - 1):
        num = DT(1)
        for r in range(k):
            num *= lam[r]
        
        total_sum = DT(0)
        for j in range(k + 1):
            denom = DT(1)
            for r in range(k + 1):
                if r != j:
                    diff = lam[r] - lam[j]
                    # Avoid division by zero for very close decay constants
                    denom *= diff if abs(diff) > DT(1e-100) else DT(1e-100)
            
            total_sum += np.exp(-lam[j] * t) / denom
        
        N[k] = N0_parent * num * total_sum

    # Stable End Product (Conservation of atoms)
    N[-1] = N0_parent - np.sum(N[:-1], dtype=DT)
    
    return N

def forward_euler_solver(dt, t_total, lam, N0_parent):
    """
    Numerical solution using Explicit (Forward) Euler method.
    """
    steps = int(np.ceil(t_total / dt))
    N = np.zeros(NUM_ISOTOPES, dtype=DT)
    N[0] = N0_parent
    
    for _ in range(steps):
        prev = N.copy()
        
        # Parent
        N[0] = prev[0] + (-lam[0] * prev[0]) * dt
        if N[0] < 0: N[0] = DT(0)
        
        # Daughters
        for k in range(1, NUM_ISOTOPES - 1):
            production = lam[k-1] * prev[k-1]
            decay = lam[k] * prev[k]
            N[k] = prev[k] + (production - decay) * dt
            if N[k] < 0: N[k] = DT(0)
            
        # Stable Product
        N[-1] = prev[-1] + (lam[NUM_ISOTOPES-2] * prev[NUM_ISOTOPES-2]) * dt
        
    return N

# --- MAIN EXECUTION ---

def run_analysis():
    print("=== PROJECT 1: URANIUM SERIES DECAY ANALYSIS ===\n")

    # --- PART 1: MASS EVOLUTION (Analytic) ---
    M_NOW_U238_KG = DT('2.892e15')
    N_now_U238 = mass_to_atoms(M_NOW_U238_KG, A_MASS[0])
    
    # Back-calculate initial atoms at t=0
    N0_U238 = N_now_U238 * np.exp(LAM[0] * T_END)
    
    N_start = bateman_solution(DT(0), LAM, N0_U238)
    N_now = bateman_solution(T_END, LAM, N0_U238)
    
    m_start = [atoms_to_mass(N_start[i], A_MASS[i]) for i in range(NUM_ISOTOPES)]
    m_now = [atoms_to_mass(N_now[i], A_MASS[i]) for i in range(NUM_ISOTOPES)]
    
    print("--- Q1: Analytic Masses (kg) ---")
    print(f"{'Isotope':<18} {'At Formation (t=0)':<25} {'Today (t=4.5e9)':<25}")
    for i in range(NUM_ISOTOPES):
        print(f"{NAMES[i]:<18} {float(m_start[i]):.6e}                        {float(m_now[i]):.6e}")
    print("\n")

    # --- PART 2: NUMERICAL vs ANALYTIC (Euler Method) ---
    print("--- Q2: Numerical Stability Analysis (Forward Euler) ---")
    dt_list = [DT('1000'), DT('10000'), DT('1000000')]
    
    ref_mass = np.sum(m_now) # Conservation of mass check
    
    for dt in dt_list:
        print(f"\n[Time Step dt = {int(dt)} years]")
        
        # Stability Check
        unstable_isotopes = []
        for i in range(NUM_ISOTOPES - 1):
            if float(LAM[i] * dt) >= 1.0:
                unstable_isotopes.append(NAMES[i])
        if unstable_isotopes:
            print(f"Warning: Solution likely unstable for: {', '.join(unstable_isotopes)}")
            
        N_num = forward_euler_solver(dt, T_END, LAM, N0_U238)
        m_num = [atoms_to_mass(N_num[i], A_MASS[i]) for i in range(NUM_ISOTOPES)]
        
        print(f"{'Isotope':<18} {'Analytic (kg)':<18} {'Numeric (kg)':<18} {'Rel Err (%)':<18}")
        for i in range(NUM_ISOTOPES):
            a = float(m_now[i])
            n = float(m_num[i])
            diff = abs(n - a) / a * 100 if a > 0 else 0.0
            print(f"{NAMES[i]:<18} {a:.4e}            {n:.4e}            {diff:.4f}%")

    # --- PART 3: BACK-CALCULATION FROM PB-206 ---
    print("\n--- Q3: Reverse Calculation (Pb-206 -> U-238) ---")
    
    M_PB206_TODAY = DT('6.44e18')
    N_Pb206_today = mass_to_atoms(M_PB206_TODAY, A_MASS[-1])
    
    # Calculate what fraction of U-238 has decayed into Pb-206 by today
    # We use Bateman to find a[-1] which is the coefficient for the last daughter
    N_normalized = bateman_solution(T_END, LAM, DT(1.0))
    fraction_converted = N_normalized[-1]
    
    # Calculate Initial U-238 required to produce this much Pb-206
    N0_required = N_Pb206_today / fraction_converted
    M_U238_start_calc = atoms_to_mass(N0_required, A_MASS[0])
    
    # Calculate remaining U-238 today from that start
    M_U238_today_calc = M_U238_start_calc * np.exp(-float(LAM[0]) * float(T_END))
    
    print(f"Given Pb-206 mass:   {float(M_PB206_TODAY):.6e} kg")
    print(f"Calculated U-238(0): {float(M_U238_start_calc):.6e} kg")
    print(f"Calculated U-238(T): {float(M_U238_today_calc):.6e} kg")

if __name__ == "__main__":
    run_analysis()
