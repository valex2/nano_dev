import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
from scipy.stats import norm

# Pre-defined fluids and their properties (density in kg/m³, viscosity in Pa·s)
fluids = {
    "DMF": {"density": 944, "viscosity": 0.00092},
    "3M 7200": {"density": 1400, "viscosity": 0.00061},
    "Hexanes": {"density": 655, "viscosity": 0.000297},
    "Acetone": {"density": 784, "viscosity": 0.000306},
    "50:50 Hexanes:Acetone": {"density": 720, "viscosity": 0.000301},
    "Ethanol": {"density": 789, "viscosity": 0.001074},
    "Water": {"density": 1000, "viscosity": 0.00089},
}

# Particle density (example: iron oxide nanoparticles ~5000 kg/m³)
particle_density = 5000

def sedimentation_time(d, rho_p, rho_s, RCF, eta, h):
    g = 9.81
    v = (d**2 * (rho_p - rho_s) * RCF * g) / (18 * eta)
    t = h / v
    return t

def particle_diameter(t, rho_p, rho_s, RCF, eta, h):
    g = 9.81
    d = np.sqrt((18 * eta * h) / ((rho_p - rho_s) * RCF * g * t))
    return d

# Compute required RCF given sedimentation time and particle size percentile
def required_RCF(d, rho_p, rho_s, t, eta, h):
    g = 9.81
    RCF = (18 * eta * h) / (d**2 * (rho_p - rho_s) * g * t)
    return RCF

# Streamlit GUI
st.title("Nanoparticle Centrifugation Simulator")

fluid_choice = st.selectbox("Select Fluid", list(fluids.keys()))
fluid_props = fluids[fluid_choice]

col1, col2, col3 = st.columns(3)
with col1:
    RCF = st.number_input("RCF (g)", value=10000, step=500)
with col2:
    t_cent_min = st.number_input("Centrifugation Time (min)", value=10, step=1)
    t_cent = t_cent_min * 60  # convert to seconds

with col3:
    tube_type = st.selectbox("Tube Type", ["50ml Eppendorf", "15ml Eppendorf", "2ml Eppendorf"])
    if tube_type == "50ml Eppendorf":
        fill_ml = st.selectbox("Fill Level (ml)", [5, 10, 15, 20, 25, 30, 40, 50])
        h_cm = fill_ml / 5
    elif tube_type == "15ml Eppendorf":
        fill_ml = st.selectbox("Fill Level (ml)", [0.8, 1.6, 2.4, 3.2, 4.0, 8.0, 12.0, 15.0])
        h_cm = fill_ml / 0.8
    elif tube_type == "2ml Eppendorf":
        fill_ml = st.selectbox("Fill Level (ml)", [0.1, 0.5, 1.0, 1.5, 2.0])
        h_cm = fill_ml * 2
    h = h_cm / 100  # convert cm to m

st.markdown("---")
st.subheader("Initial Particle Size Distribution")
mean_diameter_nm = st.number_input("Mean Diameter (nm)", value=50.0)
sd_diameter_nm = st.number_input("Std Deviation (nm)", value=10.0)

# Convert particle diameter from nm to m for calculations
mean_diameter_m = mean_diameter_nm * 1e-9
sd_diameter_m = sd_diameter_nm * 1e-9

# Simulate particle distribution
diameters_nm = np.linspace(mean_diameter_nm - 4*sd_diameter_nm, mean_diameter_nm + 4*sd_diameter_nm, 500)
diameters_m = diameters_nm * 1e-9
distribution = norm.pdf(diameters_nm, mean_diameter_nm, sd_diameter_nm)

# Compute critical diameter based on centrifugation conditions
critical_diameter_m = particle_diameter(t_cent, particle_density, fluid_props["density"], RCF, fluid_props["viscosity"], h)
critical_diameter_nm = critical_diameter_m * 1e9

pellet_fraction = diameters_nm >= critical_diameter_nm
supernatant_fraction = diameters_nm < critical_diameter_nm

# Plotting
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(diameters_nm, distribution, label="Particle Size Distribution", color='blue')
ax.fill_between(diameters_nm, distribution, where=pellet_fraction, color='red', alpha=0.5, label='Pellet')
supernatant = ~pellet_fraction
ax.fill_between(diameters_nm, distribution, where=supernatant, color='green', alpha=0.5, label='Supernatant')

# Critical diameter line
ax.axvline(critical_diameter_nm, color='black', linestyle='--', label=f"Critical Diameter ({critical_diameter_nm:.2f} nm)")

ax.set_xlabel("Particle Diameter (nm)")
ax.set_ylabel("Probability Density")
ax.set_title("Nanoparticle Centrifugation Simulation")
ax.legend()

percentile = st.slider("Desired Percentile to Precipitate (%)", min_value=0, max_value=100, value=95)
required_diameter_m = norm.ppf((100 - percentile) / 100, mean_diameter_m, sd_diameter_m)
required_diameter_nm = required_diameter_m * 1e9

ax.axvline(required_diameter_nm, color='orange', linestyle='--', label=f"Mean Diameter ({mean_diameter_nm:.2f} nm)")

st.pyplot(fig)

# Additional functionality
st.markdown("---")
st.subheader("Compute RCF or Time for Desired Sedimentation")

col4, col5 = st.columns(2)

with col4:
    given_time_min = st.number_input("Given Time (min)", value=10, step=1, key='time')
    given_time_sec = given_time_min * 60
    req_RCF = required_RCF(required_diameter_m, particle_density, fluid_props["density"], given_time_sec, fluid_props["viscosity"], h)
    st.write(f"Required RCF: {req_RCF:.0f} g")

with col5:
    given_RCF = st.number_input("Given RCF (g)", value=10000, step=500, key='rcf')
    req_time_sec = sedimentation_time(required_diameter_m, particle_density, fluid_props["density"], given_RCF, fluid_props["viscosity"], h)
    req_time_min = req_time_sec / 60
    st.write(f"Required Time: {req_time_min:.2f} min")
