import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

st.set_page_config(page_title="Acid-Base Titration Simulator", layout="centered")
st.title("🔬 Acid-Base Titration Simulator")

st.sidebar.header("Titration Setup")
acid_type = st.sidebar.selectbox("Select Acid Type", ["Strong Acid (HCl)", "Weak Acid (CH3COOH)"])
base_type = st.sidebar.selectbox("Select Base Type", ["Strong Base (NaOH)", "Weak Base (NH3)"])

acid_conc = st.sidebar.number_input("Acid Concentration (mol/L)", min_value=0.01, value=0.1, step=0.01)
acid_vol = st.sidebar.number_input("Acid Volume (mL)", min_value=1.0, value=25.0)
base_conc = st.sidebar.number_input("Base Concentration (mol/L)", min_value=0.01, value=0.1, step=0.01)
indicator = st.sidebar.selectbox("Choose Indicator", ["Phenolphthalein", "Methyl Orange", "Bromothymol Blue"])

step_size = st.sidebar.slider("Base Addition Step (mL)", min_value=0.1, max_value=5.0, value=0.5, step=0.1)

indicator_range = {
    "Phenolphthalein": (8.2, 10.0),
    "Methyl Orange": (3.1, 4.4),
    "Bromothymol Blue": (6.0, 7.6)
}

indicator_color_map = {
    "Phenolphthalein": {"low": "colorless", "high": "pink"},
    "Methyl Orange": {"low": "red", "high": "yellow"},
    "Bromothymol Blue": {"low": "yellow", "high": "blue"}
}

pKa_acetic = 4.76
Kb_ammonia = 1.8e-5

V_base = np.arange(0, 50 + step_size, step_size)
pH = []

for V in V_base:
    V_total = acid_vol + V
    if acid_type.startswith("Strong") and base_type.startswith("Strong"):
        moles_H = acid_conc * acid_vol / 1000
        moles_OH = base_conc * V / 1000
        excess = moles_H - moles_OH
        if excess > 0:
            ph = -np.log10(excess / V_total * 1000)
        elif excess < 0:
            ph = 14 + np.log10(abs(excess) / V_total * 1000)
        else:
            ph = 7
    elif acid_type.startswith("Weak") and base_type.startswith("Strong"):
        Ka = 10**(-pKa_acetic)
        moles_HA = acid_conc * acid_vol / 1000
        moles_OH = base_conc * V / 1000
        if moles_OH < moles_HA:
            moles_A = moles_OH
            moles_HA -= moles_OH
            ph = pKa_acetic + np.log10(moles_A / moles_HA)
        elif moles_OH == moles_HA:
            ph = pKa_acetic
        else:
            excess_OH = (moles_OH - moles_HA) / V_total * 1000
            ph = 14 + np.log10(excess_OH)
    elif acid_type.startswith("Strong") and base_type.startswith("Weak"):
        moles_H = acid_conc * acid_vol / 1000
        moles_B = base_conc * V / 1000
        if moles_B == 0:
            ph = -np.log10(acid_conc)
        else:
            OH_conc = np.sqrt(Kb_ammonia * (moles_B / V_total * 1000))
            ph = 14 + np.log10(OH_conc)
    else:
        ph = 7

    pH.append(ph)

# Determine indicator color change
low, high = indicator_range[indicator]
colors = []
for val in pH:
    if val < low:
        colors.append(indicator_color_map[indicator]["low"])
    elif val > high:
        colors.append(indicator_color_map[indicator]["high"])
    else:
        colors.append("orange")

st.subheader("📈 Drop-by-Drop Titration Control")
st.write("Use the slider below to add base incrementally.")
drop_index = st.slider("Volume of Base Added (mL)", min_value=0.0, max_value=float(V_base[-1]), value=0.0, step=step_size)

fig, ax = plt.subplots()
ax.set_xlim(0, 50)
ax.set_ylim(0, 14)
ax.set_xlabel("Volume of Base Added (mL)")
ax.set_ylabel("pH")
ax.set_title(f"Titration Curve ({indicator})")
ax.grid(True)

plot_index = int(drop_index / step_size) + 1
ax.scatter(V_base[:plot_index], pH[:plot_index], c=colors[:plot_index], s=10)
ax.axhline(low, color='gray', linestyle='--', linewidth=0.8)
ax.axhline(high, color='gray', linestyle='--', linewidth=0.8)
st.pyplot(fig)

# Conductivity approximation
st.subheader("🔌 Conductometric Titration")
conductivity = []
for i, V in enumerate(V_base):
    if acid_type.startswith("Strong") and base_type.startswith("Strong"):
        if V < acid_vol:
            cond = 10 - 0.1 * V
        elif V == acid_vol:
            cond = 5
        else:
            cond = 5 + 0.15 * (V - acid_vol)
    else:
        cond = 5 + 0.02 * (V - acid_vol)
    conductivity.append(cond)

fig2, ax2 = plt.subplots()
ax2.plot(V_base, conductivity, color='purple')
ax2.set_xlabel("Volume of Base Added (mL)")
ax2.set_ylabel("Conductivity (a.u.)")
ax2.set_title("Conductometric Titration Curve")
ax2.grid(True)
st.pyplot(fig2)

# Export results
if st.button("📤 Export Results to CSV"):
    df_export = pd.DataFrame({
        "Volume_Base_mL": V_base,
        "pH": pH,
        "Conductivity": conductivity
    })
    df_export.to_csv("titration_results.csv", index=False)
    st.success("✅ Results exported to titration_results.csv")

st.markdown("""
### Notes:
- Drop-by-drop pH visualization is now supported.
- Choose your step size and control titration interactively.
- Indicator color reflects pH values in real-time.
""")
