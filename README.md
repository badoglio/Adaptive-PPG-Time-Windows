# Adaptive PPG Windows 🌀

> **Adaptive Biosegnal (PPG) Analysis via Variable Temporal Windows and Dynamic Overlap Guided by Morphological Complexity**

This repository contains the complete implementation of a scientific pipeline for extracting and analyzing morphological features from Photoplethysmogram (PPG) signals. Unlike traditional analysis systems that rely on static window intervals, this framework applies a **variable-width ensemble averaging** algorithm that dynamically adapts in real time to the transition velocity of the subject's vascular state.

---

## 📖 Physiological Rationale

In real-world biosegnal analysis, a central bottleneck is the arbitrary trade-off between **signal stationarity** and **temporal resolution**.

### 1. The Limitation of Traditional Fixed Windows
Classic windows (such as the 5-minute interval recommended by the *1996 Task Force* for Heart Rate Variability - HRV) were designed to guarantee a stable spectral estimate in low frequency (LF band at $0.04 \text{ Hz}$). However, the real physiological stress response of the Autonomic Nervous System (ANS) is neither stationary nor synchronous across subsystems:
*   **Vagal (Parasympathetic) Activity:** Reacts in fractions of a second or a few seconds (HF band: $0.15 - 0.40 \text{ Hz}$).
*   **Sympathetic Activity:** Exhibites much slower dynamics, with latencies of $3 - 5 \text{ seconds}$ and responses spanning tens of seconds to minutes (LF band: $0.04 - 0.15 \text{ Hz}$).

Averaging these two pathways inside a single fixed window (such as 30s or 5m) smooths out rapid transients (blunting sympathetic spikes) and mixes uncorrelated slow and fast dynamics.

### 2. From Technical Constraints to Physiology
Instead of forcing human physiology to match a rigid spectral window, this project inverts the paradigm: **the dynamics of the observed feature (the instant vascular state complexity) dictate segment length $T_w$ and overlap ratio $O_w$**.
*   **During stationary phases (low complexity):** The window expands ($T_w \to T_{max}$, e.g. 30-40 beats). The ensemble average spans more cycles, suppressing random noise and maximizing the Signal-to-Noise Ratio ($SNR$).
*   **During transient phases (high complexity, e.g. acute stress or vasoconstriction):** The window contracts rapidly ($T_w \to T_{min}$, e.g. 3-5 beats). The averaging window narrows to trace the instant evolution of the vascular state.

---

## 🧠 Algorithmic Pipeline

To implement window adaptation while avoiding the **circular reasoning trap** (*"calculating window size based on a feature that itself requires a window to be computed"*), the pipeline adopts a **Two-Stage Bounded Adaptive Framework**:

```
[ Raw PPG Signal ]
        │
        ▼ (io_loader & preprocessing)
[ Segmented Normalized Beats [0, 1] ]
        │
        ▼
[ Phase 1: Compute Morphological Entropy H_morph ]  ◄── Ex-Ante Non-Stationarity Index
        │
        ├───────────────────────────────────────┐
        ▼ (Stage 2: Linear Mapping)             ▼ (Spectral Decoupling)
[ Target Window Tw & Overlap Ow ]       [ Critical Sensitivity Bounds T_crit^k ]
        │                                       │
        └───────────────────┬───────────────────┘
                            ▼
              [ Effective Windows T_w^k ]
                            │
                            ▼ (feature_extraction)
             [ Feature Extraction & Benchmark ]
```

### Stage A: Instant Complexity Estimation ($H_{morph}[n]$)
Over the matrix of segmented, amplitude-normalized $[0,1]$ and duration-resampled (128 samples) beats, we compute the **Morphological Entropy ($H_{morph}$)** over the past $N_{past}$ beats (default $10$). The algorithm evaluates consecutive Root Mean Squared (RMS) differences between successive beat profiles:
*   $H_{morph} \approx 0$ (Stationary State): Successive beats are highly matching. Vasomotor tone is stable.
*   $H_{morph} \to 1$ (Transient State): The morphology varies rapidly from beat to beat (indicating vasoconstriction or compliance adjustments).

### Stage B: Sigmoidal/Linear Window ($T_w$) and Overlap ($O_w$) Mapping
The target window size (in beats) and overlap percentages are mapped dynamically:
$$T_w[n] = T_{max} - (T_{max} - T_{min}) \cdot H_{morph}[n]$$
$$O_w[n] = O_{min} + (O_{max} - O_{min}) \cdot H_{morph}[n]$$
During rapid transients ($H_{morph} \to 1$), the window contracts and the overlap increases up to $80-90\%$ to track transitions smoothly.

### Stage C: Spectral Decoupling and Critical Sensitivity ($T_{crit}^k$)
To resolve the time-frequency uncertainty (overly short windows render numerical derivatives unstable), we apply category-specific critical bounds $T_{crit}^k$:
$$T_w^k[n] = \max\left( T_{target}[n], \; T_{crit}^k \right)$$

| Feature Category | Examples | Theoretical $T_{min}$ | Stationary $T_{max}$ | Critical Bound $T_{crit}^k$ | Physiological Meaning |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Macro-Morphological** | *Pulse Slope Index*, *Pulse Width 50%*, *Area Ratio* | $3 \text{ beats}$ | $30 \text{ s}$ | $\ge 3 \text{ beats}$ | Primary arterial vasomotor tone and compliance. |
| **Time-Volume** | *Systolic Peak Time*, *Pulse Transit Duration* | $3 \text{ beats}$ | $30 \text{ s}$ | $\ge 3 \text{ beats}$ | Ventricular ejection velocity and propagation. |
| **Derivatives (VPG / SDPTG)** | $b/a$ or $d/a$ ratios from second derivative | $5 \text{ beats}$ | $60 \text{ s}$ | $\ge 10 \text{ beats}$ | Vascular stiffness (SDPTG requires more cycles for numerical stability). |

---

## 🛠️ Software Architecture & Modules

The codebase is built as a modular pipeline:

1.  📂 **[io_loader.py](file:///f:/__WB/AntiGravity%20Projects/Adaptive%20PPG%20Windows/io_loader.py):** I/O Manager. Loads CSV/TXT (Pandas), MAT (SciPy), and EDF/BDF (pyedflib) files from disk or in-memory streams, returning a unified `SignalData` object.
2.  𝄠 **[preprocessing.py](file:///f:/__WB/AntiGravity%20Projects/Adaptive%20PPG%20Windows/preprocessing.py):** Zero-phase Butterworth bandpass filter ($0.5 - 30 \text{ Hz}$), linear detrending, beat segmentation (BioSPPy), min-max amplitude scaling, and 128-point interpolation.
3.  🌀 **[adaptive_engine.py](file:///f:/__WB/AntiGravity%20Projects/Adaptive%20PPG%20Windows/adaptive_engine.py):** Central engine. Computes $H_{morph}$ trends, maps target windows/overlaps, applies $T_{crit}$ constraints, and generates index-slicing intervals.
4.  📊 **[feature_extraction.py](file:///f:/__WB/AntiGravity%20Projects/Adaptive%20PPG%20Windows/feature_extraction.py):** Extracts 21 morphology metrics (time, area, VPG, SDPTG $a,b,c,d,e$ points, and Aging Index). Compares variable windows to standard 30-second fixed window baselines.
5.  📈 **[visualization.py](file:///f:/__WB/AntiGravity%20Projects/Adaptive%20PPG%20Windows/visualization.py):** Plotly engine. Generates synchronized subplots for matched zooming, and sets up Streamlit download controls (CSV, JSON).
6.  🚀 **[app.py](file:///f:/__WB/AntiGravity%20Projects/Adaptive%20PPG%20Windows/app.py):** Streamlit interface. Coordinates execution stages and manages Cascade Caching (`@st.cache_data`) for real-time interaction.

---

## 🚀 Getting Started

### Prerequisites
The project requires Python 3.10+ and the following packages:
*   `numpy`
*   `pandas`
*   `scipy`
*   `pyedflib`
*   `biosppy`
*   `peakutils`
*   `plotly`
*   `streamlit`

### 1. Clone the repository and configure virtualenv
```bash
git clone https://github.com/your-username/adaptive-ppg-windows.git
cd adaptive-ppg-windows

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

### 2. Install dependencies
```bash
pip install numpy pandas scipy pyedflib biosppy peakutils plotly streamlit
```

### 3. Launch the dashboard
```bash
streamlit run app.py
```

### 4. Run Isolated Module Tests
Each file contains a test block generating synthetic data. Test them in isolation by running:
```bash
python io_loader.py
python preprocessing.py
python adaptive_engine.py
python feature_extraction.py
python visualization.py
```

---

## 📊 Sample Dataset

To help users test the dashboard, a sample clinical recording is included in the `sample_data/` folder:
- **File**: `sample_data/bidmc_01_Signals.csv`
- **Signal Column**: `PLETH` (which represents the PPG photoplethysmogram)
- **Sampling Frequency**: $125 \text{ Hz}$

### Dataset Source & Citation
This sample is extracted from the **BIDMC PPG and Respiration Dataset** hosted publicly on **PhysioNet**:
> *Pimentel, M. A. F., Charlton, P. H., & Clifton, D. A. (2016). BIDMC PPG and Respiration Dataset (version 1.0.0). PhysioNet. https://doi.org/10.13026/C2GG69.*
> 
> *Pimentel, M. A. F., et al. "Towards a Robust Estimation of Respiratory Rate from Pulse Oximeters." IEEE Transactions on Biomedical Engineering, 64(8), pp. 1914-1923, 2016.*

