import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# ─────────────────────────────────────────────
# EDIT THESE TO MATCH YOUR OUTPUT FILE NAMES
PMF_FILE  = "/media/supremeleader/Pantera/simulation/analysis_2024/average_water_molecule/pmf.xvg"       # output of -o flag in gmx wham
HIST_FILE = "/media/supremeleader/Pantera/simulation/analysis_2024/average_water_molecule/histo.xvg"     # output of -hist flag in gmx wham

# Plot labels — change to your system name
SYSTEM_TITLE = "Potential of Mean Force"
SYSTEM_SUBTITLE = "TRIC Channel"          # e.g. "6IYX Eukaryotic TRIC"
X_LABEL = "ξ (nm)"
# ─────────────────────────────────────────────


def read_xvg(filename):
    """
    Read a GROMACS .xvg file, skipping comment/header lines (@ and #).
    Returns a 2D numpy array where each column is a data series.
    """
    data = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('#') or line.startswith('@'):
                continue
            if line:
                data.append([float(x) for x in line.split()])
    return np.array(data)


def plot_pmf_and_histograms(pmf_file, hist_file):

    # ── Load PMF data ──────────────────────────────────────────────
    print(f"Reading PMF file: {pmf_file}")
    pmf_data = read_xvg(pmf_file)
    pmf_x = pmf_data[:, 0]   # reaction coordinate (nm)
    pmf_y = pmf_data[:, 1]   # PMF values (kcal/mol)

    # ── Load Histogram data ────────────────────────────────────────
    print(f"Reading histogram file: {hist_file}")
    hist_data = read_xvg(hist_file)
    hist_x = hist_data[:, 0]           # position (nm)
    hist_cols = hist_data[:, 1:]       # one column per umbrella window

    n_windows = hist_cols.shape[1]
    print(f"  Found {n_windows} umbrella windows in histogram file")

    # ── Color cycle (matches tutorial style) ──────────────────────
    colors = ['black', 'red', 'green', 'blue', 'magenta',
              'cyan', 'orange', 'purple', 'brown', 'pink',
              'gray', 'olive', 'navy', 'gold', 'teal']

    # ── Create figure with two panels ─────────────────────────────
    fig = plt.figure(figsize=(8, 8))
    gs  = gridspec.GridSpec(2, 1, height_ratios=[1, 1], hspace=0.08)

    ax1 = fig.add_subplot(gs[0])   # PMF panel
    ax2 = fig.add_subplot(gs[1])   # Histogram panel

    # ── Top panel: PMF ────────────────────────────────────────────
    ax1.plot(pmf_x, pmf_y, color='black', linewidth=1.5)
    ax1.set_ylabel("PMF (kcal mol$^{-1}$)", fontsize=12)
    ax1.set_xlim(pmf_x.min(), pmf_x.max())
    ax1.set_title(f"{SYSTEM_TITLE}\n{SYSTEM_SUBTITLE}", fontsize=13)
    ax1.tick_params(labelbottom=False)   # hide x tick labels on top panel
    ax1.axhline(0, color='gray', linewidth=0.5, linestyle='--')
    ax1.grid(False)

    # ── Bottom panel: Histograms ───────────────────────────────────
    for i in range(n_windows):
        color = colors[i % len(colors)]
        ax2.plot(hist_x, hist_cols[:, i], color=color, linewidth=1.0)

    ax2.set_xlabel(X_LABEL, fontsize=12)
    ax2.set_ylabel("count", fontsize=12)
    ax2.set_xlim(pmf_x.min(), pmf_x.max())
    ax2.grid(False)

    # ── Align x-axes of both panels ───────────────────────────────
    ax1.set_xlim(ax2.get_xlim())

    # ── Save and show ─────────────────────────────────────────────
    output_name = "pmf_histogram_plot.png"
    plt.savefig(output_name, dpi=300, bbox_inches='tight')
    print(f"\n✔ Plot saved as: {output_name}")
    plt.show()


# ── Run ───────────────────────────────────────────────────────────
if __name__ == "__main__":
    try:
        plot_pmf_and_histograms(PMF_FILE, HIST_FILE)
    except FileNotFoundError as e:
        print(f"\nERROR: {e}")
        print("Make sure pmf.xvg and histo.xvg are in the same directory as this script.")
        print("Or edit PMF_FILE and HIST_FILE at the top of the script.")


