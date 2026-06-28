# DMD limitations recreation using the exact equation you shared

This project recreates the Dynamic Mode Decomposition example using this equation:

\[
f(x,t) = f_1(x,t) + f_2(x,t)
\]

where

\[
f_1(x,t) = \mathrm{sech}(x+6-t)e^{i2.3t}
\]

\[
f_2(x,t) = 2\mathrm{sech}(x)\tanh(x)e^{i2.8t}
\]

The two angular frequencies are 2.3 and 2.8.

The code does **not** contain pre-generated results. The `outputs/` folder is created only after you run the script.

---

## How to run in VS Code

Open this folder in VS Code. Then open the terminal inside VS Code.

### 1. Create a virtual environment

Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```

Windows Command Prompt:

```cmd
python -m venv .venv
.venv\Scripts\activate
```

Mac/Linux:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

### 2. Install required packages

```bash
pip install -r requirements.txt
```

### 3. Run the project

```bash
python run_all.py
```

---

## What files will be generated?

After running, the code creates an `outputs/` folder.

Important files:

```text
outputs/rank_error_summary.csv
outputs/svd_energy_table.csv
outputs/01_f_total_real.png
outputs/02_f_total_magnitude.png
outputs/03_svd_singular_values.png
outputs/04_svd_cumulative_energy.png
outputs/rank_2/reconstruction_line_comparison.png
outputs/rank_2/dmd_eigenvalues.png
outputs/rank_2/dmd_mode_magnitudes.png
outputs/rank_2/dmd_eigenvalue_frequency_table.csv
outputs/rank_5/...
outputs/rank_10/...
```

---

## What the code is doing

### 1. It creates the data matrix

The code evaluates `f(x,t)` on a grid.

Rows = different `x` positions.  
Columns = different time snapshots.

So the data matrix looks like:

```text
X = [snapshot at t1, snapshot at t2, snapshot at t3, ...]
```

Each column is the full spatial profile at one time.

### 2. It makes X1 and X2

DMD assumes the system moves from one snapshot to the next:

```text
x2 ≈ A x1
x3 ≈ A x2
x4 ≈ A x3
```

So the code forms:

```text
X1 = [snapshot 1, snapshot 2, ..., snapshot m-1]
X2 = [snapshot 2, snapshot 3, ..., snapshot m]
```

Then DMD tries to find a linear operator `A` such that:

```text
X2 ≈ A X1
```

### 3. It uses reduced-order DMD

The full operator `A` would be very large. Instead, the code uses SVD and projects the data into a smaller rank-r space.

The reduced DMD operator is:

```text
A_tilde = U_r^* X2 V_r Sigma_r^{-1}
```

This is the small matrix whose eigenvalues and modes are used to reconstruct the system.

### 4. It tests ranks 2, 5, and 10

These ranks are used because the article discusses how low rank can fail for translation-like behavior. A rank-2 model may not capture the translating part well, even though there are only two visible frequencies.

---

## How to understand the results

### `rank_error_summary.csv`

This gives the reconstruction error for each rank.

Lower error means DMD reconstructed the data more accurately.

Expected interpretation:

- Rank 2 should usually be the weakest.
- Rank 5 should improve.
- Rank 10 should improve more.

The reason is that DMD does not represent pure spatial translation very efficiently. Even though `f1` has one angular frequency, the moving spatial shape `sech(x+6-t)` requires multiple DMD modes.

### `svd_energy_table.csv`

This shows how much data energy is captured by each SVD rank.

High SVD energy does not always mean perfect DMD reconstruction. SVD tells you about data compression. DMD also needs time evolution to be represented well by the eigenvalues.

### `dmd_eigenvalue_frequency_table.csv`

Look at this column:

```text
omega_imag_angular_frequency
```

The important frequencies should be close to:

```text
2.3 and 2.8
```

But because the first component translates in space, DMD may introduce additional modes/eigenvalues to approximate it.

### `dmd_eigenvalues.png`

This plot shows the DMD eigenvalues in the complex plane.

For a stable oscillatory signal, eigenvalues should lie near the unit circle. Their angle is related to frequency.

### `dmd_mode_magnitudes.png`

This shows the spatial patterns of the DMD modes.

Modes are not the same as the original `f1` and `f2`. They are DMD's learned spatial structures that evolve exponentially/oscillatorily in time.

---

## Notes

The file `scripts/recreate_article_equation.py` has the grid and rank settings at the top:

```python
NX = 400
NT = 200
X_MIN = -10.0
X_MAX = 10.0
T_MIN = 0.0
T_MAX = 4.0 * np.pi
RANKS = [2, 5, 10]
```

If your article shows different grid values, change only these lines.
