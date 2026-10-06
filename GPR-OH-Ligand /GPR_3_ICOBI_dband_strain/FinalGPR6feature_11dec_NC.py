#!/usr/bin/env python3

from pathlib import Path
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from numpy.linalg import LinAlgError
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, RBF
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import LeaveOneOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.inspection import permutation_importance
import matplotlib as mpl
import shap

#global setting of fonts for math formula
mpl.rcParams.update({
    "mathtext.fontset": "stix",
    "font.family": "STIXGeneral",
})

TICK_FONTSIZE = 14
AXIS_LABEL_FONTSIZE = 14
HIST_COLOR = "steelblue"

LABEL = (
    r"$\boldsymbol{\Delta\Delta}\,"
    r"\mathbf{G}_{\mathbf{OH}^{\boldsymbol{*}}}^{\mathbf{Ligand\ effect}}$ (eV)"
)

#read csv file

INFILE = "cleaned_6feature_Ligand_OH_tablebased_SI.csv"
FEATURE_COLS = ["ICOBI-L1", "d-band center", "Strain"]

TARGET_COL   = "E_OH(ligand)"
META_COLS    = ["Metal", "Support", "Termination", "Number of Metal Layers"]

SEED     = 42
RESTARTS = 10
ALPHA    = 0.1

OUTDIR = Path("final_clean_output")
OUTDIR.mkdir(exist_ok=True)

warnings.filterwarnings("ignore", category=UserWarning)
plt.style.use("default")
plt.rcParams.update({"axes.grid": False})

#for different color in plot belonging to different metal
def metal_color(name):
    s = str(name).lower()
    if "au" in s: return "red"
    if "ag" in s: return "skyblue"
    if "pt" in s: return "green"
    return "gray"

#parity plot
def parity_plot(ax, y_true, y_pred, metals):
    for metal in np.unique(metals):
        idx = (metals == metal)
        ax.scatter(
            y_true[idx], y_pred[idx],
            label=metal,
            color=metal_color(metal),
            edgecolors="k",
            alpha=0.85, s=70
        )

    lo = float(min(y_true.min(), y_pred.min())) - 0.1
    hi = float(max(y_true.max(), y_pred.max())) + 0.1

    ax.plot([lo, hi], [lo, hi], "k--", lw=1.5)

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("DFT " + LABEL, fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
    ax.set_ylabel("Predicted " + LABEL, fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
    ax.tick_params(labelsize=TICK_FONTSIZE)
    ax.legend(title="Metal", fontsize=12, title_fontsize=12)

#histogram plot
def residual_hist(ax, values, title):
    bins = np.arange(-0.8, 0.8 + 0.2, 0.2)
    ax.hist(values, bins=bins, color=HIST_COLOR, edgecolor='k', alpha=0.9)

    ax.axvline(-0.2, color='k', linestyle='--', lw=1.5)
    ax.axvline(+0.2, color='k', linestyle='--', lw=1.5)
    ax.axvline(0.0,  color='k', linestyle=':',  lw=1.5)

    ax.set_xlim(-0.8, 0.8)
    ax.set_xlabel("Residual (true − pred), eV", fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
    ax.set_ylabel("Count", fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")
    ax.tick_params(labelsize=TICK_FONTSIZE)
    ax.set_title(title, fontsize=AXIS_LABEL_FONTSIZE, fontweight="bold")

#MEtrics
def metrics(y, yp):
    return (
        float(mean_absolute_error(y, yp)),
        # float(mean_squared_error(y, yp, squared=False)),
        float(mean_squared_error(y, yp)),
        float(r2_score(y, yp))
    )
#GPR MODEL
def make_pipeline(alpha):
    kernel = ConstantKernel(1.0, (1e-2, 1e3)) * RBF(1.0, (1e-3, 1e3))
    return Pipeline([
        ("scaler", StandardScaler()),
        ("gpr", GaussianProcessRegressor(
            kernel=kernel,
            alpha=alpha,
            normalize_y=True,
            n_restarts_optimizer=RESTARTS,
            random_state=SEED
        ))
    ])


def main():

    df = pd.read_csv(INFILE, skipinitialspace=True)
    df.columns = df.columns.str.strip()
    df = df.dropna(subset=[TARGET_COL]).reset_index(drop=True)

    X = df[FEATURE_COLS].fillna(df[FEATURE_COLS].median(numeric_only=True)).to_numpy()
    y = df[TARGET_COL].astype(float).to_numpy()
    metals = df["Metal"].astype(str).to_numpy()

    n = len(y)

    # Insample GPR model 
    model = make_pipeline(ALPHA)
    model.fit(X, y)

    Xs = model.named_steps["scaler"].transform(X)
    y_pred_in, y_std_in = model.named_steps["gpr"].predict(Xs, return_std=True)
    res_in = y - y_pred_in

#OUtput in csv for insample
    pd.DataFrame({
        "y_true": y,
        "y_pred": y_pred_in,
        "y_std": y_std_in,
        "residual": res_in
    }).to_csv(OUTDIR / "predictions_in_sample.csv", index=False)

#parityplot for insample
    fig, ax = plt.subplots(figsize=(6, 6))
    parity_plot(ax, y, y_pred_in, metals)
    fig.tight_layout()
    fig.savefig(OUTDIR / "parity_in_sample.png", dpi=300)
    plt.close()

#reidual histogram for insample
    fig, ax = plt.subplots(figsize=(7, 4))
    residual_hist(ax, res_in, "Residual Histogram (In-sample)")
    fig.tight_layout()
    fig.savefig(OUTDIR / "hist_in_sample.png", dpi=300)
    plt.close()

    mae_in, rmse_in, r2_in = metrics(y, y_pred_in)


    # LOOCV 
    loo = LeaveOneOut()
    y_pred_loo = np.full(n, np.nan)

    for train_idx, test_idx in loo.split(X):
        fold = make_pipeline(ALPHA)
        fold.fit(X[train_idx], y[train_idx])
        Xt = fold.named_steps["scaler"].transform(X[test_idx])
        m = fold.named_steps["gpr"].predict(Xt)
        y_pred_loo[test_idx] = float(m[0])

    res_loo = y - y_pred_loo

#loo output saved 
    pd.DataFrame({
        "y_true": y,
        "y_pred": y_pred_loo,
        "residual": res_loo
    }).to_csv(OUTDIR / "predictions_loocv.csv", index=False)

#parityplot for LOocv
    fig, ax = plt.subplots(figsize=(6, 6))
    parity_plot(ax, y, y_pred_loo, metals)
    fig.tight_layout()
    fig.savefig(OUTDIR / "parity_loocv.png", dpi=300)
    plt.close()

#histogram for loocv
    fig, ax = plt.subplots(figsize=(7, 4))
    residual_hist(ax, res_loo, "Residual Histogram (LOOCV)")
    fig.tight_layout()
    fig.savefig(OUTDIR / "hist_loocv.png", dpi=300)
    plt.close()

    mae_loo, rmse_loo, r2_loo = metrics(y, y_pred_loo)


    # PFI on full datatse
    pfi = permutation_importance(
        model, X, y,
        n_repeats=20,
        scoring="neg_mean_absolute_error",
        random_state=SEED
    )

    df_pfi = (
        pd.DataFrame({
            "Feature": FEATURE_COLS,
            "PermutationFeatureImportance": np.abs(pfi.importances_mean)
        })
        .sort_values("PermutationFeatureImportance", ascending=False)
    )

    df_pfi.to_csv(OUTDIR / "pfi_importances.csv", index=False)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.barh(
        df_pfi["Feature"],
        df_pfi["PermutationFeatureImportance"],
        color="darkorange",
        edgecolor="k"
    )

    ax.set_xlabel(
        "Permutation Feature Importance",
        fontsize=AXIS_LABEL_FONTSIZE,
        fontweight="bold"
    )
    ax.set_ylabel(
        "Feature",
        fontsize=AXIS_LABEL_FONTSIZE,
        fontweight="bold"
    )

    ax.tick_params(labelsize=TICK_FONTSIZE)
    ax.invert_yaxis()

    fig.tight_layout()
    fig.savefig(OUTDIR / "pfi_importances.png", dpi=300)
    plt.close()


    # SHAP
    try:
        rng = np.random.default_rng(SEED)

        bg_size = min(100, len(X))
        ex_size = min(200, len(X))

        bg_idx = rng.choice(len(X), size=bg_size, replace=False)
        ex_idx = rng.choice(len(X), size=ex_size, replace=False)

        # background & explanation sets
        X_bg = X[bg_idx].astype(float, copy=False)
        X_s  = X[ex_idx].astype(float, copy=False)

        # prediction function
        f_predict = lambda data: model.predict(data)

        # SHAP explainer
        explainer = shap.KernelExplainer(f_predict, X_bg, link="identity")

        # SHAP values for the explanation subset
        shap_values = explainer.shap_values(X_s, nsamples=200)

        # ---- SAVE RAW SHAP VALUES FOR ALL EXPLAINED SAMPLES ----
        df_shap_full = pd.DataFrame(
            shap_values,
            columns=FEATURE_COLS
        )
        df_shap_full.insert(0, "sample_index", ex_idx)
        df_shap_full.to_csv(OUTDIR / "shap_values_full.csv", index=False)

        # mean absolute SHAP
        shap_abs = np.mean(np.abs(shap_values), axis=0)

        (pd.DataFrame({"feature": FEATURE_COLS, "mean_abs_shap": shap_abs})
           .sort_values("mean_abs_shap", ascending=False)
           .to_csv(OUTDIR / "shap_importances.csv", index=False))

        # SHAP beeswarm plot
        shap.summary_plot(
            shap_values, X_s, feature_names=FEATURE_COLS, show=False
        )
        plt.tight_layout()
        plt.savefig(OUTDIR / "shap_beeswarm.png", dpi=300)
        plt.close()

    except Exception as e:
        (OUTDIR / "shap_importances.csv").write_text(
            "SHAP failed\n" + str(e) + "\n"
        )


    #  Summary
    (OUTDIR / "summary.txt").write_text(
        f"In-sample: MAE={mae_in:.6f}, RMSE={rmse_in:.6f}, R2={r2_in:.6f}\n"
        f"LOOCV:     MAE={mae_loo:.6f}, RMSE={rmse_loo:.6f}, R2={r2_loo:.6f}\n"
    )

    print("DONE. Saved outputs in:", OUTDIR)


if __name__ == "__main__":
    main()
