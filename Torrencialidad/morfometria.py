# %%

import joblib
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.svm import SVC
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.preprocessing import StandardScaler


def train_and_plot_models(df, combinations):
    """
    Entrena un modelo SVM para clasificar cuencas torrenciales basado en
    distintas combinaciones de características y genera gráficos de las
    predicciones.

    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame que contiene las características de las cuencas y la
        variable objetivo 'FFR' que indica si la cuenca es torrencial (1)
        o no (0).

    combinations : list of tuple
        Lista de tuplas que contienen las combinaciones de características
        (columnas) a ser utilizadas para entrenar los modelos SVM.

    Returns:
    --------
    models : dict
        Diccionario que contiene los modelos SVM entrenados. Las claves son
        índices numéricos que corresponden a las combinaciones en el orden
        de entrada.

    scaler : StandardScaler
        Objeto `StandardScaler` ajustado con los datos de entrada, utilizado
        para normalizar las características antes de entrenar los modelos.
    """

    df["FFR"] = (df["FFR"] == "T").astype(int)
    X = df[["A (sq. km)", "Su", "M", "Rh"]]
    y = df["FFR"]
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    models = {}
    fig, axes = plt.subplots(1, len(combinations), figsize=(18, 6))
    sns.set(style="whitegrid")

    for i, ((x_col, y_col), ax) in enumerate(zip(combinations, axes)):
        X_plot = X_scaled[:, [X.columns.get_loc(x_col), X.columns.get_loc(y_col)]]
        model = SVC(kernel="linear", C=10, gamma=0.1)
        model.fit(X_plot, y)
        models[i] = model

        xx, yy = np.meshgrid(
            np.linspace(X_plot[:, 0].min() - 1, X_plot[:, 0].max() + 1, 100),
            np.linspace(X_plot[:, 1].min() - 1, X_plot[:, 1].max() + 1, 100),
        )
        Z = model.predict(np.c_[xx.ravel(), yy.ravel()])
        Z = Z.reshape(xx.shape)

        color_torrencial = "darkorange"
        color_notorrencial = "cornflowerblue"
        cmap_custom = ListedColormap([color_notorrencial, color_torrencial])

        ax.contourf(xx, yy, Z, alpha=0.3, cmap=cmap_custom)
        ax.scatter(
            X_plot[y == 0, 0],
            X_plot[y == 0, 1],
            c=color_notorrencial,
            label="No torrencial",
            edgecolor="k",
            s=100,
            alpha=0.8,
        )
        ax.scatter(
            X_plot[y == 1, 0],
            X_plot[y == 1, 1],
            c=color_torrencial,
            label="Torrencial",
            edgecolor="k",
            s=100,
            alpha=0.8,
        )

        x_col_idx = X.columns.get_loc(x_col)
        y_col_idx = X.columns.get_loc(y_col)

        if i == 0:  # First plot
            x_ticks_original = np.arange(0, 1001, 200)
        elif i == 1:  # Second plot
            x_ticks_original = np.arange(2, 6, 1)
        else:  # Third plot
            x_ticks_original = np.arange(0, 0.81, 0.2)

        y_ticks_original = np.arange(0, 0.31, 0.1)

        x_ticks_scaled = scaler.transform(
            np.column_stack(
                [
                    (
                        x_ticks_original
                        if j == x_col_idx
                        else np.zeros_like(x_ticks_original)
                    )
                    for j in range(X.shape[1])
                ]
            )
        )[..., x_col_idx]
        y_ticks_scaled = scaler.transform(
            np.column_stack(
                [
                    (
                        y_ticks_original
                        if j == y_col_idx
                        else np.zeros_like(y_ticks_original)
                    )
                    for j in range(X.shape[1])
                ]
            )
        )[..., y_col_idx]

        ax.set_xticks(x_ticks_scaled)
        ax.set_xticklabels(
            [f"{x:.1f}" if x < 1 else f"{x:.0f}" for x in x_ticks_original], fontsize=18
        )
        ax.set_yticks(y_ticks_scaled)
        ax.set_yticklabels([f"{y:.1f}" for y in y_ticks_original], fontsize=18)

        ax.set_xlabel(x_col, fontsize=20)
        ax.set_ylabel(y_col, fontsize=20)
        ax.tick_params(axis="both", which="major", labelsize=18)

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles, labels, loc="center", bbox_to_anchor=(0.5, -0.08), fontsize=20, ncol=2
    )
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

    return models, scaler


df = pd.read_csv("torrencialidad.csv")
combinations = [("A (sq. km)", "Rh"), ("Su", "Rh"), ("M", "Rh")]
models, scaler = train_and_plot_models(df, combinations)


joblib.dump(models, "morfo-arango.pkl")
joblib.dump(scaler, "scaler-arango.pkl")
