from typing import Tuple
import pandas as pd
from fcsparser import parse
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
import numpy as np
from sklearn.preprocessing import StandardScaler
import hdbscan
import re

channels = {
    "B8": "B8-A",
    "R1": "R1-A",
    "B4": "B4-A"
}


def plot_df_sns_with_hue(df, channels, hue_col=None, title=None):
    """
    Scatter plot con Seaborn, opcionalmente coloreando según una columna hue.

    Args:
        df: DataFrame con los datos
        channels: diccionario con claves 'B8' y 'R1'
        hue_col: columna para colorear puntos (ej. 'bead')
        filename: nombre del archivo donde se guarda el plot
        title: título del plot
    """
    plt.figure(figsize=(6,5))
    sns.scatterplot(
        data=df,
        x=channels["B8"],
        y=channels["R1"],
        hue=hue_col,
        palette={0:'blue', 1:'red'} if hue_col else None,
        s=3,
        alpha=0.7
    )
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(channels["B8"])
    plt.ylabel(channels["R1"])
    if title:
        plt.title(title)
    if hue_col:
        plt.legend(title=hue_col)
    plt.tight_layout()

def plot_custom(df: pd.DataFrame):
    # for idx, valor in enumerate(df[channels["B8"]]):
    #     print(f"Punto {idx}: {valor}")
    plot_df_sns_with_hue(df, channels)

def load_fcs(path):
    """
    Carga un FCS y devuelve un DataFrame con los eventos.
    """
    meta, data = parse(path, reformat_meta=True)
    for i, name in enumerate(meta['$P1N'] if '$P1N' in meta else []):
        meta[f'$P{i+1}N'] = re.sub(r'^[-_\s]+', '', name)
    return meta, data



def gate_diagonal_percentile_df(df: pd.DataFrame, xchan: str = "FSC-A", ychan: str = "FSC-H",
                                x_low_pct: float = 5, x_high_pct: float = 99,
                                epsilon: float = 0.25, y_offset_factor: float = 0.2,
                                plot: bool = False) -> pd.DataFrame:
    """
    Gating diagonal basado en percentiles.
    
    Args:
        df: DataFrame con los datos
        xchan, ychan: nombres de los canales
        x_low_pct, x_high_pct: percentiles para los extremos de x
        epsilon: ancho relativo de la diagonal
        y_offset_factor: desplazamiento vertical proporcional a x
        plot: si True, guarda un scatter plot mostrando los seleccionados
        filename: archivo donde se guarda el plot
    
    Returns:
        df_gated: DataFrame solo con los eventos dentro de la diagonal
    """

    if(plot):
        plot_df_sns_with_hue(df, channels)
    x = df[xchan]
    y = df[ychan]

    x_min = np.percentile(x, x_low_pct)
    x_max = np.percentile(x, x_high_pct)

    offset_min = y_offset_factor * x_min
    offset_max = y_offset_factor * x_max

    eps_start = epsilon * 2      # Ancho inicial
    eps_end   = epsilon          # Ancho final

    eps_x = eps_start + (x - x_min) * (eps_end - eps_start) / (x_max - x_min)

    y_lower = (1 - eps_x) * x - (offset_min + (x - x_min)*(offset_max-offset_min)/(x_max-x_min))
    y_upper = (1 + eps_x) * x - (offset_min + (x - x_min)*(offset_max-offset_min)/(x_max-x_min))
   
    # Selección de puntos dentro de la banda
    mask = (y >= y_lower) & (y <= y_upper)
    df_gated = df[mask].copy()

    if plot:
        df_copy = df.copy()
        df_copy['inside_gate'] = 0
        df_copy.loc[mask, 'inside_gate'] = 1

        plt.figure(figsize=(6,5))
        sns.scatterplot(data=df_copy, x=xchan, y=ychan, hue='inside_gate',
                        palette={0:'gray', 1:'red'}, alpha=0.7, s=5)
        plt.xlabel(xchan)
        plt.ylabel(ychan)
        plt.title("Diagonal Gate")
        plt.tight_layout()

    return df_gated

def apply_threshold(df: pd.DataFrame, plot = False) -> pd.DataFrame:
    """
    Aplica el threshold sobre channels["B8"]
    En x lo aplica igual que la app de escritorio y tambien aplica un threshold en y
    
    Args:
        :param df: Description
        :type df: pd.DataFrame
        :param plot: Description
    Returns:
        Dataframe con solo los eventos que pasan el threshold
    """

    channel_values = df[channels["B8"]].dropna()
    if channel_values.empty:
        valores_max = 1
    else:
        available = len(channel_values)
        if available >= 20:
            top_n = 10
        elif available >= 10:
            top_n = 5
        else:
            top_n = available
        top_n = max(1, top_n)
        # Promedia los eventos más intensos para evitar que outliers aíslen el umbral
        top_values = channel_values.nlargest(top_n)
        valores_max = top_values.mean()
        if not np.isfinite(valores_max) or valores_max <= 0:
            valores_max = channel_values.max()


    min_val = 1
    valores_max_lineal = np.log10(valores_max)
    valores_min_lineal = np.log10(min_val)

    corte = valores_min_lineal + 0.7 * (valores_max_lineal - valores_min_lineal)
    valor = 10 ** corte
    y_max = df[channels["R1"]].max()
    thresh_val_y = y_max * 0.05
    thresh_val = valor
    df_thresh = df[(df[channels["B8"]] > thresh_val) & (df[channels["R1"]] > thresh_val_y)].copy()

    if plot:
        plot_df_sns_with_hue(df, channels, title="Original Data")
        plot_df_sns_with_hue(df_thresh, channels, title="After Threshold")

    return df_thresh

def custom_hdbscan(df: pd.DataFrame, eps=0.2, min_samples=50, scale="log", plot = False) -> pd.DataFrame:
    """
    Aplica HDBSCAN clustering y retorna el cluster más grande.
    
    Args:
        df: DataFrame con columnas B8 y R1 (canales).
        eps: Parámetro de selección de clusters (por defecto 0.2).
        min_samples: Muestras mínimas para core point (por defecto 50).
        scale: Escala del plot ("log" o "linear", por defecto "log").
        plot: Si True, muestra visualización del clustering.
    
    Returns:
        DataFrame con solo los eventos del cluster más grande. Vacío si no hay clusters.
    """
    
    x_col = channels["B8"]
    y_col = channels["R1"]
    X = df[[x_col, y_col]].values
    X_scaled = StandardScaler().fit_transform(X)

    clusterer = hdbscan.HDBSCAN(min_cluster_size=min_samples,
                                min_samples=min_samples,
                                cluster_selection_epsilon=eps,
                                cluster_selection_method='eom',
                                allow_single_cluster=True)
    labels = clusterer.fit_predict(X_scaled)

    df_filtered = df.copy()
    df_filtered['cluster'] = labels

    # Seleccionar solo el cluster más grande, ignorando outliers (-1)
    cluster_counts = df_filtered[df_filtered['cluster'] != -1]['cluster'].value_counts()
    if cluster_counts.empty:
        print("No clusters detected, returning empty DataFrame")
        return pd.DataFrame(columns=df_filtered.columns)
    
    largest_cluster_label = cluster_counts.idxmax()
    df_largest = df_filtered[df_filtered['cluster'] == largest_cluster_label].copy()

    if plot:
        plt.figure(figsize=(6,5))
        sns.scatterplot(
            x=df_filtered[x_col],
            y=df_filtered[y_col],
            hue=df_filtered['cluster'],
            palette="tab10",
            alpha=0.3,
            s=5,
            legend=False
        )
        sns.scatterplot(
            x=df_largest[x_col],
            y=df_largest[y_col],
            color='red',
            s=10,
            alpha=0.9,
            label='Largest cluster'
        )
        plt.xscale(scale)
        plt.yscale(scale)
        plt.xlabel(x_col)
        plt.ylabel(y_col)
        plt.title("Largest Cluster after HDBSCAN")
        plt.tight_layout()

    return df_largest

def getImf(df: pd.DataFrame) -> float:
    """
    Calcula el IMF medio de todos los eventos de un DataFrame usando el canal B4.

    Args:
        df: DataFrame con los datos

    Returns:
        imf_medio: valor medio del IMF en el canal B4
    """

    imf_column = channels["B4"]
    if imf_column is None:
        raise ValueError("No se encontró el canal B4 en el diccionario 'channels'")

    if imf_column not in df.columns:
        raise ValueError(f"La columna '{imf_column}' no existe en el DataFrame")

    # for idx, valor in enumerate(df[imf_column]):
    #     print(f"Punto {idx}: {valor}")


    imf_medio = df[imf_column].median()
    return imf_medio

def _infer_channels_from_citometer(citometer: str) -> Tuple[str, str, str]:
    if citometer == "Aurora":
        return "B8-A", "R1-A", "B4-A"
    if citometer in ["Calibur", "FACSCalibur"]:
        return "FL3-Height", "FL4-Height", "FL2-Height"
    if citometer in ["Aria", "FACSAriaII"]:
        return "PerCP-A", "APC-A", "PE-A"
    if citometer in ["FACSCantoII"]:
        return "PerCP-Cy5-5-A", "APC-A", "PE-A"
    # TOLEDO
    # if citometer in ["BD FACSLyric"]:
    #     return "- _PerCP Cy5-5-A", "- _APC-A", "- _PE-A"
    if citometer in ["BD FACSLyric"]:
        return "PerCP-Cy5.5-A", "APC-A", "PE-A"
    return "", "", ""

def setChannelsManual(b8: str, r1: str, b4: str):
    channels["B8"] = b8
    channels["B4"] = b4
    channels["R1"] = r1

    return channels


def setChannels(meta) -> str:
    claves_intentadas = ['$CYT', '$CYTNAME', 'CYT', 'CYTOMETER', 'INSTRUMENT']

    # intentar claves directas (sensibles a mayúsculas según parse)
    citometro = None
    for k in claves_intentadas:
        if k in meta:
            citometro = meta[k]
            break

    # si no lo encontramos, buscar por substring en las claves
    if not citometro:
        for k, v in meta.items():
            lk = k.lower()
            if 'cyt' in lk or 'cytometer' in lk or 'instrument' in lk or 'instru' in lk:
                citometro = v
                break

    b8, r1, b4 = _infer_channels_from_citometer(citometro)
    channels["B8"] = b8
    channels["B4"] = b4
    channels["R1"] = r1

    return citometro
