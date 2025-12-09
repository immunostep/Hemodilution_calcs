from ._operations import load_fcs, getImf, custom_hdbscan, apply_threshold, gate_diagonal_percentile_df, setChannels
from typing import Tuple

def getHemodilutionEventsImf(fcsFilePath: str, plot = False) -> Tuple[int, float]:
    meta, df = load_fcs(fcsFilePath)
    cit = setChannels(meta)

    if cit not in ["calibur", "FACSCalibur"]:
        df_diagonal = gate_diagonal_percentile_df(df, plot=plot)
    else:
        df_diagonal = df.copy()
    df_threshold = apply_threshold(df_diagonal, plot=plot)
    final = custom_hdbscan(df_threshold, plot=plot)

    event_count = final.shape[0]
    imf = getImf(final)
    
    return event_count, imf