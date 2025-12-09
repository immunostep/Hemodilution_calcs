import os
import pytest
import json
from hemodilution_calcs.main import getHemodilutionEventsImf

TEST_FILES_DIR = os.path.join(os.path.dirname(__file__), "test_files")
RESULTS_JSON = os.path.join(TEST_FILES_DIR, "values.json")

# Cargar resultados esperados
with open(RESULTS_JSON, "r") as f:
    resultados_esperados = json.load(f)

# Diccionario para acceso rápido por nombre de fichero
esperados_dict = {r["file"]: r for r in resultados_esperados}

# Listado de ficheros .fsc
fsc_files = [f for f in os.listdir(TEST_FILES_DIR) if f.endswith(".fcs")]

# Si no hay ficheros .fsc, pytest lo marcará como SKIPPED dentro del test
@pytest.mark.parametrize("fsc_file", fsc_files)
def test_algoritmo(fsc_file):
    if not fsc_file:
        pytest.skip("No hay ficheros .fsc para testear")

    file_path = os.path.join(TEST_FILES_DIR, fsc_file)
    
    # Ejecutamos el algoritmo
    num_eventos, imf = getHemodilutionEventsImf(file_path)
    
    # Obtenemos los valores esperados del JSON
    esperados = esperados_dict.get(fsc_file)
    assert esperados is not None, f"No hay resultados esperados para {fsc_file}"

    #Test opcional, comprobar los numeros de eventos.
    
    # Comprobamos IMF con tolerancia del 5%
    tolerancia = 0.05
    diferencia = abs(imf - esperados["IMF"])
    assert diferencia / esperados["IMF"] <= tolerancia, f"{fsc_file}: IMF fuera de tolerancia"
