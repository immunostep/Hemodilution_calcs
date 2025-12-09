# Hemodilution_calcs

Proyecto para procesamiento y gating de ficheros FCS (citometría de flujo). Incluye utilidades para cargar ficheros FCS, definir gates (diagonal, densidad, HDBSCAN, FlowCal), calcular IMF (canal `B4`) y exportar/visualizar resultados. Especializado para el kit de hemodilución

## Estructura del repositorio

- `hemodiluton_calcs/` - Código principal.
  - `hemodiluton_calcs/main.py` - Punto de entrada.
  - `hemodiluton_calcs/_operations/` - Operaciones de procesamiento y gating.

## Requisitos

- Python 3.12+ recomendado. (Usado para pruebas: 3.12.12)
- Dependencias:
  - `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`
  - `fcsparser`
  - `scikit-learn`, `hdbscan`
  - `FlowCal`

Revisa `/requirements.txt` para dependencias específicas del subpaquete.

## Instalación (uso local desde otro proyecto)

Si vas a usar esta librería en otro proyecto localmente, instala directamente desde la ruta del proyecto. Hay dos opciones comunes:

- Instalar desde la ruta del repositorio (no editable):

```powershell
# desde el proyecto consumidor
pip install "ruta/a/hemodilution_calcs"
```

- Instalar en modo editable (recomendado durante desarrollo):

```powershell
# desde el proyecto consumidor
pip install -e "ruta/a/hemodilution_calcs"
```

- Otra manera seria estar en la carpeta directamente
```powershell
cd "ruta/a/hemodilution_calcs"
pip install .

# O de manera editable
pip install -e .
```

Instalarlo de manera editable hace que:
- Se cree un enlace simbólico, es decir, los cambios en el codigo estarán aplicados directamente (ideal para desarrollo)

Si no se instala editable, pip empaqueta el proyecto y lo instala (se copia)

## Uso como librería (ejemplo)

Una vez instalado, importa el paquete/funciones desde Python:

```python
from hemodilution_calcs import getHemodilutionEventsImf

eventos, imf = getHemodilutionEventsImf("data/1_Aurora.fcs")

print("Eventos detectados:", eventos)
print("IMF:", imf)
```
La librería expone únicamente la función de alto nivel getHemodilutionEventsImf, encargada de:
- Cargar el archivo .fcs
- Aplicar los filtros y el gating
- Ejecutar el clustering HDBSCAN
- Calcular el IMF final
