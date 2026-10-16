import openmc
import glob

# Buscar autimaticamente el archivo statepoint generado en la carpeta
statepoint_files = glob.glob("statepoint.*.h5")
if not statepoint_files:
    print("No se encontró ningún archivo statepoin.h5. Corrite la simulación?")
    exit()

sp_filename = statepoint_files[0]
print(f"Leyendo archivo de resultados: {sp_filename}\n")

# Abrir rl archivo statepoint
with openmc.StatePoint(sp_filename) as sp:
    
    # 1. Obtener el k-effective combinado
    k_combined = sp.k_combined
    print(f"Factor de multiplicación (k_eff): {k_combined.nominal_value:.5f} +/- {k_combined.std_dev: .5f}")

    # 2. Obtener el Tally que definimos (por su nombre o ID)
    # Open MC asigna nombres a los tallies si se los damos
    tally = sp.get_tally(name="Fuel Flux Tally")

    print("\n--Resultados del Tally: Flujo y Fisiones en el combustible---")
    print(tally.mean)


    