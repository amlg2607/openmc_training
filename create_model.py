import openmc

print("--Iniciando la construcción del modelo openmc desde cero--")

#==========================
# PASO 1: DEFINICION DE MATERIALES Y BIBLIOTECAS NUCLEARES
#==========================
#explicación física:
#los materiales nucleares requieren especificar la densidad y la composición
# isotópica (fracciones atómicas o en peso). OpenMC busca las secciones eficaces
# correspondientes (secciones transversales de interacción neutrónica) en formato HD5
# utilizando la variable de entorno OPENMC_CROSS_SECTIONS.
# también se añade la ley de dispersión témica S((a,b)) para modelar cómo los
# neutrones térmicos interactuan con átomos ligados en moléculas (ej. H en H2O).

print("[1/4] configurando materiales...")

#-- Combustible: Dióxido de uránio (UO2) enriquecido al 2.4% en U-235---
fuel = openmc.Material(name="Uranium Dioxide Fuel")
fuel.add_nuclide("U234", 0.0003) #trazas de U-234
fuel.add_nuclide("U235", 0.024) #URÁNIO FISIBLE (2.4%)
fuel.add_nuclide("U238", 0.9757) #URÁNIO FERTIL (97.57%)
fuel.add_element("O", 2.0) #Oxígeno/estequimetria del óxido
fuel.set_density("g/cm3", 10.297) #densidad del UO2 sintetizada en pastilla

# Moderador y refrigerante: agua ligera (H2O 296 K)
water= openmc.Material(name="Light Water moderator")
water.add_element("H", 2.0) #hidrógeno
water.add_element("O", 1.0) #oxígeno
water.set_density("g/cm3", 0.9982) #densidad del agua a 296 K
# Aplicacion de S((a,b)) térmica crucial para la moderación correcta en hidrógeno ligero:
water.add_s_alpha_beta("c_H_in_H2O") #ley de dispersión térmica para H en H2O

# Agrupar los materiales en um contenedor oficial de OpenMC
materials = openmc.Materials([fuel, water])
materials.export_to_xml() #exportar a XML para que OpenMC lo lea
print("-- Archivo 'materials.xml'generado correctamente!")

#==================================================
# PASO 2: DEFINICION DE GEOMETRIA Y CONDICIONES DE FRONTERA
#==================================================
# Explicación física:
# La geometria se construye mediante superficies analíticas (cilindros, planos).
# Las regiones espaciales se definen mediante operaciones lógicas de conjuntos 
# (+ para el exterior da superfície, - para el interior).
# Las condiciones de frontera ('boundary_type') determinan que le ocurre al neutrón
# al alcanzar los límites: 'vacuum'(escapa y se pierde) o 'reflective' (rebota)

print("[2/4] construyendo geometria y superficies...")
# -----superficies-------
# Cilindro interior que delimita la pastilla de combustible (radio = 0.42 cm)
fuel_outer_radius = openmc.Cylinder(r=0.42, name="Fuel Outer Surface")
# Cilindro exterior que delimita el canal del moderador/celda unitária radio=1.26 cm
# se usa frontera reflectiva para simular un arreglo infinito repetido de celdas
pitch_radius = openmc.Cylinder(r=1.26, name="Unit Cell Boundary", boundary_type="reflective")

#planos axiales superiores e inferiores para limitar la altura del núcleo (Z=-50 cm a +50 cm)
# Aquí se configura la condición de vacío en las tapas superior e inferior
z_lower = openmc.ZPlane(z0=-50.0, name="Lower Boundary", boundary_type="vacuum")
z_upper = openmc.ZPlane(z0=50.0, name="Upper Boundary", boundary_type="vacuum")

#------Regiones (áreas especiales limitadas por superficies)-------
# Región de combustible: dentro del cilindro de combustible y entre los planos axiales
fuel_region = -fuel_outer_radius & +z_lower & -z_upper
# Región de moderador: entre el cilindro de combustible y el radio de la celda,confinado axialmente
water_region = +fuel_outer_radius & -pitch_radius & +z_lower & -z_upper

#------Celdas (Asociaciónde regiones con sus respectivos materiales) -------
fuel_cell = openmc.Cell(name="Fuel Cell", fill=fuel)
fuel_cell.region = fuel_region
water_cell = openmc.Cell(name="Water Moderator Cell", fill=water)
water_cell.region = water_region

#------Universo y geometria raiz------
root_universe = openmc.Universe(cells=[fuel_cell, water_cell], name="Root Universe")
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml() #exportar a XML para que OpenMC lo lea
print("-- Archivo 'geometry.xml' generado correctamente!")

#========================================================
# PASO 3: PARÁMETROS DE SIMULACIÓN Y FUENTE DE NEUTRONES (SETTINGS)
#========================================================
# Explicación física:
# OpenMC utiliza el metodo Monte Carlo de k-eigenvalor para problemas críticos.
# Se configuran lotes (batches), generaciones de neutrones y el espacio de la
# fuente inicial de fisión distribuida uniformemente dentro del volumen de combustible.

print("[3/4] configurando parámetros de Monte Carlo (settings)...")

settings = openmc.Settings()
settings.batches = 100 #número total de lotes (batches)
settings.inactive = 10 #lotes inactivos (no se cuentan para estabilizar la distribuición de la fuente)
settings.particles = 10000 #número de partículas por lote

# definicion espacial de la fuente inicial de neutrones (acotada al volumen del combustible)
lower_left = (-0.42, -0.42, -50.0) #coordenadas de la esquina inferior izquierda del cilindro
upper_right = (0.42, 0.42, 50.0) #coordenadas de la esquina superior derecha del cilindro
source_box = openmc.stats.Box(lower_left, upper_right, only_fissionable=True) #fuente solo en combustible
settings.source = openmc.Source(space=source_box) #asignar la fuente al settings
settings.export_to_xml() #exportar a XML para que OpenMC lo lea
print("-- Archivo 'settings.xml' generado correctamente!")

# =============================================
# PASO 4: TALLIES (REGISTROS DE MEDICIÓN FÍSICA)
# =============================================
# Explicación física:
# Los tallies permiten estimar magnitudes físicas como flujos neutronicos, tasas
# de reacción o potencias en zonas especificas del modelo. Para este tutorial
# inicial básico, dejamos la estructura limpia de tallies

print("[4/4] configurando tallies (registros de medición)...")
tallies = openmc.Tallies() #inicializar contenedor de tallies

# Definir un filtro espacial (para medir específicamente en la celda de combustible)
cell_filter = openmc.CellFilter([fuel_cell])

# Definir un filtro de energia con grupos estándar (ej. límites térmico, epitérmico y rápido)
# Fronteras en eV: Térmico (<0.625 eV), Epitérmico(0.625 eV - 100 keV), Rápido (>100 keV)
energy_filter = openmc.EnergyFilter([0.0, 0.625, 1.0e5, 20.0e6])

# Crear el Tally de Flujo
flux_tally = openmc.Tally(name="Fuel Flux Tally")
flux_tally.filters = [cell_filter,energy_filter]
flux_tally.scores = ['flux', 'fission'] #Quiero medir el flujo y tasa de fisiones

# Añadir el Tally de Flujo
tallies.append(flux_tally)

tallies.export_to_xml() #exportar a XML para que OpenMC lo lea
print("-- Archivo 'tallies.xml' generado correctamente!")

print("--Modelo creado exitosamente con el Tally de flujo y fisiones!--")

# =========================================================================
# PASO 5: VISUALIZACIÓN GEOMÉTRICA (PLOTS)
# =========================================================================
print("[5/5] Configurando gráficos de geometría...")

# Creamos un objeto de ploteo 2D
plot = openmc.Plot()
plot.filename = 'geometry_plot'
plot.basis = 'xy'                      # Corte en el plano XY
plot.width = (3.0, 3.0)                # Ancho y alto de la ventana de visualización en cm
plot.pixels = (400, 400)               # Resolución de la imagen
plot.color_by = 'material'             # Colorear por material asignado
plot.colors = {
    fuel: 'red',                       # El combustible se verá rojo
    water: 'blue'                      # El moderador se verá azul
}

# Agrupamos y exportamos a plots.xml
plots = openmc.Plots([plot])
plots.export_to_xml()
print(" -> Archivo 'plots.xml' generado correctamente.")
