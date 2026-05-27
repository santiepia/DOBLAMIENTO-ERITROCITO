import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from matplotlib.backends.backend_pdf import PdfPages

# ============================================================================
# 1. DEFINIR LA FUNCIÓN Y(x) de Cassini
# ============================================================================
def Y_cassini(x, A, B, C):
    """
    Función de Cassini modificada para el perfil de la figura
    y(x) = B * sqrt( sqrt(C^4 + 4*A^2*x^2) - A^2 - x^2 )
    """
    interior = np.sqrt(C**4 + 4.0 * A**2.0 * x**2.0) - A**2.0 - x**2.0
    
    # Si el interior es negativo, la función no está definida (y=0)
    with np.errstate(invalid='ignore'):
        y = B * np.sqrt(np.maximum(interior, 0))
    
    return y


# ============================================================================
# 1.5 FUNCIÓN PARA CALCULAR LA ENERGÍA DE DOBLAMIENTO
# ============================================================================
def bending_energy(A, B, C, N=500):
    """
    Calcula la energía usando malla no uniforme:
    - Más puntos cerca del borde (donde la derivada es grande)
    - Menos puntos en el centro (donde la función es suave)
    """
    x_max = np.sqrt(A**2 + C**2)
    if x_max <= 1e-10:
        return 0.0
    
    # Malla no uniforme: concentrar puntos cerca de x_max
    # Usamos transformación: x = x_max * (t)^p
    # p > 1 concentra puntos cerca de x_max
    # p < 1 concentra puntos cerca de 0
    
    if B < 0.6:  # Formas con borde afilado (flying saucer)
        p = 2.5  # Concentración extra cerca del borde
    else:
        p = 1.0  # Malla uniforme
    
    t = np.linspace(0, 1, N+1)
    x = x_max * (t ** p)
    
    # dx no es constante, necesitamos los diferenciales para Simpson
    # Para malla no uniforme, usamos Simpson compuesto en cada subintervalo
    
    # Calcular y(x)
    y = np.array([Y_cassini(xi, A, B, C) for xi in x])
    
    # Derivada numérica (diferencias centrales en malla no uniforme)
    yp = np.zeros(N+1)
    for i in range(1, N):
        if y[i] > 0:
            # Fórmula para malla no uniforme
            dx_left = x[i] - x[i-1]
            dx_right = x[i+1] - x[i]
            yp[i] = ( (dx_left**2) * (y[i+1] - y[i]) + 
                      (dx_right**2) * (y[i] - y[i-1]) ) / (dx_left * dx_right * (dx_left + dx_right))
    
    yp[0] = (y[1] - y[0]) / (x[1] - x[0])
    yp[N] = (y[N] - y[N-1]) / (x[N] - x[N-1])
    
    # Integración con malla no uniforme (regla del trapecio generalizada para simplicidad)
    E = 0.0
    for i in range(N):
        if y[i] <= 1e-10 or y[i+1] <= 1e-10:
            continue
        
        # Aproximación del integrando en cada subintervalo (promedio)
        # Cálculo de radios en ambos extremos (simplificado)
        
        # Para no complicar, usamos valores en el punto medio
        x_mid = (x[i] + x[i+1]) / 2
        y_mid = Y_cassini(x_mid, A, B, C)
        
        if y_mid <= 1e-10:
            continue
        
        # Derivada en el punto medio
        dx_mid = (x[i+1] - x[i]) / 2
        yp_mid = (Y_cassini(x_mid + dx_mid, A, B, C) - Y_cassini(x_mid - dx_mid, A, B, C)) / (2*dx_mid)
        
        # Segunda derivada (aproximada)
        ypp_mid = (Y_cassini(x_mid + dx_mid, A, B, C) - 2*y_mid + Y_cassini(x_mid - dx_mid, A, B, C)) / (dx_mid**2)
        
        # Radios
        if abs(ypp_mid) > 1e-10:
            R1 = (1 + yp_mid**2)**1.5 / abs(ypp_mid)
        else:
            R1 = 1e10
        
        if abs(yp_mid) > 1e-10:
            R2 = y_mid * np.sqrt(1 + yp_mid**2) / abs(yp_mid)
        else:
            R2 = y_mid
        
        curv_media = 1.0/R1 + 1.0/R2
        integrando_mid = (curv_media**2) * y_mid * np.sqrt(1 + yp_mid**2)
        
        # Contribución del subintervalo (regla del trapecio adaptada)
        E += integrando_mid * (x[i+1] - x[i])
    
    E = 2.0 * np.pi * E
    
    return E

# ============================================================================
# 2. CREAR DIRECTORIO PARA GUARDAR LAS GRÁFICAS
# ============================================================================
os.makedirs("graficas_filtradas", exist_ok=True)

# ============================================================================
# 3. LEER EL ARCHIVO FILTRADO
# ============================================================================
print("Leyendo archivo filtrado...")

# Intentar leer el archivo (formato tabulado o CSV)
try:
    # Primero intentar leer como CSV
    df = pd.read_csv('filtro.csv')
    print("Archivo leído como CSV")
except:
    try:
        # Si no, leer como texto tabulado
        df = pd.read_csv('filtro.txt', sep='\t')
        print("Archivo leído como TXT tabulado")
    except:
        # Si falla, leer formato de Fortran
        df = pd.read_csv('filtro.txt', delim_whitespace=True, comment='#')
        print("Archivo leído con formato whitespace")

print(f"Total de combinaciones filtradas: {len(df)}")

# ============================================================================
# 3.5 CALCULAR LA ENERGÍA PARA CADA COMBINACIÓN
# ============================================================================
print("\nCalculando energía de doblamiento para cada forma...")

energias = []
for idx, row in df.iterrows():
    A = row['A']
    B = row['B']
    C = row['C']
    
    E = bending_energy(A, B, C)
    energias.append(E)
    
    if (idx + 1) % 10 == 0:
        print(f"  Procesadas {idx + 1} / {len(df)} formas...")

df['Energia'] = energias

# Guardar DataFrame con la nueva columna de energía
df.to_csv('filtro_con_energia.csv', index=False)
print(f"\n✅ Energías calculadas. Archivo guardado: 'filtro_con_energia.csv'")

# ============================================================================
# 4. PARÁMETROS DE LA GRÁFICA
# ============================================================================
NUM_PUNTOS = 1000  # Resolución de la curva

# ============================================================================
# 5. NUEVA GRÁFICA: ENERGÍA vs. B
# ============================================================================
print("\nGenerando gráfica: Energía de doblamiento vs. Parámetro B...")

plt.figure(figsize=(10, 7))

# Crear scatter plot coloreado por A (o por C)
scatter = plt.scatter(df['B'], df['Energia'], 
                      c=df['A'], cmap='viridis', 
                      s=30, alpha=0.7, edgecolors='black', linewidth=0.5)

# Barra de color para A
cbar = plt.colorbar(scatter)
cbar.set_label('Parámetro A', fontsize=12)

plt.xlabel('Parámetro B (altura/escala vertical)', fontsize=14)
plt.ylabel('Energía de doblamiento (E)', fontsize=14)
plt.title('Energía de curvatura vs. Parámetro B\n(formas que cumplen condiciones fisiológicas)', fontsize=14)
plt.grid(True, alpha=0.3)

# Guardar gráfica
plt.savefig('graficas_filtradas/energia_vs_B.png', dpi=150, bbox_inches='tight')
plt.savefig('graficas_filtradas/energia_vs_B.pdf', bbox_inches='tight')
plt.close()

print("✅ Gráfica 'Energía vs. B' guardada en 'graficas_filtradas/energia_vs_B.png'")

# ============================================================================
# 5.5 GRÁFICA ADICIONAL: ENERGÍA vs. B (con regresión o líneas de tendencia)
# ============================================================================
print("\nGenerando gráfica adicional: Energía vs. B (con líneas por valor de A)...")

plt.figure(figsize=(12, 8))

# Agrupar por A y graficar líneas
valores_A = sorted(df['A'].unique())
colores = plt.cm.plasma(np.linspace(0, 1, len(valores_A)))

for i, A_val in enumerate(valores_A):
    subset = df[df['A'] == A_val]
    if len(subset) > 1:
        # Ordenar por B para la línea
        subset_ordenado = subset.sort_values('B')
        plt.plot(subset_ordenado['B'], subset_ordenado['Energia'], 
                 'o-', color=colores[i], linewidth=1.5, markersize=4, 
                 label=f'A = {A_val:.2f}', alpha=0.8)

plt.xlabel('Parámetro B', fontsize=14)
plt.ylabel('Energía de doblamiento (E)', fontsize=14)
plt.title('Energía de curvatura vs. Parámetro B (agrupado por A)', fontsize=14)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
plt.grid(True, alpha=0.3)
plt.tight_layout()

plt.savefig('graficas_filtradas/energia_vs_B_por_A.png', dpi=150, bbox_inches='tight')
plt.close()

print("✅ Gráfica 'Energía vs. B (por A)' guardada")

# ============================================================================
# 6. GENERAR UN PDF CON TODAS LAS GRÁFICAS DE FORMAS
# ============================================================================
print("\nGenerando PDF con todas las gráficas de formas...")

with PdfPages('graficas_filtradas/todas_formas.pdf') as pdf:
    for idx, row in df.iterrows():
        A = row['A']
        B = row['B']
        C = row['C']
        E = row['Energia']
        
        b_limite = np.sqrt(A**2 + C**2)
        x = np.linspace(-b_limite, b_limite, NUM_PUNTOS)
        y = Y_cassini(x, A, B, C)
        
        plt.figure(figsize=(8, 6))
        plt.plot(x, y, 'b-', linewidth=2)
        plt.plot(x, -y, 'r-', linewidth=2)
        plt.fill_between(x, -y, y, color='lightblue', alpha=0.5)
        
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(f'A={A:.4f}, B={B:.4f}, C={C:.4f}\nEnergía = {E:.4e}')
        plt.axhline(y=0, color='black', linewidth=0.5)
        plt.axvline(x=0, color='black', linewidth=0.5)
        plt.grid(True, alpha=0.3)
        plt.axis('equal')
        
        pdf.savefig()
        plt.close()
        
        if (idx + 1) % 20 == 0:
            print(f"  Procesadas {idx + 1} / {len(df)} formas para PDF...")

print("✅ PDF guardado en 'graficas_filtradas/todas_formas.pdf'")

# ============================================================================
# 7. MOSTRAR UN RESUMEN EN PANTALLA
# ============================================================================
print("\n=== RESUMEN DE COMBINACIONES FILTRADAS ===")
print(df.describe())

# Mostrar las 5 formas con menor energía (más flexibles)
print("\n=== 5 FORMAS CON MENOR ENERGÍA (MÁS FLEXIBLES) ===")
df_menor_energia = df.nsmallest(5, 'Energia')
print(df_menor_energia[['A', 'B', 'C', 'Energia']].to_string(index=False))

# Mostrar las 5 formas con mayor energía (más rígidas)
print("\n=== 5 FORMAS CON MAYOR ENERGÍA (MÁS RÍGIDAS) ===")
df_mayor_energia = df.nlargest(5, 'Energia')
print(df_mayor_energia[['A', 'B', 'C', 'Energia']].to_string(index=False))

# ============================================================================
# 8. GRÁFICA RESUMEN (todas las formas superpuestas)
# ============================================================================
print("\nGenerando gráfica resumen (formas superpuestas)...")

plt.figure(figsize=(10, 8))
for idx, row in df.iterrows():
    A = row['A']
    B = row['B']
    C = row['C']
    
    b_limite = np.sqrt(A**2 + C**2)
    x = np.linspace(-b_limite, b_limite, NUM_PUNTOS)
    y = Y_cassini(x, A, B, C)
    
    # Colorear por energía
    plt.plot(x, y, 'b-', alpha=0.2, linewidth=0.5)

plt.xlabel('x')
plt.ylabel('y')
plt.title('Superposición de todas las formas filtradas')
plt.axhline(y=0, color='black', linewidth=0.5)
plt.axvline(x=0, color='black', linewidth=0.5)
plt.grid(True, alpha=0.3)
plt.axis('equal')
plt.savefig('graficas_filtradas/superposicion.png', dpi=150, bbox_inches='tight')
plt.close()

print("✅ Gráfica de superposición guardada")

# ============================================================================
# 9. GRÁFICA: ENERGÍA vs. VOLUMEN (información adicional)
# ============================================================================
print("\nGenerando gráfica: Energía vs. Volumen...")

plt.figure(figsize=(10, 7))
scatter = plt.scatter(df['Volumen'], df['Energia'], 
                      c=df['B'], cmap='plasma', 
                      s=30, alpha=0.7, edgecolors='black', linewidth=0.5)

cbar = plt.colorbar(scatter)
cbar.set_label('Parámetro B', fontsize=12)

plt.xlabel('Volumen (μm³)', fontsize=14)
plt.ylabel('Energía de doblamiento (E)', fontsize=14)
plt.title('Energía vs. Volumen para formas filtradas', fontsize=14)
plt.grid(True, alpha=0.3)
plt.savefig('graficas_filtradas/energia_vs_volumen.png', dpi=150, bbox_inches='tight')
plt.close()

print("✅ Gráfica 'Energía vs. Volumen' guardada")

print("\n=== PROCESO COMPLETADO ===")
print("\nArchivos generados:")
print("  - filtro_con_energia.csv (datos con energía)")
print("  - graficas_filtradas/energia_vs_B.png")
print("  - graficas_filtradas/energia_vs_B_por_A.png")
print("  - graficas_filtradas/energia_vs_volumen.png")
print("  - graficas_filtradas/todas_formas.pdf")
print("  - graficas_filtradas/superposicion.png")

