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
def bending_energy(A, B, C, N=300):
    """
    Calcula la energía de curvatura para una GUV
    E = 2π ∫ (1/R₁ + 1/R₂)² · Y · √(1 + Y'²) dx
    """
    x_max = np.sqrt(A**2 + C**2)
    if x_max <= 1e-10:
        return 0.0
    
    x = np.linspace(0, x_max, N+1)
    dx = x[1] - x[0]
    
    y = np.zeros(N+1)
    yp = np.zeros(N+1)
    
    # Calcular Y(x) para cada punto
    for i, xi in enumerate(x):
        y[i] = Y_cassini(xi, A, B, C)
    
    # Derivada numérica (diferencias centrales)
    for i in range(1, N):
        if y[i] > 0:
            yp[i] = (Y_cassini(x[i]+dx, A, B, C) - Y_cassini(x[i]-dx, A, B, C)) / (2*dx)
    
    yp[0] = (y[1] - y[0]) / dx
    yp[N] = (y[N] - y[N-1]) / dx
    
    integrando = np.zeros(N+1)
    
    for i in range(1, N):
        if y[i] <= 1e-10:
            continue
        
        # Segunda derivada (para R1)
        ypp = (yp[i+1] - yp[i-1]) / (2*dx)
        
        # Radio R1 (curvatura en el plano del meridiano)
        if abs(ypp) > 1e-10:
            R1 = (1 + yp[i]**2)**1.5 / abs(ypp)
        else:
            R1 = 1e10
        
        # Radio R2 (curvatura perpendicular)
        if abs(yp[i]) > 1e-10:
            R2 = y[i] * np.sqrt(1 + yp[i]**2) / abs(yp[i])
        else:
            R2 = y[i]
        
        # Curvatura media y término del integrando
        curv_media = 1.0/R1 + 1.0/R2
        integrando[i] = (curv_media**2) * y[i] * np.sqrt(1 + yp[i]**2)
    
    # Integración por Simpson
    suma_impares = sum(integrando[1:N:2])
    suma_pares = sum(integrando[2:N-1:2])
    
    E = 2.0 * np.pi * (dx/3.0) * (integrando[0] + integrando[N] + 4*suma_impares + 2*suma_pares)
    
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
# 4. ORDENAR EL DATAFRAME POR B (DE MENOR A MAYOR)
# ============================================================================
df_ordenado = df.sort_values('B').reset_index(drop=True)
print(f"\n✅ DataFrame ordenado por parámetro B (de menor a mayor)")
print(f"   B mínimo: {df_ordenado['B'].min():.4f}")
print(f"   B máximo: {df_ordenado['B'].max():.4f}")

# ============================================================================
# 5. PARÁMETROS DE LA GRÁFICA
# ============================================================================
NUM_PUNTOS = 500  # Resolución de la curva

# ============================================================================
# 6. GENERAR PDF CON TODAS LAS GRÁFICAS (ORDENADAS POR B)
# ============================================================================
print("\nGenerando PDF con todas las gráficas (ordenadas por B de menor a mayor)...")

with PdfPages('graficas_filtradas/todas_formas_ordenadas_por_B.pdf') as pdf:
    for idx, row in df_ordenado.iterrows():
        A = row['A']
        B = row['B']
        C = row['C']
        E = row['Energia']
        
        b_limite = np.sqrt(A**2 + C**2)
        x = np.linspace(-b_limite, b_limite, NUM_PUNTOS)
        y = Y_cassini(x, A, B, C)
        
        # Crear figura con tamaño adecuado
        plt.figure(figsize=(8, 6))
        plt.plot(x, y, 'b-', linewidth=2, label='Perfil superior')
        plt.plot(x, -y, 'r-', linewidth=2, label='Perfil inferior')
        plt.fill_between(x, -y, y, color='lightblue', alpha=0.5)
        
        # Configurar gráfica
        plt.xlabel('x (eje horizontal)', fontsize=12)
        plt.ylabel('y (altura)', fontsize=12)
        plt.title(f'Forma #{idx+1}:  A={A:.4f}, B={B:.4f}, C={C:.4f}\nEnergía = {E:.4e}', fontsize=12)
        plt.axhline(y=0, color='black', linewidth=0.5)
        plt.axvline(x=0, color='black', linewidth=0.5)
        plt.grid(True, alpha=0.3)
        plt.axis('equal')
        plt.legend(loc='upper right', fontsize=8)
        
        # Agregar anotación con el orden y el valor de B
        plt.annotate(f'B = {B:.4f}', xy=(0.02, 0.95), xycoords='axes fraction', 
                     fontsize=10, bbox=dict(boxstyle="round,pad=0.3", facecolor='yellow', alpha=0.7))
        
        pdf.savefig()
        plt.close()
        
        if (idx + 1) % 20 == 0:
            print(f"  Procesadas {idx + 1} / {len(df_ordenado)} formas para PDF...")

print("✅ PDF guardado en 'graficas_filtradas/todas_formas_ordenadas_por_B.pdf'")

# ============================================================================
# 7. MOSTRAR RESUMEN DE VALORES DE B EN EL PDF
# ============================================================================
print("\n=== DISTRIBUCIÓN DEL PARÁMETRO B EN EL PDF ===")
print(f"Número total de formas: {len(df_ordenado)}")
print(f"B - Mínimo: {df_ordenado['B'].min():.6f}")
print(f"B - Máximo: {df_ordenado['B'].max():.6f}")
print(f"B - Media: {df_ordenado['B'].mean():.6f}")
print(f"B - Desviación estándar: {df_ordenado['B'].std():.6f}")

# Mostrar cuántas formas hay en cada rango de B
print("\n=== FORMAS POR RANGO DE B ===")
bins = [1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
for i in range(len(bins)-1):
    count = df_ordenado[(df_ordenado['B'] >= bins[i]) & (df_ordenado['B'] < bins[i+1])].shape[0]
    print(f"B ∈ [{bins[i]:.1f}, {bins[i+1]:.1f}): {count} formas")

# ============================================================================
# 8. GRÁFICA: ENERGÍA vs. B (con los datos ordenados)
# ============================================================================
print("\nGenerando gráfica: Energía de doblamiento vs. Parámetro B...")

plt.figure(figsize=(10, 7))

# Crear scatter plot coloreado por A
scatter = plt.scatter(df_ordenado['B'], df_ordenado['Energia'], 
                      c=df_ordenado['A'], cmap='viridis', 
                      s=40, alpha=0.7, edgecolors='black', linewidth=0.5)

# Barra de color para A
cbar = plt.colorbar(scatter)
cbar.set_label('Parámetro A', fontsize=12)

plt.xlabel('Parámetro B (altura/escala vertical)', fontsize=14)
plt.ylabel('Energía de doblamiento (E)', fontsize=14)
plt.title('Energía de curvatura vs. Parámetro B\n(ordenado por B de menor a mayor)', fontsize=14)
plt.grid(True, alpha=0.3)

# Guardar gráfica
plt.savefig('graficas_filtradas/energia_vs_B_ordenado.png', dpi=150, bbox_inches='tight')
plt.close()

print("✅ Gráfica 'Energía vs. B' guardada")

# ============================================================================
# 9. GRÁFICA: ENERGÍA vs. B (con líneas por valor de A)
# ============================================================================
print("\nGenerando gráfica: Energía vs. B (con líneas por valor de A)...")

plt.figure(figsize=(12, 8))

# Agrupar por A y graficar líneas
valores_A = sorted(df_ordenado['A'].unique())
colores = plt.cm.plasma(np.linspace(0, 1, len(valores_A)))

for i, A_val in enumerate(valores_A):
    subset = df_ordenado[df_ordenado['A'] == A_val]
    if len(subset) > 1:
        # Ya está ordenado por B
        plt.plot(subset['B'], subset['Energia'], 
                 'o-', color=colores[i], linewidth=1.5, markersize=4, 
                 label=f'A = {A_val:.2f}', alpha=0.8)

plt.xlabel('Parámetro B', fontsize=14)
plt.ylabel('Energía de doblamiento (E)', fontsize=14)
plt.title('Energía de curvatura vs. Parámetro B (agrupado por A)\n(ordenado por B)', fontsize=14)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
plt.grid(True, alpha=0.3)
plt.tight_layout()

plt.savefig('graficas_filtradas/energia_vs_B_por_A_ordenado.png', dpi=150, bbox_inches='tight')
plt.close()

print("✅ Gráfica 'Energía vs. B (por A)' guardada")

# ============================================================================
# 10. MOSTRAR LAS 5 FORMAS CON MENOR Y MAYOR ENERGÍA
# ============================================================================
print("\n=== 5 FORMAS CON MENOR ENERGÍA (MÁS FLEXIBLES) ===")
df_menor_energia = df_ordenado.nsmallest(5, 'Energia')
print(df_menor_energia[['A', 'B', 'C', 'Energia']].to_string(index=False))

print("\n=== 5 FORMAS CON MAYOR ENERGÍA (MÁS RÍGIDAS) ===")
df_mayor_energia = df_ordenado.nlargest(5, 'Energia')
print(df_mayor_energia[['A', 'B', 'C', 'Energia']].to_string(index=False))

print("\n=== PROCESO COMPLETADO ===")
print("\nArchivos generados:")
print("  - filtro_con_energia.csv (datos con energía)")
print("  - graficas_filtradas/todas_formas_ordenadas_por_B.pdf (PDF ordenado por B)")
print("  - graficas_filtradas/energia_vs_B_ordenado.png")
print("  - graficas_filtradas/energia_vs_B_por_A_ordenado.png")