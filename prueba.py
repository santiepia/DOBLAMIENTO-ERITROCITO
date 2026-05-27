import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ============================================================================
# 1. LEER EL ARCHIVO CON LOS DATOS
# ============================================================================
print("Leyendo archivo filtro_con_energia.csv...")

try:
    # Intentar leer como CSV
    df = pd.read_csv('filtro_con_energia.csv')
    print("Archivo leído como CSV")
except:
    try:
        # Si no, leer como texto tabulado
        df = pd.read_csv('filtro_con_energia.txt', sep='\t')
        print("Archivo leído como TXT tabulado")
    except:
        # Si falla, leer formato whitespace
        df = pd.read_csv('filtro_con_energia.csv', delim_whitespace=True)
        print("Archivo leído con formato whitespace")

print(f"Total de combinaciones en el archivo: {len(df)}")
print(f"Columnas disponibles: {list(df.columns)}")

# ============================================================================
# 2. FILTRAR POR A = 4.6 (con cierta tolerancia)
# ============================================================================
A_objetivo = 4.6
tolerancia = 0.01

df_filtrado = df[(df['A'] >= A_objetivo - tolerancia) & 
                  (df['A'] <= A_objetivo + tolerancia)]

print(f"\nCombinaciones con A ≈ {A_objetivo}: {len(df_filtrado)}")

if len(df_filtrado) == 0:
    print(f"\n⚠️ No se encontraron datos para A = {A_objetivo}")
    print("Valores de A disponibles en el archivo:")
    print(sorted(df['A'].unique()))
    exit()

# ============================================================================
# 3. ORDENAR POR B PARA UNA LÍNEA SUAVE
# ============================================================================
df_ordenado = df_filtrado.sort_values('B')

# ============================================================================
# 4. CREAR LA GRÁFICA
# ============================================================================
plt.figure(figsize=(10, 7))

# Puntos individuales (scatter)
plt.scatter(df_ordenado['B'], df_ordenado['Energia'], 
            color='blue', s=60, alpha=0.7, edgecolors='black', 
            label='Datos calculados', zorder=3)

# Línea que conecta los puntos (tendencia)
plt.plot(df_ordenado['B'], df_ordenado['Energia'], 
         'b-', linewidth=2, alpha=0.8, label='Tendencia', zorder=2)

# Marcar el punto de mínima energía (más flexible)
idx_min = df_ordenado['Energia'].idxmin()
B_min = df_ordenado.loc[idx_min, 'B']
E_min = df_ordenado.loc[idx_min, 'Energia']

plt.scatter(B_min, E_min, color='red', s=150, marker='*', 
            edgecolors='black', linewidth=1.5,
            label=f'Mínimo: B={B_min:.3f}, E={E_min:.4f}', zorder=4)

# Configurar gráfica
plt.xlabel('Parámetro B (altura/escala vertical)', fontsize=14)
plt.ylabel('Energía de doblamiento (E)', fontsize=14)
plt.title(f'Energía de curvatura vs. Parámetro B\n(A = {A_objetivo:.2f})', fontsize=14)
plt.grid(True, alpha=0.3)
plt.legend(fontsize=11)

# Agregar anotación con el mínimo
plt.annotate(f'Mínimo en B = {B_min:.3f}', 
             xy=(B_min, E_min),
             xytext=(B_min + 0.1, E_min + 1),
             arrowprops=dict(arrowstyle='->', color='red', lw=1.5),
             fontsize=11, color='red')

# Ajustar límites para mejor visualización (opcional)
# plt.xlim(min(df_ordenado['B']) - 0.05, max(df_ordenado['B']) + 0.05)
# y_min = df_ordenado['Energia'].min()
# y_max = df_ordenado['Energia'].max()
# plt.ylim(y_min - 1, y_max + 1)

plt.tight_layout()

# Guardar la gráfica
plt.savefig('grafica_A_fijo_4.6.png', dpi=150, bbox_inches='tight')
plt.savefig('grafica_A_fijo_4.6.pdf', bbox_inches='tight')

# Mostrar la gráfica en pantalla
plt.show()

# ============================================================================
# 5. MOSTRAR INFORMACIÓN ADICIONAL
# ============================================================================
print("\n=== RESUMEN PARA A = 4.6 ===")
print(df_ordenado[['A', 'B', 'C', 'Energia']].to_string(index=False))

print(f"\n=== ESTADÍSTICAS ===")
print(f"B - Mínimo: {df_ordenado['B'].min():.4f}")
print(f"B - Máximo: {df_ordenado['B'].max():.4f}")
print(f"B - Media: {df_ordenado['B'].mean():.4f}")
print(f"Energía - Mínimo: {E_min:.6f} en B = {B_min:.4f}")
print(f"Energía - Máximo: {df_ordenado['Energia'].max():.6f}")

print("\n✅ Gráfica guardada como 'grafica_A_fijo_4.6.png'")