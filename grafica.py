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
# 4. PARÁMETROS DE LA GRÁFICA
# ============================================================================
NUM_PUNTOS = 500  # Resolución de la curva


# ============================================================================
# 6. OPCIONAL: GENERAR UN PDF CON TODAS LAS GRÁFICAS
# ============================================================================
print("\nGenerando PDF con todas las gráficas...")

with PdfPages('graficas_filtradas/todas_formas.pdf') as pdf:
    for idx, row in df.iterrows():
        A = row['A']
        B = row['B']
        C = row['C']
        
        b_limite = np.sqrt(A**2 + C**2)
        x = np.linspace(-b_limite, b_limite, NUM_PUNTOS)
        y = Y_cassini(x, A, B, C)
        
        plt.figure(figsize=(8, 6))
        plt.plot(x, y, 'b-', linewidth=2)
        plt.plot(x, -y, 'r-', linewidth=2)
        plt.fill_between(x, -y, y, color='lightblue', alpha=0.5)
        
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title(f'A={A:.4f}, B={B:.4f}, C={C:.4f}')
        plt.axhline(y=0, color='black', linewidth=0.5)
        plt.axvline(x=0, color='black', linewidth=0.5)
        plt.grid(True, alpha=0.3)
        plt.axis('equal')
        
        pdf.savefig()
        plt.close()

print("✅ PDF guardado en 'graficas_filtradas/todas_formas.pdf'")

# ============================================================================
# 7. OPCIONAL: MOSTRAR UN RESUMEN EN PANTALLA
# ============================================================================
print("\n=== RESUMEN DE COMBINACIONES FILTRADAS ===")
print(df.describe())

# ============================================================================
# 8. OPCIONAL: GRÁFICA RESUMEN (todas las formas superpuestas)
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
    
    plt.plot(x, y, 'b-', alpha=0.3, linewidth=0.5)

plt.xlabel('x')
plt.ylabel('y')
plt.title('Superposición de todas las formas filtradas')
plt.axhline(y=0, color='black', linewidth=0.5)
plt.axvline(x=0, color='black', linewidth=0.5)
plt.grid(True, alpha=0.3)
plt.axis('equal')
plt.savefig('graficas_filtradas/superposicion.png', dpi=150, bbox_inches='tight')
plt.close()

print("✅ Gráfica de superposición guardada en 'graficas_filtradas/superposicion.png'")

print("\n=== PROCESO COMPLETADO ===")