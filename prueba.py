import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# ============================================================================
# 1. FUNCIÓN DIRECTA Y(x)
# ============================================================================
def Y_cassini(x, A, B, C):
    interior = np.sqrt(C**4 + 4.0 * A**2.0 * x**2.0) - A**2.0 - x**2.0
    with np.errstate(invalid='ignore'):
        y = B * np.sqrt(np.maximum(interior, 0))
    return y


# ============================================================================
# 2. FUNCIÓN INVERSA: x(y) resolviendo Y(x) = y
# ============================================================================
def X_cassini(y, A, B, C):
    """
    Calcula x dado y resolviendo Y_cassini(x) = y
    Usa brentq (búsqueda de raíz) porque la función es monótona
    """
    if y <= 0:
        return 0.0
    
    x_max = np.sqrt(A**2 + C**2)
    if x_max <= 0:
        return 0.0
    
    # Valor máximo de Y (en x = x_max/2 aproximadamente)
    y_max = Y_cassini(x_max/2, A, B, C)
    
    if y >= y_max:
        return x_max / 2
    
    # Buscar raíz de Y(x) - y = 0
    def f(x):
        return Y_cassini(x, A, B, C) - y
    
    try:
        # La función es monótona decreciente en [0, x_max]
        x_root = brentq(f, 0, x_max)
        return x_root
    except (ValueError, RuntimeError):
        return 0.0


# ============================================================================
# 3. DERIVADA DE x RESPECTO A y (para la integración en y)
# ============================================================================
def dX_dy(y, A, B, C, delta=1e-6):
    """
    Derivada numérica de x respecto a y
    dx/dy ≈ (x(y+δ) - x(y-δ)) / (2δ)
    """
    if y <= 0:
        return 0.0
    
    x_plus = X_cassini(y + delta, A, B, C)
    x_minus = X_cassini(max(0, y - delta), A, B, C)
    
    return (x_plus - x_minus) / (2 * delta)


# ============================================================================
# 4. SEGUNDA DERIVADA DE x RESPECTO A y (opcional, para mejor precisión)
# ============================================================================
def d2X_dy2(y, A, B, C, delta=1e-6):
    """
    Segunda derivada numérica de x respecto a y
    d²x/dy² ≈ (x(y+δ) - 2x(y) + x(y-δ)) / δ²
    """
    if y <= 0:
        return 0.0
    
    x_plus = X_cassini(y + delta, A, B, C)
    x = X_cassini(y, A, B, C)
    x_minus = X_cassini(max(0, y - delta), A, B, C)
    
    return (x_plus - 2*x + x_minus) / (delta**2)


# ============================================================================
# 5. ENERGÍA DE CURVATURA: MÉTODO DE CANHAM (integración en y en el borde)
# ============================================================================
def bending_energy_canham(A, B, C, N=200, umbral_pendiente=10.0):
    """
    Implementa el método de Canham:
    - Integración en x en la región central (donde |dy/dx| < umbral)
    - Integración en y cerca del borde (donde |dy/dx| → ∞)
    """
    x_max = np.sqrt(A**2 + C**2)
    if x_max <= 1e-10:
        return 0.0
    
    # ========================================================================
    # PASO 1: Encontrar el punto de transición x_trans
    # donde |dy/dx| supera el umbral
    # ========================================================================
    x_test = np.linspace(0, x_max, 200)
    y_test = Y_cassini(x_test, A, B, C)
    dy_dx_test = np.gradient(y_test, x_test)
    
    x_trans = x_max
    for i in range(len(x_test)):
        if abs(dy_dx_test[i]) > umbral_pendiente:
            x_trans = x_test[i]
            break
    
    y_trans = Y_cassini(x_trans, A, B, C)
    
    # ========================================================================
    # PASO 2: Integración en x (desde 0 hasta x_trans)
    # ========================================================================
    N_x = N // 2
    x_central = np.linspace(0, x_trans, N_x + 1)
    dx = x_central[1] - x_central[0]
    
    E_x = 0.0
    for i in range(1, N_x):
        xi = x_central[i]
        yi = Y_cassini(xi, A, B, C)
        if yi <= 1e-10:
            continue
        
        # Derivadas numéricas
        yp = (Y_cassini(xi + dx, A, B, C) - Y_cassini(xi - dx, A, B, C)) / (2*dx)
        ypp = (Y_cassini(xi + dx, A, B, C) - 2*yi + Y_cassini(xi - dx, A, B, C)) / (dx**2)
        
        # Radios
        if abs(ypp) > 1e-10:
            R1 = (1 + yp**2)**1.5 / abs(ypp)
        else:
            R1 = 1e10
        
        if abs(yp) > 1e-10:
            R2 = yi * np.sqrt(1 + yp**2) / abs(yp)
        else:
            R2 = yi
        
        curv = 1.0/R1 + 1.0/R2
        integrando = (curv**2) * yi * np.sqrt(1 + yp**2)
        
        # Peso de Simpson
        if i == 1 or i == N_x:
            E_x += integrando
        elif i % 2 == 0:
            E_x += 2 * integrando
        else:
            E_x += 4 * integrando
    
    E_x = 2.0 * np.pi * (dx/3.0) * E_x
    
    # ========================================================================
    # PASO 3: Integración en y (desde y_trans hasta 0)
    # ========================================================================
    N_y = N // 2
    y_borde = np.linspace(0, y_trans, N_y + 1)
    dy = y_borde[1] - y_borde[0]
    
    E_y = 0.0
    for i in range(1, N_y):
        yi = y_borde[i]
        if yi <= 0:
            continue
        
        # Calcular x(y) y sus derivadas
        xi = X_cassini(yi, A, B, C)
        xp = dX_dy(yi, A, B, C)      # dx/dy
        xpp = d2X_dy2(yi, A, B, C)   # d²x/dy²
        
        # Radios en coordenadas invertidas
        # ρ₁ = -(1 + x'²)^(3/2) / x''
        if abs(xpp) > 1e-10:
            rho1 = -(1 + xp**2)**1.5 / abs(xpp)
        else:
            rho1 = 1e10
        
        # ρ₂ = x · √(1 + x'²)
        rho2 = xi * np.sqrt(1 + xp**2)
        
        curv = 1.0/rho1 + 1.0/rho2
        integrando = (curv**2) * xi * np.sqrt(1 + xp**2)
        
        # Peso de Simpson
        if i == 1 or i == N_y:
            E_y += integrando
        elif i % 2 == 0:
            E_y += 2 * integrando
        else:
            E_y += 4 * integrando
    
    E_y = 2.0 * np.pi * (dy/3.0) * E_y
    
    # ========================================================================
    # ENERGÍA TOTAL
    # ========================================================================
    return E_x + E_y


# ============================================================================
# 6. VERSIÓN SIMPLIFICADA (solo integración en y, para comparar)
# ============================================================================
def bending_energy_en_y(A, B, C, N=300):
    """
    Versión simplificada: integración completa en y
    (útil para formas con borde muy afilado)
    """
    x_max = np.sqrt(A**2 + C**2)
    y_max = Y_cassini(x_max/2, A, B, C)
    
    if y_max <= 1e-10:
        return 0.0
    
    y = np.linspace(0, y_max, N+1)
    dy = y[1] - y[0]
    
    E = 0.0
    for i in range(1, N):
        yi = y[i]
        if yi <= 0:
            continue
        
        xi = X_cassini(yi, A, B, C)
        xp = dX_dy(yi, A, B, C)
        xpp = d2X_dy2(yi, A, B, C)
        
        if abs(xpp) > 1e-10:
            rho1 = -(1 + xp**2)**1.5 / abs(xpp)
        else:
            rho1 = 1e10
        
        rho2 = xi * np.sqrt(1 + xp**2)
        
        curv = 1.0/rho1 + 1.0/rho2
        integrando = (curv**2) * xi * np.sqrt(1 + xp**2)
        
        if i == 1 or i == N:
            E += integrando
        elif i % 2 == 0:
            E += 2 * integrando
        else:
            E += 4 * integrando
    
    E = 2.0 * np.pi * (dy/3.0) * E
    return E


# ============================================================================
# 7. PRUEBA RÁPIDA
# ============================================================================
if __name__ == "__main__":
    # Parámetros de prueba
    A = 1.0
    C = 1.05
    B_vals = np.arange(0.3, 1.1, 0.05)
    
    energias = []
    for B in B_vals:
        E = bending_energy_canham(A, B, C, N=200)
        energias.append(E)
        print(f"B={B:.2f}, E={E:.4f}")
    
    plt.figure(figsize=(8, 6))
    plt.plot(B_vals, energias, 'o-', linewidth=2)
    plt.xlabel('Parámetro B')
    plt.ylabel('Energía de curvatura')
    plt.title('Método de Canham: Energía vs B')
    plt.grid(True, alpha=0.3)
    plt.show()