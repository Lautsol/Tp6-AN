import numpy as np
import matplotlib.pyplot as plt

# Parámetros físicos
L = 10e-3  # Largo del dominio (m)
W = 1e-3   # Ancho del dominio (m)
Nx = 100   # Número de nodos en x
Ny = 10    # Número de nodos en y
dx = L / Nx  # Espaciado en x
dy = W / Ny  # Espaciado en y

alpha = 1.2e-7  # Difusividad térmica (m^2/s)
L_fusion = 334000  # Calor latente de fusión (J/kg)
rho = 1000  # Densidad del agua/hielo (kg/m^3)
cp = 2100   # Capacidad calorífica específica (J/(kg K))

# Temperaturas iniciales y de frontera
T_inicial = -10  # Temperatura inicial del hielo (°C)
T_borde_der_inicial = -10  # Temperatura inicial en el borde derecho (°C)
T_borde_der_final = 85     # Temperatura final en el borde derecho (°C)

# Inicialización de matrices
T = np.ones((Nx, Ny)) * T_inicial  # Campo de temperatura inicial
calor_latente_restante = np.zeros((Nx, Ny))  # Calor restante para cambio de fase
calor_latente_restante[:] = rho * L_fusion

# Criterio de estabilidad
dt_max = 0.5 * dx**2 / alpha

# Función para aplicar condiciones de frontera
def aplicar_condiciones_frontera(T, t):
    T[0, :] = T[1, :]  # Aislado en el borde izquierdo
    T[:, 0] = T[:, 1]  # Aislado en el borde superior
    T[:, -1] = T[:, -2]  # Aislado en el borde inferior
    T[-1, :] = T_borde_der_inicial + (T_borde_der_final - T_borde_der_inicial) * min(t / 10, 1)

# Función para realizar la simulación
def simular(dt):
    Nt = int(50 / dt)  # Número de pasos de tiempo
    T = np.ones((Nx, Ny)) * T_inicial
    calor_latente_restante = np.zeros((Nx, Ny))
    calor_latente_restante[:] = rho * L_fusion

    temperaturas_centrales = []
    tiempo = []
    
    for n in range(Nt):
        t = n * dt
        tiempo.append(t)
        aplicar_condiciones_frontera(T, t)
        T_new = T.copy()

        for i in range(1, Nx - 1):
            for j in range(1, Ny - 1):
                d2Tdx2 = (T[i + 1, j] - 2 * T[i, j] + T[i - 1, j]) / dx**2
                d2Tdy2 = (T[i, j + 1] - 2 * T[i, j] + T[i, j - 1]) / dy**2

                if T[i, j] >= 0 and calor_latente_restante[i, j] > 0:
                    calor_absorbido = rho * cp * alpha * dt * (d2Tdx2 + d2Tdy2)
                    calor_latente_restante[i, j] -= calor_absorbido
                    if calor_latente_restante[i, j] > 0:
                        T_new[i, j] = 0
                    else:
                        exceso_calor = -calor_latente_restante[i, j]
                        T_new[i, j] = 0 + exceso_calor / (rho * cp)
                else:
                    T_new[i, j] += alpha * dt * (d2Tdx2 + d2Tdy2)

        T = T_new
        temperaturas_centrales.append(T[:, Ny // 2].copy())
    return temperaturas_centrales, tiempo, T

# Valores de dt para analizar
dt_values = [dt_max / 2, dt_max, dt_max * 1.5]
plt.figure(figsize=(12, 8))

# Graficar efecto del paso de tiempo en la estabilidad
for dt in dt_values:
    temperaturas_centrales, _, _ = simular(dt)
    plt.plot(np.linspace(0, L, Nx), temperaturas_centrales[-1], label=f"dt = {dt:.1e}")

plt.title("Efecto del paso de tiempo en la estabilidad del sistema")
plt.xlabel("Posición (m)")
plt.ylabel("Temperatura (°C)")
plt.legend()
plt.grid()
plt.show()

# Simulación con un dt adecuado (dt_max / 2) para analizar evolución temporal
dt_adecuado = dt_max / 2
temperaturas_centrales, tiempo, T_final = simular(dt_adecuado)

# Graficar la evolución de la temperatura en la línea central
plt.figure(figsize=(10, 6))
for i in range(0, len(temperaturas_centrales), len(temperaturas_centrales) // 5):
    plt.plot(np.linspace(0, L, Nx), temperaturas_centrales[i], label=f"t = {tiempo[i]:.1f} s")

plt.title("Evolución de la temperatura en la línea central")
plt.xlabel("Posición (m)")
plt.ylabel("Temperatura (°C)")
plt.legend()
plt.grid()
plt.show()
