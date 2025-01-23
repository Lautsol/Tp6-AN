import numpy as np
import matplotlib.pyplot as plt

# Parámetros físicos
L = 10e-3  # Largo del dominio (m)
W = 1e-3   # Ancho del dominio (m)
Nx = 100   # Número de nodos en x
Ny = 10    # Número de nodos en y
dx = L / Nx  # Espaciado en x
dy = W / Ny  # Espaciado en y

dt = 0.01  # Paso de tiempo (s)
T_total = 100  # Tiempo total de simulación (s)
Nt = int(T_total / dt)  # Número de pasos de tiempo

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

# Inicializar calor latente restante solo en nodos que comiencen con hielo
calor_latente_restante[:] = rho * L_fusion

# Función para aplicar condiciones de frontera
def aplicar_condiciones_frontera(T, t):
    # Bordes izquierdo, superior e inferior: aislados térmicamente
    T[0, :] = T[1, :]
    T[:, 0] = T[:, 1]
    T[:, -1] = T[:, -2]

    # Borde derecho: rampa lineal de temperatura
    T[-1, :] = T_borde_der_inicial + (T_borde_der_final - T_borde_der_inicial) * min(t / 10, 1)

# Simulación
tiempo = []
temperatura_central = []

for n in range(Nt):
    t = n * dt

    # Almacenar la temperatura en la línea central
    tiempo.append(t)
    temperatura_central.append(T[:, Ny // 2].copy())

    # Aplicar condiciones de frontera
    aplicar_condiciones_frontera(T, t)

    # Calcular el nuevo campo de temperatura
    T_new = T.copy()
    for i in range(1, Nx - 1):
        for j in range(1, Ny - 1):
            d2Tdx2 = (T[i + 1, j] - 2 * T[i, j] + T[i - 1, j]) / dx**2
            d2Tdy2 = (T[i, j + 1] - 2 * T[i, j] + T[i, j - 1]) / dy**2

            # Si estamos en el cambio de fase
            if T[i, j] >= 0 and calor_latente_restante[i, j] > 0:
                # Calor absorbido en este paso
                calor_absorbido = rho * cp * alpha * dt * (d2Tdx2 + d2Tdy2)
                calor_latente_restante[i, j] -= calor_absorbido

                # Mientras quede calor latente, la temperatura no sube
                if calor_latente_restante[i, j] > 0:
                    T_new[i, j] = 0
                else:
                    # El cambio de fase termina, continuar con el aumento de temperatura
                    exceso_calor = -calor_latente_restante[i, j]  # Calor que excedió el necesario para el cambio de fase
                    T_new[i, j] = 0 + exceso_calor / (rho * cp)
            else:
                # Actualización normal fuera del punto de cambio de fase
                T_new[i, j] += alpha * dt * (d2Tdx2 + d2Tdy2)

    # Actualizar la matriz de temperatura
    T = T_new

# Graficar la temperatura en la línea central en diferentes tiempos
plt.figure(figsize=(10, 6))
for i in range(0, Nt, Nt // 5):
    plt.plot(np.linspace(0, L, Nx), temperatura_central[i], label=f"t = {i * dt:.1f} s")

plt.title("Evolución de la temperatura en la línea central")
plt.xlabel("Posición (m)")
plt.ylabel("Temperatura (°C)")
plt.legend()
plt.grid()
plt.show()
