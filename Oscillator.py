import numpy as np
import matplotlib.pyplot as plt
import math

print('Please, type in the values of w0, w and h respectively.')
print('w0 := ', end='')
w0 = float(input())
print('w := ', end='')
w = float(input())
print('h := ', end='')
h = float(input())

# sample input: 
#     w0 = 1.0
#     w = 1.01
#     h = 0.1

# The equation of the parametric oscillator stands as follows:
#     d^2theta/dt^2 + w0^2 * (1 + h * cos(2w * t)) * theta = 0

# It is suggested that we find a solution in the form of 
#     theta(t) = e^(mu * t) * cos(w * t + phi),
# which gives
#     theta'(t) = mu * e^(mu * t) * cos(w * t + phi) - e^(mu * t) * w * sin(w * t + phi)

# An approximation for Tailor series gives the following system of equations on phi:
#     (w0^2 - w^2 + mu^2 + h/2 * w0^2) * cos(phi) - 2 * w * mu * sin(phi) = 0;
#     2 w * mu * cos(phi) + (w0^2 - w^2 + mu^2 - h/2 * w0^2) * sin(phi) = 0

# The system above has a solution iff its determinant equals 0:
#    mu^4 +  2(w^2 + w0)^2 * mu^2 + (w0^2 - w^2)^2 - h^2/4 * w0^4 = 0

# If we have h <= 2 * |1 - w^2/w0^2|, then only complex solutions exist; however, we are interested in real solutions, and thus demand 
#    h > 2 * |1 - w^2/w0^2|
# Now we check, whether the above-listed condition stands, and if it does, find two possible values of mu:
# print(2 * abs(1 - (w / w0) ** 2))

if h <= 2 * abs(1 - (w / w0) ** 2):
    print('Unfortunately, this parametric oscillator has no real solutions')
else:
    b = 2 * (w * w + w0 * w0)
    c = (w * w - w0 * w0) ** 2 - h * h / 4 * w0 ** 4
    D = b * b - 4 * c
    mu = math.sqrt((-b + math.sqrt(D)) / 2)
    
    def P(mu):
        return w0 * w0 - w * w + mu * mu + h / 2 * w0 * w0
    
    def Q(mu):
        return 2 * w * mu
    
    phi = math.atan2(P(mu), Q(mu))
    
# Now we draw the two solutions, for the positive and negative values of mu respectively:
    
    def theta(t, phi, mu):
        return np.exp(mu * t) * np.cos(w * t + phi)
    
    def theta_o(t, phi, mu):
        return mu * np.exp(mu * t) * np.cos(w * t + phi) - np.exp(mu * t) * w * np.sin(w * t + phi)
    
    plt.figure(1)
    plt.subplot(121)
    plt.axis([-100, 100, -100, 100])
    t = np.arange(0.0, 200.0, 0.005)
    plt.xlabel('theta')
    plt.ylabel("theta'")
    plt.title('Solution for mu > 0')
    plt.grid(True)    
    plt.plot(theta(t, phi, mu), theta_o(t, phi, mu))
    
    plt.subplot(122)
    plt.axis([-0.1, 0.1, -0.1, 0.1])
    t = np.arange(0.0, 200.0, 0.005)
    plt.xlabel('theta')
    plt.ylabel("theta'")
    plt.title('Solution for mu < 0')
    plt.grid(True)    
    plt.plot(theta(t, phi, -mu), theta_o(t, phi, -mu))
    
    # The bifurcation diagram can be represented as a solution of the following equality:
    #     h = 2 * |1 - (w/w0)^2|
    # on w -- h plane with w and h being the parameters and w0 fixed
    
    plt.figure(2)
    plt.axis([0, 2 * w0, 0, 2])
    t = np.arange(0, 2 * w0, 0.005)
    plt.xlabel('w')
    plt.ylabel("h")
    plt.title('Bifurcation diagram for parametric oscillator with w0 = ' + str(w0))
    plt.grid(True)  
    plt.plot(t, 2 * abs(1 - (t / w0) ** 2))
    
    plt.show()