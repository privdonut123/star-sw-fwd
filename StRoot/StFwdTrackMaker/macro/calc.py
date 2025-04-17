import numpy as np
import matplotlib.pyplot as plt
deltaZ = 30 #cm



def beta_velocity( p, m =139.57 ):
    """
    Calculate the velocity factor beta = v/c
    p: momentum in MeV/c
    m: mass in MeV/c^2
    """
    return np.sqrt(1 - (m**2 / (p**2 + m**2)))

def multiple_scattering_rms_angle(p, x=0.1, X0=10, m = 139.57, z = 1.0):
    """
    Calculate the multiple scattering angle in radians.
    x: distance in cm traversed in the material
    z: charge of the particle
    p: momentum in MeV/c
    X0: radiation length in cm
    beta: velocity factor (v/c)
    13.6 is the constant for the multiple scattering formula in MeV
    """
    beta = beta_velocity(p, m)
    return (13.6 / (beta * p)) * z * np.sqrt(x / X0) * (1 + 0.038 * np.log(x / X0))

def Lperp( eta):
    theta = 2 * np.arctan(np.exp(-eta))
    S = deltaZ / np.cos(theta) 
    L = S * np.sin(theta)
    return L

def curvature( pt ):
    # pt in MeV/c
    # curvature in 1/cm
    return (300.0 / pt) * 1e-2

print( "L_perp @eta=2.5 == %0.3f" % Lperp( 2.5 ))
print( "L_perp @eta=4.0 == %0.3f" % (Lperp( 4.0 )))

print( "curvature( 0.1 ) == %0.3f [1/cm]" % (curvature( 0.1 )))
print( "curvature( 0.5 ) == %0.3f" % (curvature( 0.5 )))
print( "curvature( 1.0 ) == %0.3f" % (curvature( 1.0 )))

print( "delta phi (eta=2.5, pT=0.1) == %0.3f [rad]" % (Lperp( 2.5 ) * curvature( 0.1 )))
print( "delta phi (eta=2.5, pT=0.5) == %0.3f" % (Lperp( 2.5 ) * curvature( 0.5 )))
print( "delta phi (eta=2.5, pT=1.0) == %0.3f" % (Lperp( 2.5 ) * curvature( 1.0 )))

print( "delta phi (eta=3.5, pT=0.1) == %0.3f" % (Lperp( 3.5 ) * curvature( 0.1 )))
print( "delta phi (eta=3.5, pT=0.5) == %0.3f" % (Lperp( 3.5 ) * curvature( 0.5 )))
print( "delta phi (eta=3.5, pT=1.0) == %0.3f" % (Lperp( 3.5 ) * curvature( 1.0 )))

print( "delta phi (eta=4.0, pT=0.1) == %0.3f" % (Lperp( 4.0 ) * curvature( 0.1 )))
print( "delta phi (eta=4.0, pT=0.5) == %0.3f" % (Lperp( 4.0 ) * curvature( 0.5 )))
print( "delta phi (eta=4.0, pT=1.0) == %0.3f" % (Lperp( 4.0 ) * curvature( 1.0 )))

sigma_FST = 0.0040 #0.000835 #radians from 3 FST planes
print( "Expected momentum resolution given sigma_total = %0.3f (single disk)" % (sigma_FST))
print( "DeltaP/P (eta=2.5, pT=0.1) = %0.3f" % (sigma_FST / (Lperp( 2.5 ) * curvature( 0.1 ))))
print( "DeltaP/P (eta=2.5, pT=0.5) = %0.3f" % (sigma_FST / (Lperp( 2.5 ) * curvature( 0.5 ))))
print( "DeltaP/P (eta=2.5, pT=1.0) = %0.3f" % (sigma_FST / (Lperp( 2.5 ) * curvature( 1.0 ))))


print( "DeltaP/P (eta=4.0, pT=0.1) = %0.3f" % (sigma_FST / (Lperp( 4.0 ) * curvature( 0.1 ))))
print( "DeltaP/P (eta=4.0, pT=0.5) = %0.3f" % (sigma_FST / (Lperp( 4.0 ) * curvature( 0.5 ))))
print( "DeltaP/P (eta=4.0, pT=1.0) = %0.3f" % (sigma_FST / (Lperp( 4.0 ) * curvature( 1.0 ))))

sigma_total = 0.000835 #radians from 3 FST planes
print( "Expected momentum resolution given sigma_total = %0.3f (ideal 3 disks)" % (sigma_total))
print( "DeltaP/P (eta=2.5, pT=0.1) = %0.3f" % (sigma_total / (Lperp( 2.5 ) * curvature( 0.1 ))))
print( "DeltaP/P (eta=2.5, pT=0.5) = %0.3f" % (sigma_total / (Lperp( 2.5 ) * curvature( 0.5 ))))
print( "DeltaP/P (eta=2.5, pT=1.0) = %0.3f" % (sigma_total / (Lperp( 2.5 ) * curvature( 1.0 ))))


print( "DeltaP/P (eta=4.0, pT=0.1) = %0.3f" % (sigma_total / (Lperp( 4.0 ) * curvature( 0.1 ))))
print( "DeltaP/P (eta=4.0, pT=0.5) = %0.3f" % (sigma_total / (Lperp( 4.0 ) * curvature( 0.5 ))))
print( "DeltaP/P (eta=4.0, pT=1.0) = %0.3f" % (sigma_total / (Lperp( 4.0 ) * curvature( 1.0 ))))

eta = np.linspace( 2.5, 4.0, 100 )
pT = np.linspace( 10, 1000, 100 )
reso_eta_pT0p1 = sigma_FST / (Lperp( eta ) * curvature( 100 ))

plt.plot( eta, reso_eta_pT0p1, label='pT=0.1 GeV/c', color='blue' )
plt.show()

reso_pT_eta2p5 = sigma_FST / (Lperp( 3.0 ) * curvature( pT ))
plt.plot( pT, reso_pT_eta2p5, label='angular resolution', color='red' )
plt.xlabel('pT [GeV/c]')
plt.ylabel('DeltaP/P')
plt.title('Momentum Resolution vs pT at eta=2.5')

plt.grid()


multiple_scat = multiple_scattering_rms_angle( pT )
plt.plot( pT, multiple_scat, label='Multiple Scattering', color='green' )
plt.plot( pT, reso_pT_eta2p5 + multiple_scat, label='Total Resolution', color='orange' )
plt.legend()

# set the y-axis limit
plt.ylim(0, 1.2)

plt.show()


def combined_resolution( nFST=3, nSTGC=1, sigma_STGC = 0.001 ):
    """
    Calculate the combined resolution from the number of FST and STGC layers.
    nFST: number of FST layers
    nSTGC: number of STGC layers
    """
    return np.sqrt( (sigma_FST / np.sqrt(nFST-1))**2 + (sigma_STGC / np.sqrt(nSTGC))**2 )

print( "combined resolution (nFST=3, nSTGC=1) == %0.5f" % (combined_resolution( 3, 1 )))
print( "combined resolution (nFST=3, nSTGC=2) == %0.5f" % (combined_resolution( 3, 2 )))
print( "combined resolution (nFST=3, nSTGC=3) == %0.5f" % (combined_resolution( 3, 3 )))
print( "combined resolution (nFST=3, nSTGC=4) == %0.5f" % (combined_resolution( 3, 4 )))