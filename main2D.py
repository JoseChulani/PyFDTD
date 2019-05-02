import numpy as np
import math
import scipy.constants
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ==== Preamble ===============================================================
c0   = scipy.constants.speed_of_light
mu0  = scipy.constants.mu_0
eps0 = scipy.constants.epsilon_0
imp0 = math.sqrt(mu0 / eps0)

# ==== Inputs / Pre-processing ================================================ 
# ---- Problem definition -----------------------------------------------------
L         = 10.0
dx        = 0.1
dy        = 0.1
finalTime = L/c0*5
cfl       = .99

gridEX = np.linspace(0,      L,        num=L/dx+1, endpoint=True)
gridEY = np.linspace(0,      L,        num=L/dy+1, endpoint=True)
gridHX = np.linspace(dx/2.0, L-dx/2.0, num=L/dx,   endpoint=True)
gridHY = np.linspace(dy/2.0, L-dy/2.0, num=L/dy,   endpoint=True)

# ---- Materials --------------------------------------------------------------
#PML
pmlxStart=3.0/4.0*L
pmlxBegin=1.0/4.0*L
pmlyStart=3.0/4.0*L
pmlyBegin=1.0/4.0*L
pmlSigmaE0X=1
pmlSigmaH0X=pmlSigmaE0X*mu0/eps0
pmlSigmaE0Y=0.1
pmlSigmaH0Y=pmlSigmaE0Y*mu0/eps0


# ---- Boundary conditions ----------------------------------------------------
 
# ---- Sources ----------------------------------------------------------------
# Initial field
spread   = 1.0
center   = (L/2.0, L/2.0)

initialH = np.zeros((gridHX.size, gridHY.size))
for i in range(gridHX.size):
    for j in range(gridHY.size):
        initialH[i,j] = math.exp( 
            - ((gridHX[i]-center[0])**2 + (gridHY[j]-center[1])**2) /
            math.sqrt(2.0) / spread)




#Para poder quitar los bucles for buscar en internet "double index slicing"
vector=[1,2,3,4,5]

# ---- Output requests --------------------------------------------------------
samplingPeriod = 0.0
 
# ==== Processing =============================================================
# ---- Solver initialization --------------------------------------------------
dt = cfl * dx / c0 / math.sqrt(2)
numberOfTimeSteps = int( finalTime / dt )

if samplingPeriod == 0.0:
    samplingPeriod = dt 
nSamples  = int( math.floor(finalTime/samplingPeriod) )
probeH    = np.zeros((gridHX.size, gridHY.size, nSamples))
probeTime = np.zeros(nSamples) 

exOld = np.zeros((gridEX.size, gridEY.size))
exNew = np.zeros((gridEX.size, gridEY.size))
eyOld = np.zeros((gridEX.size, gridEY.size))
eyNew = np.zeros((gridEX.size, gridEY.size))
hzOld = np.zeros((gridHX.size, gridHY.size))
hzNew = np.zeros((gridHX.size, gridHY.size))

pmlIndexX=np.searchsorted(gridEX,3.0/4.0*L)
pmlIndexX2=np.searchsorted(gridEX,L/4.0)
pmlIndexY=np.searchsorted(gridEY,3.0/4.0*L)


if 'initialH' in locals():
    hzOld = initialH

# Determines recursion coefficients
cEx = dt / eps0 / dx
cEy = dt / eps0 / dy
cHx = dt / mu0  / dx
cHy = dt / mu0  / dy

# ---- Time integration -------------------------------------------------------
print('--- Processing starts---')
tic = time.time();

t = 0.0
for n in range(numberOfTimeSteps):
    # --- Updates E field ---
    #for i in range(1, gridEX.size-1):
    #    for j in range(1, gridEY.size-1):
    #        exNew[i][j] = exOld[i][j] + cEy * (hzOld[i][j] - hzOld[i  ][j-1])
    #        eyNew[i][j] = eyOld[i][j] - cEx * (hzOld[i][j] - hzOld[i-1][j  ])
    
    #for i in range(1, gridEY.size-1):
    exNew[1:-1, 1:-1] = exOld[1:-1, 1:-1] + cEy * (hzOld[1:, 1:] - hzOld[1:, :-1])
    eyNew[1:-1, 1:-1] = eyOld[1:-1, 1:-1] - cEx * (hzOld[1:, 1:] - hzOld[:-1, 1: ])
        
    
    # E field boundary conditions
    # ---Empieza PML ---
    #for i in range(pmlIndexX,gridEX.size-1):
    #    pmlSigmaEX=pmlSigmaE0X*pow((gridEX[i]-pmlxStart)/(L-pmlxStart),3)
    #    eyNew[i, 1: -1] = (2.0*eps0-pmlSigmaEX*dt)/(2.0*eps0+dt*pmlSigmaEX)*eyOld[i, 1:-1]- \
    #        2.0*dt/(2.0*eps0+dt*pmlSigmaEX)/dx*(hzOld[i, 1:] - hzOld[i-1, 1:])
    #    exNew[i, 1: -1] = (2.0*eps0-dt*pmlSigmaEX)/(2.0*eps0+dt*pmlSigmaEX)*exOld[i, 1:-1]+ \
    #        2.0*dt/(2.0*eps0+dt*pmlSigmaEX)/dy * (hzOld[i, 1:] - hzOld[i, :-1])
#
#
    #for i in range(pmlIndexY,gridEY.size-1):
    #    pmlSigmaEY=pmlSigmaE0Y*pow((gridEY[i]-pmlxStart)/(L-pmlxStart),3)
    #    eyNew[1:-1, i] = (2.0*eps0-pmlSigmaEY*dt)/(2.0*eps0+dt*pmlSigmaEY)*eyOld[1:-1, i]- \
    #        2.0*dt/(2.0*eps0+dt*pmlSigmaEY)/dx*(hzOld[1:, i] - hzOld[:-1, i])
    #    exNew[1:-1, i] = (2.0*eps0-dt*pmlSigmaEY)/(2.0*eps0+dt*pmlSigmaEY)*exOld[1:-1, i]+ \
    #        2.0*dt/(2.0*eps0+dt*pmlSigmaEY)/dy * (hzOld[1:, i] - hzOld[1:, i-1])
#
    #for i in range(26,1,-1):
    #    pmlSigmaEX=pmlSigmaE0X*pow((gridEX[101-i]-pmlxStart)/(L-pmlxStart),3)
    #    eyNew[i, 1: -1] = (2.0*eps0-pmlSigmaEX*dt)/(2.0*eps0+dt*pmlSigmaEX)*eyOld[i, 1:-1]- \
    #       2.0*dt/(2.0*eps0+dt*pmlSigmaEX)/dx * (hzOld[i, 1:] - hzOld[i-1, 1:])
    #    exNew[i, 1: -1] = (2.0*eps0-dt*pmlSigmaEX)/(2.0*eps0+dt*pmlSigmaEX)*exOld[i, 1:-1]+ \
    #       2.0*dt/(2.0*eps0+dt*pmlSigmaEX)/dy * (hzOld[i, 1:] - hzOld[i, :-1])



    #for i in range(26,1,-1):
    #    pmlSigmaEY=pmlSigmaE0Y*pow((gridEY[101-i]-pmlxStart)/(L-pmlxStart),3)
    #    eyNew[1:-1, i] = (2.0*eps0-pmlSigmaEY*dt)/(2.0*eps0+dt*pmlSigmaEY)*eyOld[1:-1, i]- \
    #        2.0*dt/(2.0*eps0+dt*pmlSigmaEY)/dx*(hzOld[1:, i] - hzOld[:-1, i])
    #    exNew[1:-1, i] = (2.0*eps0-dt*pmlSigmaEY)/(2.0*eps0+dt*pmlSigmaEY)*exOld[1:-1, i]+ \
    #        2.0*dt/(2.0*eps0+dt*pmlSigmaEY)/dy * (hzOld[1:, i] - hzOld[1:, i-1])



    # PEC
    exNew[ :, 0] = 0.0;
    exNew[ :,-1] = 0.0;
    eyNew[ 0, :] = 0.0;
    eyNew[-1, :] = 0.0;  

    # --- Updates H field ---
    #for i in range(gridHX.size):
    #    for j in range(gridHX.size):
    #        hzNew[i][j] = hzOld[i][j] - cHx * (eyNew[i+1][j  ] - eyNew[i][j]) +\
    #                                    cHy * (exNew[i  ][j+1] - exNew[i][j])
    hzNew[:, :] = hzOld[:, :] - cHx * (eyNew[1:, :-1  ] - eyNew[:-1, :-1]) +\
                    cHy * (exNew[:-1  ,1:] - exNew[:-1, :-1])
      
    # ---Empieza PML ---  
    #for i in range(pmlIndexX-1,gridHX.size-1):
    #    pmlSigmaHX = pmlSigmaH0X*pow((gridHX[i]-pmlxStart)/(L-pmlxStart),3)
    #    hzNew[i, :] = (2.0*mu0-dt*pmlSigmaHX)/(2.0*mu0+dt*pmlSigmaHX)*hzOld[i, :] +\
    #       2.0*dt/(2.0*mu0+dt*pmlSigmaHX)*(exNew[i, 1:]-exNew[i, :-1])/dy -\
    #            2.0*dt/(2.0*mu0+dt*pmlSigmaHX)/dx*(eyNew[i+1  ,:-1] - eyNew[i ,:-1])
#
#
    #for i in range(pmlIndexY-1,gridHY.size-1):
    #    pmlSigmaHY = pmlSigmaH0Y*pow((gridHY[i]-pmlxStart)/(L-pmlxStart),3)
    #    hzNew[:, i] = (2.0*mu0-dt*pmlSigmaHY)/(2.0*mu0+dt*pmlSigmaHY)*hzOld[:, i] +\
    #       2.0*dt/(2.0*mu0+dt*pmlSigmaHY)*(exNew[:-1, i+1]-exNew[:-1, i])/dy -\
    #            2.0*dt/(2.0*mu0+dt*pmlSigmaHY)/dx*(eyNew[1: ,i] - eyNew[:-1 ,i])
#
#
    #for j in range(25,0-1):
    #    pmlSigmaHX = pmlSigmaH0X*pow((gridHX[100-j]-pmlxStart)/(L-pmlxStart),3)
    #    hzNew[j, :] = (2.0*mu0-dt*pmlSigmaHX)/(2.0*mu0+dt*pmlSigmaHX)*hzOld[j, :] +\
    #       2.0*dt/(2.0*mu0+dt*pmlSigmaHX)*(exNew[j, 1:]-exNew[j, :-1])/dy -\
    #            2.0*dt/(2.0*mu0+dt*pmlSigmaHX)/dx*(eyNew[j+1  ,:-1] - eyNew[j ,:-1])
#
#
    #for i in range(25,0,-1):
    #    pmlSigmaHY = pmlSigmaH0Y*pow((gridHY[100-i]-pmlxStart)/(L-pmlxStart),3)
    #    hzNew[:, i] = (2.0*mu0-dt*pmlSigmaHY)/(2.0*mu0+dt*pmlSigmaHY)*hzOld[:, i] +\
    #       2.0*dt/(2.0*mu0+dt*pmlSigmaHY)*(exNew[:-1, i+1]-exNew[:-1, i])/dy -\
    #            2.0*dt/(2.0*mu0+dt*pmlSigmaHY)/dx*(eyNew[1: ,i] - eyNew[:-1 ,i])
#

  
    # --- Updates output requests ---
    probeH[:,:,n] = hzNew[:,:]
    probeTime[n] = t
    
    # --- Updates fields and time 
    exOld[:] = exNew[:]
    eyOld[:] = eyNew[:]
    hzOld[:] = hzNew[:]
    t += dt
    print ("Time step: %d of %d" % (n, numberOfTimeSteps-1))

tictoc = time.time() - tic;
print('--- Processing finished ---')
print("CPU Time: %f [s]" % tictoc)

# ==== Post-processing ========================================================

# --- Creates animation ---
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
#ax = plt.axes(xlim=(gridE[0], gridE[-1]), ylim=(-1.1, 1.1))
ax.set_xlabel('X coordinate [m]')
ax.set_ylabel('Y coordinate [m]')
line = plt.imshow(probeH[:,:,0], animated=True, vmin=-0.5, vmax=0.5)
timeText = ax.text(0.02, 0.95, '', transform=ax.transAxes)

def init():
    line.set_array(probeH[:,:,0])
    timeText.set_text('')
    return line, timeText

def animate(i):
    line.set_array(probeH[:,:,i])
    timeText.set_text('Time = %2.1f [ns]' % (probeTime[i]*1e9))
    return line, timeText

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=nSamples, interval=50, blit=True)

plt.show()

print('=== Program finished ===')
