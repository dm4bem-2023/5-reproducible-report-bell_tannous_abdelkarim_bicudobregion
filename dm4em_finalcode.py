import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import dm4bem

w_glass=0.02
w_cocncrete_ex=0.35
w_concrete_int=0.15
w_insulation=0.15
w_oak=0.04
lambda_air=0.025
lambda_concrete=1.4
lambda_insulation=0.035
lambda_oak=0.7
lambda_glass=1.05
h_glass=15
h_air=50
h_concrete=10
h_insulation=20
h_oak=10
height=2.4
density_concrete=2300
density_glass=2500
density_insulation=20
density_oak=750
cp_glass=840
cp_concrete=880
cp_insulation=1030
cp_oak=2380
n=10 # time steps
wall_surface_sun=23.05
glass_surface_sun=3.95
α_wSW=0.2
τ_gSW=0.85
α_gSW=0.01


Va = height*10.3*9.85                   # m³, volume of air
ACH = 1                     # air changes per hour
Va_dot = ACH / 3600 * Va    # m³/s, air infiltration
# ventilation & advection
Gv = 1.204 * 1000 * Va_dot

# P-controler gain
Kp = 1e4            # almost perfect controller Kp -> ∞
#Kp = 1e-3           # no controller Kp -> 0
#Kp = 0


A = np.zeros([26, 17])       # n° of branches X n° of nodes
A[0, 0], A[0,10] = -1,1                 # branch 0: <- node 0
#A[1, 1] = -1    # branch 1: node 0 -> node 1
A[1, 1] = 1    # branch 1: node 0 -> node 1
A[2, 1], A[2, 2] = -1, 1    # branch 2: node 1 -> node 2
A[3, 2], A[3, 3] = -1, 1    # branch 3: node 2 -> node 3
A[4, 4] = 1    # branch 4: node 3 -> node 4
A[5, 4], A[5, 5] = -1, 1    # branch 5: node 4 -> node 5
#A[6, 4], A[6, 6] = -1, 1    # branch 6: node 4 -> node 6
A[6, 5], A[6, 6] = -1, 1    # branch 6: node 4 -> node 6
#A[7, 5], A[7, 6] = -1, 1    # branch 7: node 5 -> node 6
#A[7, 5], A[7, 6] = -1, -1    # branch 7: node 5 -> node 6
A[7,6]=-1 ##
A[7,7]=1
A[8, 7],A[8,3] = -1,1                 # branch 8: -> node 7
A[9, 8], A[9, 3] = -1, 1    # branch 9: node 5 -> node 7
A[10, 8] = 1                # branch 10: -> node 6
A[11, 3] = 1
A[12,3]=1# branch 11: -> node 6
A[13,3], A[13,10]=-1,1
A[14,3],A[14,9]=-1,1
A[15,9], A[15,10]=-1,1
A[16,10]=1
A[17,11]=1
A[18,10],A[18,11]=1,-1
A[19,12]=1
A[20,12]=-1
A[20,13]=1
A[21,13]=-1
A[21,14]=1
A[22,14],A[22,15]=-1,1
A[23,15],A[23,10]=-1,1
A[24,16], A[25,16],A[25,0]=1,-1,1

#Divide by 2 everywhere
#we have divided the app to two rooms: left room facing sun: upper wall-window+left wall=4.95+2.9+11=18.85
#bottom wall of the left room not facing sun: 10.47
#right room upper wall facing sun: 1.55+7.6-4.95=4.2
# right room buildings no sun: 11+5.53=16.53
G=np.zeros([26,26])
G[0,0]=2*lambda_concrete*height*16.53/w_cocncrete_ex+h_concrete*height*16.53
G[1,1]=lambda_concrete*height*10.47/w_cocncrete_ex+h_concrete*height*10.47
G[2,2]=lambda_concrete*height*10.47/w_cocncrete_ex+lambda_insulation*height*10.47/w_insulation
G[3,3]=lambda_insulation*height*10.47/w_insulation+h_insulation*height*10.47
G[4,4]=h_concrete*height*18.85
G[5,5]=lambda_concrete*height*18.85/w_cocncrete_ex
G[6,6]=lambda_concrete*height*18.85/w_cocncrete_ex+lambda_insulation*height*18.85/w_insulation
G[7,7]=lambda_concrete*height*18.85/w_cocncrete_ex
G[8,8]=h_insulation*height*18.85
G[9,9]=h_glass*height*2.5
G[10,10]=lambda_glass*height*2.5
G[11,11]=Gv
G[12,12]=Kp #perfect controller
G[13,13]=lambda_oak*height*1
G[14,14]=lambda_concrete*height*10.3/w_concrete_int+h_concrete*height*10.3
G[15,15]=lambda_concrete*height*10.3/w_concrete_int+h_concrete*height*10.3
G[16,16]=Kp
G[17,17]=lambda_glass*height*1.45
G[18,18]=h_glass*height*1.45
G[19,19]=h_concrete*height*5.65
G[20,20]=lambda_concrete*height*5.65/w_cocncrete_ex
G[21,21]=lambda_concrete*height*5.65/w_cocncrete_ex+lambda_insulation*height*5.65/w_insulation
G[22,22]=lambda_concrete*height*5.65/w_cocncrete_ex
G[23,23]=h_insulation*height*5.65
G[24,24]=lambda_concrete*height*16.53/w_cocncrete_ex+h_concrete*height*16.53
G[25,25]=lambda_concrete*height*16.53/w_cocncrete_ex+lambda_insulation*height*16.53/w_insulation

#Capacties
#Matrix C
C=np.zeros([17,17])
C[1,1]=cp_concrete*density_concrete*height*w_cocncrete_ex*10.47
C[2,2]=cp_insulation*density_insulation*height*w_insulation*10.47
C[5,5]=cp_concrete*density_concrete*height*w_cocncrete_ex*18.85
C[6,6]=cp_insulation*density_insulation*height*w_insulation*18.85
C[8,8]=cp_glass*density_glass*height*w_glass*2.5
C[9,9]=cp_oak*density_oak*height*w_oak
C[8,8]=cp_glass*density_glass*height*w_glass*1.45
C[13,13]=cp_concrete*density_concrete*height*w_cocncrete_ex*5.65
C[14,14]=cp_insulation*density_insulation*height*w_insulation*5.65
C[16,16]=cp_concrete*density_concrete*height*w_cocncrete_ex*16.53
C[0,0]=cp_insulation*density_insulation*height*w_insulation*16.53

b=np.zeros([26])
b[[1,4,10,11,12,16,17,19,24]]=1

f=np.zeros([17])
f[[4,7,8,11,12,15]]=1

y=np.zeros([17])
y[[3,10]]=1

#################################################################
#Steady State

[As, Bs, Cs, Ds] = dm4bem.tc2ss(A, G, b, C, f, y)
print('As = \n', As, '\n')
print('Bs = \n', Bs, '\n')
print('Cs = \n', Cs, '\n')
print('Ds = \n', Ds, '\n')

b=np.zeros([26]) # temperature sources
b[[1,4,10,17,19,24]]=10 # outdoor temperature
b[[11,12,16]] = 20            # indoor set-point temperature

f = np.zeros(17)         # flow-rate sources

θ = np.linalg.inv(A.T @ G @ A) @ (A.T @ G @ b + f)
print(f'θ = {θ} °C')
q = G @ (-A @ θ + b)
print(f'q = {q} W')

#####################################
bT=np.array([10,10,10,20,20,20,10,10,10])
fQ=np.array([0,0,0,0,0,0])
u = np.hstack([bT, fQ])
print(f'u = {u}')

yss = (-Cs @ np.linalg.inv(As) @ Bs + Ds) @ u
print(f'yss = {yss} °C')


print(f'Max error between DAE and state-space, room 1: \
{max(abs(θ[3] - yss)):.2e} °C')
print(f'Max error between DAE and state-space, room 2: \
{max(abs(θ[10] - yss)):.2e} °C')


#################################################################################################
#Dynamic simulation
λ = np.linalg.eig(As)[0]    # eigenvalues of matrix As
print('Time constants: \n', -1 / λ, 's \n')
print('2 x Time constants: \n', -2 / λ, 's \n')
dtmax = 2 * min(-1. / λ)
print(f'Maximum time step: {dtmax:.2f} s = {dtmax / 60:.2f} min')

# time step
dt = np.floor(dtmax / 60) * 20   # s
print(f'dt = {dt} s = {dt / 60:.0f} min')

# settling time
time_const = np.array([int(x) for x in sorted(-1 / λ)])
print('4 * Time constants: \n', 4 * time_const, 's \n')

t_settle = 4 * max(-1 / λ)
print(f'Settling time: \
{t_settle:.0f} s = \
{t_settle / 60:.1f} min = \
{t_settle / (3600):.2f} h = \
{t_settle / (3600 * 24):.2f} days')

# Step response
# -------------
# Find the next multiple of 3600 s that is larger than t_settle
duration = np.ceil(t_settle / 3600) * 3600
n = int(np.floor(duration / dt))    # number of time steps
t = np.arange(0, n * dt, dt)        # time vector for n time steps

print(f'Duration = {duration} s')
print(f'Number of time steps = {n}')
# pd.DataFrame(t, columns=['time'])


#Input vector

u = np.zeros([15, n])                # u = [To To To Tisp Φo Φi Qa Φa]
u[0:3,:] = 10 * np.ones([3, n])    # To = 10 for n time steps
u[3:6,:] = 20 * np.ones([3, n])      # Tisp = 20 for n time steps
u[6:9, :] = 10 * np.ones([3, n])    # To = 10 for n time steps
print(u)

# pd.DataFrame(u)

#####################
#Time integration
n_s = As.shape[0]                      # number of state variables
θ_exp = np.zeros([n_s, t.shape[0]])    # explicit Euler in time t
θ_imp = np.zeros([n_s, t.shape[0]])    # implicit Euler in time t

I = np.eye(n_s)                        # identity matrix

for k in range(n - 1):
    θ_exp[:, k + 1] = (I + dt * As) @\
        θ_exp[:, k] + dt * Bs @ u[:, k]
    θ_imp[:, k + 1] = np.linalg.inv(I - dt * As) @\
        (θ_imp[:, k] + dt * Bs @ u[:, k])

#explicit and  implicit Euler methods
y_exp = Cs @ θ_exp + Ds @  u
y_imp = Cs @ θ_imp + Ds @  u
#Results and Plot
fig, ax = plt.subplots()
ax.plot(t / 3600, y_exp.T, t / 3600, y_imp.T)
ax.set(xlabel='Time, $t$ / h',
       ylabel='Temperatue, $θ_i$ / °C',
       title='Step input: outdoor temperature $T_o$')
ax.legend(['Explicit', 'Implicit'])
ax.grid()
plt.show()

#The value the indoor temperature obtained after the settling time is almost equal to the value obtained in steady-state.
print('Steady-state indoor temperature obtained with:')
print(f'- DAE model: Room 1: {float(θ[3]):.4f}, Room 2: {float(θ[10]):.4f} °C')
print(f'- state-space model: Room 1: {float(yss[0]):.4f}, Room 2: {float(yss[1]):.4f} °C')
print(f'- steady-state response to step input: Room 1: {float(y_exp[0, -2]):.4f}, Room 2: {float(y_exp[1, -2]):.4f} °C')



#4. Simulate response to weather
#Define start and end time.
start_date = '01-03 12:00:00'
end_date = '02-05 18:00:00'
start_date = '2000-' + start_date
end_date = '2000-' + end_date
print(f'{start_date} \tstart date')
print(f'{end_date} \tend date')

#Inputs
#Read weather data
filename = 'FRA_Lyon.074810_IWEC.epw'
[data, meta] = dm4bem.read_epw(filename, coerce_year=None)
weather = data[["temp_air", "dir_n_rad", "dif_h_rad"]]
del data

#select from the start to the end date:
weather.index = weather.index.map(lambda t: t.replace(year=2000))
weather = weather.loc[start_date:end_date]

#Solar irradiance on the walls
#calculate the: direct irradiance, W/m²,diffuse irradiance, W/m²,reflected irradiance, W/m²
surface_orientation = {'slope': 90,
                       'azimuth': 0,
                       'latitude': 45}
albedo = 0.2
rad_surf = dm4bem.sol_rad_tilt_surf(
    weather, surface_orientation, albedo)
# pd.DataFrame(rad_surf)

#total solar irradiance Etot
rad_surf['Φtot'] = rad_surf.sum(axis=1)

#Resample the weather data
data = pd.concat([weather['temp_air'], rad_surf['Φtot']], axis=1)
data = data.resample(str(dt) + 'S').interpolate(method='linear')
data = data.rename(columns={'temp_air': 'To'})
# pd.DataFrame(data)

#Other inputs
#consider indoor temperature Ti,sp=20∘C and the auxiliary heat flow Q˙a=0W
data['Ti'] = 20 * np.ones(data.shape[0])
data['Qa'] = 0 * np.ones(data.shape[0])
# pd.DataFrame(data)

#Input vector in time
# input vector
To = data['To']
Ti = data['Ti']
Φo1 = α_wSW * 18.85 *height * data['Φtot']
Φo2= α_wSW * 4.2 * height * data['Φtot']
Φa1=α_gSW * 2.5 * height * data['Φtot']
Φa2=α_gSW * 1.45 *height * data['Φtot']
Φi1 = τ_gSW * α_wSW *height *2.5 * data['Φtot']
Φi2 = τ_gSW * α_wSW * height* 1.45 * data['Φtot']
Qa = data['Qa']
#Φa = α_gSW * glass_surface_sun * data['Φtot']

#Heat sources 4,7,8,11,12,15
u = pd.concat([To, To, To, Ti,Ti,Ti,To,To,To, Φo1, Φo2,Φi1,Φi2, Φa1, Φa2 ], axis=1)
u.columns.values[[9, 10, 11, 12, 13, 14]] = ['Φo1', 'Φo2','Φi1','Φi2', 'Φa1', 'Φa2']
# pd.DataFrame(u)

#Initial conditions
θ_exp = 20 * np.ones([As.shape[0], u.shape[0]])

#Time integration
for k in range(u.shape[0] - 1):
    θ_exp[:, k + 1] = (I + dt * As) @ θ_exp[:, k]\
        + dt * Bs @ u.iloc[k, :]

y_exp = Cs @ θ_exp + Ds @ u.to_numpy().T
q_HVAC = Kp * (data['Ti'] - y_exp[0, :])
data['θi_exp1'] = y_exp.T[:,0]
data['θi_exp2']=y_exp.T[:,1]
data['q_HVAC'] = q_HVAC.T

fig, axs = plt.subplots(2, 1)

data[['To', 'θi_exp1']].plot(ax=axs[0],
                            xticks=[],
                            ylabel='Temperature, $θ$ / °C')
axs[0].legend(['$θ_{outdoor}$', '$θ_{indoor}$'],
              loc='upper right')

data[['Φtot', 'q_HVAC']].plot(ax=axs[1],
                              ylabel='Heat rate, $q$ / W')
axs[1].set(xlabel='Time')
axs[1].legend(['$Φ_{total}$', '$q_{HVAC}$'],
             loc='upper right')
plt.show()


#imulation in free-running with weather data using Euler explicit method of integration.
# a) Indoor and outdoor temperatures. b) Solar and HVAC heat flow rates.
t = dt * np.arange(data.shape[0])   # time vector

fig, axs = plt.subplots(2, 1)
# plot outdoor and indoor temperature
axs[0].plot(t / 3600 / 24, data['To'], label='$θ_{outdoor}$')
axs[0].plot(t / 3600 / 24, y_exp[0, :], label='$θ_{indoor}$')
axs[0].set(ylabel='Temperatures, $θ$ / °C',
           title='Simulation for weather')
axs[0].legend(loc='upper right')

# plot total solar radiation and HVAC heat flow
axs[1].plot(t / 3600 / 24, data['Φtot'], label='$Φ_{total}$')
axs[1].plot(t / 3600 / 24, q_HVAC, label='$q_{HVAC}$')
axs[1].set(xlabel='Time, $t$ / day',
           ylabel='Heat flows, $q$ / W')
axs[1].legend(loc='upper right')

fig.tight_layout()