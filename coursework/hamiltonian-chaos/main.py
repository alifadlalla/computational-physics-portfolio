from functions import*

K_values = [0.1, 0.25, 0.5, 0.9, 2, 4]  # Different values of K
M = 20  # Size of initial grid
N = 2000  # Number of iterations
L = 1000  # Size of array for phase space
extension = 'svg' # for saving the plots onle pdf or svg are allowed as input


# phase space for k=0.9
iterate(M,N,L,K=0.9,ext=extension,lyapov=False)

# lopp over different K
#for K in K_values:
#    iterate(M,N,L,K,ext=extension,lyapov=False)

# phase space k=0.9 with lyapunov
#iterate(M,N,L,K=0.9,ext=extension,lyapov=True)
