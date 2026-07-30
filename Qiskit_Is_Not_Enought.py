import numpy as np

i=1j

ket_0 = np.array([[1], [0]]);
ket_1 = np.array([[0], [1]]);

bra_0 = ket_0.conj().T;
bra_1 = ket_1.conj().T;

ket_00 = np.kron(ket_0, ket_0);
ket_01 = np.kron(ket_0, ket_1);
ket_10 = np.kron(ket_1, ket_0);
ket_11 = np.kron(ket_1, ket_1);

bra_00 = ket_00.conj().T;
bra_01 = ket_01.conj().T;
bra_10 = ket_10.conj().T;
bra_11 = ket_11.conj().T;

#====== Compuertas Unitarias

Id2= np.array([[1, 0], [0, 1]]); 

Pauli_X= np.array([[0, 1], [1, 0]]);
Pauli_Y= np.array([[0, -i], [i, 0]]);
Pauli_Z= np.array([[1, 0], [0, -1]]);

Ph_S = np.array([[1, 0], [0, i]]);
Ph_T = np.array([[1, 0], [0, np.exp(i*np.pi/4)]])

def PhaseR_X(angX):
    theta=angX;
    R_Xp = np.array([[np.cos(theta/2), -i*np.sin(theta/2)], [-i*np.sin(theta/2), np.cos(theta/2)]]);
    return R_Xp
def PhaseR_Y(angY):
    betta=angY;
    R_Yp = np.array([[np.cos(betta/2), -1*np.sin(betta/2)], [np.sin(betta/2), np.cos(betta/2)]]);
    return R_Yp
def PhaseR_Z(angZ):
    gamma=angZ;
    R_Zp = np.array([[np.exp(-i*gamma/2), 0], [0, np.exp(i*gamma/2)]])
    return R_Zp

Hadamard= 1/np.sqrt(2)*(np.array([[1, 1], [1, -1]]));

#========== Compuertas de sistemas Bipartitos

cx = np.array([[1, 0, 0, 0], 
               [0, 1, 0, 0],
               [0, 0, 0, 1], 
               [0, 0, 1, 0]]);

cz = np.array([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 1, 0],
               [0, 0, 0, -1]]);

cy = np.array([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 0, -i],
               [0, 0, i, 0]]);

ch = np.array([[1, 0, 0, 0],
               [0, 1, 0, 0],
               [0, 0, 1/np.sqrt(2), 1/np.sqrt(2)],
               [0, 0, 1/np.sqrt(2), -1/np.sqrt(2)]]);

swap = np.array([[1, 0, 0, 0],
                 [0, 0, 1, 0],
                 [0, 1, 0, 0],
                 [0, 0, 0, 1]]);

#========== Funciones Adicionales

def rho(vec):
    return np.kron(vec, vec.conj().T);

#=========== Post-Procesamiento

def statevector(stvec):
    flat = stvec.flatten()
    n_qubits = int(np.log2(len(flat)))
    terms = []
    for idx, amplitude in enumerate(flat):
        if not np.isclose(amplitude, 0):
            binary = format(idx, f'0{n_qubits}b')
            terms.append(f"({amplitude:.4f})|{binary}〉")
    print(" + ".join(terms))

#=============== Interfaz

def NewCircuit(noq, iqest, nocb):
    #globals().nocb
    if (iqest == 0):
        q = [ket_0.copy() for _ in range(noq)]
    if (iqest == 1):
        q = [ket_1.copy() for _ in range(noq)]
    q[0] = Id2 @ q[0]
    return q

#=================== Test Area

print("Script succesfully started")
noq = int(input("Number of qubits of your circuit: "))
iqest = int(input("Initial qubits state (0 or 1): "))
nocb = int(input("Number of classical bits of your circuit: "))
nc=NewCircuit(noq, iqest, nocb)
#res = Hadamard@nc[2]
statevector(nc[1])