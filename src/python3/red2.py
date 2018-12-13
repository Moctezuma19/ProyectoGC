from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
import Bio.pairwise2
import matplotlib.pyplot as plt
import numpy as np
import random as r
import somoclu
'''
Genomica Computacional 2019-1
Proyecto
Rivera Lopez Jorge Erick
312121645
version 2.5
'''
diccio = {0:'A',1:'T',2:'C',3:'G'}
'''
Metodo que genera una secuencia de una longitud
especifica.
'''
def generaSecuencia(longitud) :
    sec = ''
    for i in range(longitud) :
        sec += diccio[int(np.random.choice(4,1))]
    return sec
'''
    Funcion que realiza una mutacion de sustitucion
    en una secuencia de nucleotidos.
    secuencia - secuencia  
    proba1 - real variable
    proba2 - probabilidad fija
    no_mut - numero de nucleotidos a mutar
'''
def mutaS(secuencia, no_mut) :
    aux = ''
    elem = np.random.choice(len(secuencia),no_mut)
    for i in range(len(secuencia)) :
        if i in elem :
            aux += diccio[int(np.random.choice(4,1))]
            continue
        else :
            aux += secuencia[i]
    return aux
'''
    Funcion que realiza una mutacion de delecion
    en una secuencia de nucleotidos.
    secuencia - secuencia
    proba1 - real variable
    proba2 - probabilidad fija
    no_mut - numero de nucleotidos a mutar
'''
def mutaD(secuencia, no_mut) :
    aux = ''
    elem = np.random.choice(len(secuencia),no_mut)
    i = 0
    while i < len(secuencia) :
        if not i in elem :
            aux += secuencia[i]
        i += 1
    return aux

def abre() :
    with open('sequence.fasta','rU') as entrada:
        dna = list(SeqIO.parse(entrada,'fasta'))
        return dna

poblacion = abre()
S = []
I = []
R = []
MUT = []
MUTAA = []
GC = []
GC2 = []
'''
Funcion que grafica la distribucion SIR en un lapso de Tiempo
'''
def graficaDistribucionSIRTiempo() :
    cgc = range(len(S))
    plt.plot(cgc,S,'r',cgc,I,'g',cgc,R,'b')
    ##plt.title("Distribucion de SIR en el tiempo")
    plt.xlabel('Tiempo')
    plt.ylabel('Cantidad de organismos')
    plt.legend(['Susceptible', 'Infectado', 'Recuperado'])
    plt.show()
'''
Funcion que grafica la distribucion de cambios de nucleotidos
'''
def graficaMutacionNT() :
    cgc = range(len(MUT))
    plt.bar(cgc,MUT)
    #plt.title("Cantidad de cambios en NT en cada organismo")
    plt.xlabel('Organismo')
    plt.ylabel('Cantidad de mutaciones(S)')
    plt.show()
'''
Funcion que grafica la distribucion de cambios de Aminoacidos
finales en un organismo
'''
def graficaMutacionAAFinal() :
    cgc = range(len(MUTAA))
    plt.bar(cgc,MUTAA)
    #plt.title("Distribucion de mutacion AA Final")
    plt.xlabel('Organismo')
    plt.ylabel('Cantidad de cambios en AA final')
    plt.show()
'''
Funcion para graficar el contenido GC inicial
'''
def graficaContenidoGCI() :
    cgc = range(len(GC))
    plt.bar(cgc,GC)
    #plt.title("Contenido GC en los organismos(TI)")
    plt.xlabel('Organismo')
    plt.ylabel('Contenido GC')
    plt.show()
'''
Funcion para graficar el contenido GC final
'''
def graficaContenidoGCF() :
    cgc = range(len(GC2))
    plt.bar(cgc,GC2)
    #plt.title("Contenido GC en los organismos(TF)")
    plt.xlabel('Organismo')
    plt.ylabel('Contenido GC')
    plt.show()

def graficaDistribucionGC() :
    GC.sort()
    d = {}
    for x in GC :
        try :
            d[x] += 1
        except :
            d.setdefault(x,1)
    cgc = range(len(d.keys()))
    plt.bar(cgc,list(d.values()))
    #plt.title("Distribucion de Contenido GC (TI)")
    plt.xlabel('Contenido GC')
    plt.ylabel('Cantidad de organismos')
    plt.xticks(cgc,list(d.keys()), rotation = 20)
    plt.show()

def graficaDistribucionGCF() :
    GC2.sort()
    d = {}
    for x in GC2 :
        try :
            d[x] += 1
        except :
            d.setdefault(x,1)
    cgc = range(len(d.keys()))
    plt.bar(cgc,list(d.values()))
    #plt.title("Distribucion de Contenido GC (TF)")
    plt.xlabel('Contenido GC')
    plt.ylabel('Cantidad de organismos')
    k =list(d.keys())
    plt.xticks([0,int(len(k)/3),int(len(k)/3)*2,len(k)-1],[k[0],k[int(len(k)/4)],k[int(len(k)/4)*2],k[-1]])
    plt.show()

def graficaDistribucionMNT() :
    MUT.sort()
    d = {}
    for x in MUT :
        try :
            d[x] += 1
        except :
            d.setdefault(x,1)
    cgc = range(len(d.keys()))
    plt.bar(cgc,list(d.values()))
    #plt.title("Distribucion de numero de mutaciones total")
    plt.xlabel('Cantidad de cambios en NT')
    plt.ylabel('Cantidad organismos')
    k=list(d.keys())
    plt.xticks([0,int(len(k)/3),int(len(k)/3)*2,len(k)-1],[k[0],k[int(len(k)/4)],k[int(len(k)/4)*2],k[-1]])
    plt.show()

def graficaDistribucionMAA() :
    MUTAA.sort()
    d = {}
    for x in MUTAA :
        try :
            d[x] += 1
        except :
            d.setdefault(x,1)
    cgc = range(len(d.keys()))
    plt.bar(cgc,list(d.values()))
    #plt.title("Distribucion de numero de cambios en AA (TF)")
    plt.xlabel('Cantidad de cambios en AA')
    plt.ylabel('Cantidad de organismos')
    k =list(d.keys())
    plt.xticks([0,int(len(k)/3),int(len(k)/3)*2,len(k)-1],[k[0],k[int(len(k)/4)],k[int(len(k)/4)*2],k[-1]])
    plt.show()   

'''
Funcion para normalizar los resultados
'''
def normaliza() :
    prom = 0.0

    for x in GC2 :
        prom += x

    prom /= len(GC2)

    for x in GC2 :
        x /= prom

    prom = 0.0

    for x in MUT :
        prom += x

    prom /= len(MUT)

    for x in MUT:
        x /= prom

    prom = 0.0

    for x in MUTAA :
        prom += x

        prom /= len(MUTAA)
    for x in MUTAA :
        x /= prom
'''
Funcion para una grafica SOM
'''
def grafica_som():
    normaliza()
    lista = []
    for i in range(len(GC2)) :
        lista.append([MUT[i],MUTAA[i],GC2[i]])

    data = np.float32(np.array(lista))
    som = somoclu.Somoclu(100, 100, data=data)
    som.train()
    som.view_component_planes()

'''
Funcion cuenta el contenido GC de una secuencia
'''
def factorGC(secuencia):
    gc = 0
    for x in secuencia :
        if x == 'G' or x == 'C' :
            gc += 1

    return gc/len(secuencia)
'''
Clase nodo para una red
'''
class Nodo :
    def __init__(self,id,secuencia):
        self.id = id
        self.estado = 1 #se marca susceptible
        self.theta = 0.7 # resistencia a la enfermedad
        self.theta2 = 0.18 #
        self.cambios_nt = 0 #cuenta cambios en los nucleotidos
        self.p = 0.7 # probabilidad de contacto
        self.vecs = {} # vecinos
        self.secuencia = secuencia
        self.tempR = 0
        self.tempI = 0
        self.ultimoTI = 0

    def __hash__(self) :
        return self.id
    def __eq__(self,nodo) :
        self.id = nodo.id

    def connectedWith(self,nodo) :
        return nodo in self.vecs.values()
    '''
    Funcion que actualiza el estado de un nodo
    '''
    def actualiza(self) :
        if self.estado == 2 :
            for e in self.vecs.values() :
                rnd = r.random()
                if rnd < self.p :
                    rnd2 = r.random()
                    if rnd2 > e.theta and e.estado == 1 :
                        e.estado = 2
                        #----->e.estado = 3 # adquiere una inmunidad temporal
            if r.random() < self.theta:
                self.secuencia = mutaS(self.secuencia,1)
                self.cambios_nt += 1
                
            self.tempI += 1
            if self.tempI > 5 :
                if r.random() > self.theta2 :
                    return #self.tempI = 0
                else :
                    self.ultimoTI = self.tempI
                    self.tempI = 0
                    self.estado = 3
        
        if self.estado == 3 :
            if self.tempR > self.ultimoTI :
                self.estado = 1
                self.tempR = 0
            else :
                self.tempR +=1


'''
Clase para crear una red
'''
class RLE :

    def __init__(self) :
        self.nodos = []
        self.s = 1 # susceptible
        self.i = 2 # infectado
        self.r = 3 # recuperado
    '''
    Funcion que crea una red libre de escala con un orden
    especificado
    '''
    def createScaleFreeNetwork(self,nodeNumber) :
        for x in range(nodeNumber) :
            self.nodos.append(None)
        n1 = Nodo(0,str(poblacion[0].seq))
        n2 = Nodo(1,str(poblacion[1].seq))
        n1.vecs.setdefault(n2,n2)
        n2.vecs.setdefault(n1,n1)
        self.nodos[0] = n1
        self.nodos[1] = n2
        totalDegree = 2
        totalNodes = 2
        for i in range(2,nodeNumber) :
            nNode = Nodo(i,str(poblacion[i].seq))
            while len(nNode.vecs) == 0 :
                print("Conectando nodo " + str(i))
                for j in range(totalNodes) :
                    cand = self.nodos[j]
                    pr = (len(cand.vecs.values()) * 1.0) / totalDegree
                    ran = r.random()
                    if ran < pr and not nNode.connectedWith(cand) :
                        nNode.vecs.setdefault(cand,cand)
                        cand.vecs.setdefault(nNode,nNode)
                        totalDegree += 2
            self.nodos[i] = nNode
            totalNodes += 1
    '''
    Funcion que crea una red aleatoriocon un orden
    especificado y una probabilidad de conexion dada
    '''            
    def createRandomNetwork(self, nodeNumber,connectivity) :
        for x in range(nodeNumber) :
            self.nodos.append(None)
        n = Nodo(0,str(poblacion[0].seq))
        self.nodos[0] = n;
        totalNodes = 1;
        for i in range(1,nodeNumber) :
            nNode = Nodo(i,str(poblacion[i].seq))
            while len(nNode.vecs) == 0 :
                print("Conectando nodo " + str(i))
                for j in range(0,totalNodes) :
                    ran = r.random();
                    cand  = self.nodos[j]
                    if ran < connectivity and not nNode.connectedWith(cand):                      
                        nNode.vecs.setdefault(cand,cand)
                        cand.vecs.setdefault(nNode,nNode)
            self.nodos[i] = nNode;
            totalNodes += 1;
    
    '''
    Funcion que crea una red aleatoria con nodos de una red
    previemente creada y una probabilidad de conexion dada
    '''            
    def convierteRandomNetwork(self,connectivity) :
        for x in self.nodos :
            x.vecs.clear()

        totalNodes = 1;
        for i in range(1,len(self.nodos)) :
            nNode = self.nodos[i]
            while len(nNode.vecs) == 0 :
                print("Reconectando nodo " + str(i))
                for j in range(0,totalNodes) :
                    ran = r.random();
                    cand  = self.nodos[j]
                    if ran < connectivity and not nNode.connectedWith(cand):                      
                        nNode.vecs.setdefault(cand,cand)
                        cand.vecs.setdefault(nNode,nNode)
            totalNodes += 1;
        '''
    Funcion que crea una red libre de escala con nodos de una red
    previemente creada y una probabilidad de conexion dada
    '''    
    def convierteScaleFreeNetwork(self) :
        for x in self.nodos :
            x.vecs.clear()
        n1 = self.nodos[0]
        n2 = self.nodos[1]
        n1.vecs.setdefault(n2,n2)
        n2.vecs.setdefault(n1,n1)
        totalDegree = 2
        totalNodes = 2
        for i in range(2,len(self.nodos)) :
            nNode = self.nodos[i]
            while len(nNode.vecs) == 0 :
                print("Reconectando nodo " + str(i))
                for j in range(totalNodes) :
                    cand = self.nodos[j]
                    pr = (len(cand.vecs.values()) * 1.0) / totalDegree
                    ran = r.random()
                    if ran < pr and not nNode.connectedWith(cand) :
                        nNode.vecs.setdefault(cand,cand)
                        cand.vecs.setdefault(nNode,nNode)
                        totalDegree += 2
            totalNodes += 1      
    '''
    Funcion para contar cambios finales en aminoacidos
    seq1 - secuencia original
    seq2 - secuencia mutada
    '''
    def cuentaCambiosFinalAA(self,seq1,seq2) :
        mayor = ''
        menor = ''
        if len(str(seq1)) < len(str(seq2)) :
            mayor = str(seq2)
            menor = str(seq1)
        elif len(str(seq1)) > len(str(seq2)) :
            menor = str(seq2)
            mayor = str(seq1)
        else :
            menor = str(seq2)
            mayor = str(seq1)

        cuenta = 0
        for i in range(len(menor)):
            if str(menor)[i] != str(mayor)[i] :
                cuenta += 1
        return cuenta + len(mayor) - len(menor)

    '''
    Metodo para contar cambios en aminoacidos
    '''
    def alinea(self, lista1) :
        for i in range(len(lista1)) :
            #x = Bio.pairwise2.align.localxx(lista1[i], self.nodos[i].secuencia)
            a = Seq(lista1[i],generic_dna).translate()
            b = Seq(self.nodos[i].secuencia,generic_dna).translate()
            MUTAA.append(self.cuentaCambiosFinalAA(a,b))
    '''
    Funcion para contar el factor de GC de cada organismo
    '''
    def cuentaFactorGCF(self) :
        for x in self.nodos :
            GC2.append(factorGC(x.secuencia))

    '''
    Funcion para contar el factor GC de cada organismo
    '''
    def cuentaFactorGCI(self,copia) :
        for x in copia :
            GC.append(factorGC(x))
    '''
    Funcion que acutaliza la lista de cuenta de cambios
    de nucleotidos
    '''
    def cuentaMutaciones(self) :
        for x in self.nodos :
            MUT.append(x.cambios_nt)
    '''
    Funcion que actualiza la lista de numero
    de organismos susceptibles, infectados y recuperados
    '''
    def obtenInformacion(self) :
        SS = 0
        II = 0
        RR = 0
        for e in self.nodos :
            if e.estado == 1 :
                SS += 1
            elif e.estado == 2 :
                II += 1
            else :
                RR += 1
        S.append(SS)
        I.append(II)
        R.append(RR)


    '''
    Funcion que esparce una enfermedad en una red, 
    en un tiempo especificado.
    '''
    def esparceEnfermedad(self, tiempo) :
        self.nodos[0].estado = 2 # se infecta a un nodo
        lista = []
        respaldo = []
        for x in self.nodos :
            respaldo.append(x.secuencia)
        for i in range(1,tiempo + 1) :
            
            if int(tiempo/3) == i :
                self.convierteScaleFreeNetwork()
            if int(tiempo/3) * 2 == i :
                self.convierteRandomNetwork(0.1)
            
            for n in self.nodos :
                n.actualiza()
            self.obtenInformacion()
        self.cuentaMutaciones()
        self.alinea(respaldo)
        self.cuentaFactorGCI(respaldo)
        self.cuentaFactorGCF()
        return lista

    def createRegularNetwork(self, nodeNumber) :
        print(nodeNumber)
        print("Creando red regular . . .")
        for x in range(nodeNumber) :
            self.nodos.append(Nodo(x,str(poblacion[x].seq)))
        #creando anillo
        for x in range(nodeNumber - 1) :
            self.nodos[x].vecs.setdefault(self.nodos[x + 1],self.nodos[x + 1])
            self.nodos[x + 1].vecs.setdefault(self.nodos[x],self.nodos[x])
        self.nodos[nodeNumber - 1].vecs.setdefault(self.nodos[0],self.nodos[0])
        self.nodos[0].vecs.setdefault(self.nodos[nodeNumber -1],self.nodos[nodeNumber-1])
        #creando grafica regular
        cont = 0
        print(" . . . ")
        aux = nodeNumber
        '''
        if nodeNumber % 2 != 0 :
            aux -= 1
            '''
        while cont < nodeNumber-4 :
            print(cont)
            self.nodos[cont].vecs.setdefault(self.nodos[cont+2],self.nodos[cont+2])
            self.nodos[cont + 1].vecs.setdefault(self.nodos[cont+3],self.nodos[cont+3])
            cont += 4
        print("Red Regular creada")




red = RLE()
red.createRegularNetwork(len(poblacion))
#red.createScaleFreeNetwork(len(poblacion))
#red.createRandomNetwork(len(poblacion),0.1)
red.esparceEnfermedad(100)
graficaDistribucionSIRTiempo()
graficaMutacionNT()
graficaMutacionAAFinal()
graficaContenidoGCI()
graficaContenidoGCF()
graficaDistribucionGC()
graficaDistribucionGCF()
graficaDistribucionMNT()
graficaDistribucionMAA()
grafica_som()