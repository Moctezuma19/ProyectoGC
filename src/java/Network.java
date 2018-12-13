/**
* Programa para generar Redes extendida y modificada por mí.
* @version Octubre 2016
* @author andres y Rivera Lopez Jorge Erick
*/
import java.util.Random;
import java.util.ArrayList;
import java.util.HashMap;
import java.io.IOException;
import static java.lang.Math.*;

public class Network {

    final static int S = 1;//"Susceptible";
    final static int I = 2;//"Infectado";
    final static int R = 3;//"Recuperado";

    private Node [] nodos;

    Random r;

    int maxDegree;
    int minDegree;

    public Network() {
        r = new Random();
    }
    /**
    * Metodo que recibe un numero entero mayor o igual a 2 que sera el numero de nodos de la red aleatoria
    * @param nodeNumber - número de nodos
    */
    public void createScaleFreeNetwork(int nodeNumber) {
        this.nodos = new Node [nodeNumber];
        
        Node n1 = new Node(0);
        Node n2 = new Node(1);
        
        n1.vecs.put(n2,n2);
        n2.vecs.put(n1,n1);
        nodos[0] = n1;
        nodos[1] = n2;
        
        int totalDegree = 2;

        int totalNodes = 2;
        
        for (int i = 2; i < nodeNumber; i++) {

            Node nNode = new Node(i);

            while (nNode.vecs.isEmpty()) {
                System.out.println("Conectando nodo " + i);
                for (int j = 0; j < totalNodes; j++) {
                    
                        Node cand = nodos[j];
                    
                        double pr = (cand.vecs.size() * 1.0) / totalDegree;
                        double ran = r.nextFloat();
                        if (ran < pr) {
                            nNode.vecs.put(cand,cand);
                            cand.vecs.put(nNode,nNode);
                            totalDegree+=2;
                        }
                    
                }
            }
            nodos[i] = nNode;
            totalNodes++;
        }
    }
    /**
    * Metodo para crear una red Libre de Escala, que recibe el numero de nodos y la probabilidad de conectarlos
    * @param nodeNumbre - número de nodos
    * @param connectivity - probabilidad
    */
    public void createErdosRennyNetwork(int nodeNumber, double connectivity) {

        this.nodos = new Node[nodeNumber];

        Node n = new Node(0);
        nodos[0] = n;
        int totalNodes = 1;
        
        for (int i = 1; i < nodeNumber; i++) {

            Node nNode = new Node(i);
            
            while (nNode.vecs.isEmpty()) {
                System.out.println("Conectando nodo " + i);
                for (int j = 0; j < totalNodes; j++) {

                    double ran = r.nextDouble();

                    if (ran < connectivity) {
                        
                        Node cand  = nodos[j];
                        
                        nNode.vecs.put(cand,cand);
                        cand.vecs.put(nNode,nNode);                                
                    }
                }
            }
            nodos[i] = nNode;
            totalNodes++;
        }

    }

    /**
    * Metodo que devuelve la una lista con la informacion de cada tiempo sobre el numero de S, I y R.
    * @param tiempo - tiempo en el que se esparcera la enfermedad
    * @return Una Lista con la informacion requerida
    */
    public ArrayList<int[]> esparceEnfermedad(int tiempo){
        nodos[0].estado = Network.I; // se infecta al primer nodo
        ArrayList<int[]> lista = new ArrayList<int[]>();
        for (int i = 1; i <= tiempo; i++) {
            for(Node e : nodos){
                e.actualiza();
            }
            int[] ar = obtenInformacion();
            lista.add(ar);
        }
        return lista;
    }
    /**
    * Método para obtener la información en un tiempo sobre el número de infectados, suceptibles y recuperados
    * @return un arreglo con el número de S, I y R de la red
    */
    public int[] obtenInformacion(){
        int SS = 0;
        int II = 0;
        int RR = 0;
        for(Node e : nodos){
            if(e.estado() == Network.S){
                SS++;
            } else if(e.estado() == Network.I){
                II++;
            } else {
                RR++;
            }
        }
        int[] info = {SS,II,RR};
        return info;

    }

 /**
 * Extendiendo las propiedades de cada nodo
 */   
    class Node{

        int id;
        private int estado;//estado S, I o R
        private int ultimoTI;// tiempo de infectado la ultima ocasion
        private double theta;//resistencia
        private int tempI;//tiempo de infectado
        private int tempR;//tiempo recuperado
        private double p; //probabilidad de interaccion
        public HashMap<Node,Node> vecs;//lista de vecinos
        public double rnd;
        public double rnd2;
        // Todos los nodos en un principio son susceptibles
        public Node(int id){
            this.id = id;
            this.estado = Network.S;
            this.ultimoTI=0;
            this.theta = 0.6;
            this.tempI = 0;
            this.tempR = 0;
            this.p = 0.7;
            vecs = new HashMap<Node,Node>();
        }
        
        public boolean equals(Node n){
            return this.id == n.id;
        }
        
        public boolean connectedWith(Node n){
            return vecs.containsValue(n);
        }
        //recibe el tiempo activo de la enfermedad
        public void actualiza(){
            if (this.estado == Network.I) {
                for(Node e : this.vecs.values()){
                    rnd = r.nextDouble();
                    if(rnd < p){
                        rnd2 = r.nextDouble();
                        if(rnd2 > e.theta && e.estado == Network.S){
                            e.estado = Network.I;
                            e.tempI = 0;
                        }
                    }
                }
                this.tempI++;
                if(this.tempI > 5){ 
                    this.ultimoTI = this.tempI;
                    this.tempI = 0;
                    this.estado = Network.R;
                }
            }
            // si esta recuperado, y lleva mas tiempo en este estado que el tiempo que ha estado infectado anteriormente 
            // su estado cambia a susceptible
            if(this.estado == Network.R){
                if(this.tempR > ultimoTI){
                    this.estado = Network.S;
                    this.tempR = 0;
                } else {
                    this.tempR++;
                }
            }

        }
        // Metodo que devuelve el estado de un nodo
        public int estado(){
            return this.estado;
        }

    }
}
