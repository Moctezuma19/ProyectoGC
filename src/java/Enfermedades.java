
/**
* Programa que ejecuta el esparcimiento del Virus de la Gripe en una Red Libre de Escala y Aleatoria
* Computación Distribuida
* Práctica 2 - Modelo SIR
* @author Rivera Lopez Jorge Erick
* @version Octubre 2016
*/
import java.util.ArrayList;
import java.io.IOException;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.File;

public class Enfermedades{

	public static void main(String[] args) {

		Network f = new Network();
		Network e = new Network();

		try {

			f.createScaleFreeNetwork(100000);
			ArrayList<int[]> g1 = f.esparceEnfermedad(30);

			BufferedWriter bw = new BufferedWriter(new FileWriter(new File("./EnfermedadInfo.txt")));
			bw.write("VIRUS DE LA GRIPE");
			
			bw.newLine();
			
			bw.write("Red Libre de escala");
			bw.newLine();
			for(int[] i : g1){
				bw.write(i[0]+"      "+i[1]+" "+i[2]);
				bw.newLine();
			}
			
			e.createErdosRennyNetwork(100000,0.01);
			ArrayList<int[]> g2 = e.esparceEnfermedad(30);

			bw.newLine();
			
			bw.write("Red Aleatoria");
			bw.newLine();
			for(int[] i : g2){
				bw.write(i[0]+"      "+i[1]+" "+i[2]);
				bw.newLine();
			}
			
			bw.close();

		} catch(IOException k){
			System.out.println("Error  de I|O");
		}
	}
}