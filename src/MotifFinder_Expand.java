import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

public class MotifFinder_Expand {
	
	public static final double windowSize = 20;
	
	public static double aminoAcidRatio;
	public static void main(String[] args){
		aminoAcidRatio = Double.parseDouble(args[1]); //parse amino acid ratio from command line argument
		String inputfileName = args[0]; //parse fasta file from command line arguments
		Map<String,String> ProteinNameToSequences = readProteinSequences(inputfileName); //read sequences and protein names into Map
		Map<String,Set<Integer>> ProteinsToMotifs = FindMotifs(ProteinNameToSequences); //calculate amino acid ratio in protein set
		
		PrintSeqsWithMotifs(ProteinsToMotifs, ProteinNameToSequences, inputfileName); //write proteins that pass aa ratio threshold
		
	}
	
	/*
	 * function printing the sequences with the motifs identified
	 */
	
	private static void PrintSeqsWithMotifs(Map<String,Set<Integer>> ProteinsToMotifs,Map<String,String> ProteinNameToSequences, String inputFileName) {
		String outputFileName = inputFileName.substring(0,inputFileName.length()-6) + "_Result_Ratio_"+aminoAcidRatio+".txt";
		int counterProteins = 0;
		
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFileName)));
			for(Entry<String, Set<Integer>> e: ProteinsToMotifs.entrySet()){
				if(e.getValue().size()>0) {
					counterProteins++;
					out.write(e.getKey()+"\n");
					out.flush();
					int counter = 1;
					for(int i = 0; i < ProteinNameToSequences.get(e.getKey()).length(); i++) {
						if (counter == 50) {
							out.write("\t");
							for(int j = i-49; j <= i; j++) {
								if(e.getValue().contains(j)){
									out.write("*");
								}
								else {
									out.write(" ");
								}
							}
							out.write("\n");
							out.flush();
							out.write((i-49)+"\t");
							for(int j = i-49; j <= i; j++) {
								out.write(ProteinNameToSequences.get(e.getKey()).charAt(j)+"");
							}
							out.write("\n");
							out.flush();
							counter = 0;
						}
						counter++;
					}
					int startIndexReminder =  (int)Math.floor((double)(ProteinNameToSequences.get(e.getKey()).length()/50))*50;
					out.write("\t");
					for(int i = startIndexReminder; i < ProteinNameToSequences.get(e.getKey()).length(); i++) {
						
						if(e.getValue().contains(i)){
							out.write("*");
						}
						else {
							out.write(" ");
						}
					}
					out.write("\n");
					out.flush();
					out.write(startIndexReminder+"\t");
					for(int i = startIndexReminder; i < ProteinNameToSequences.get(e.getKey()).length(); i++) {
						out.write(ProteinNameToSequences.get(e.getKey()).charAt(i)+"");
					}
					out.write("\n\n");
					out.flush();
				}
				
			}
			out.close();
		}catch(Exception e) {
			e.printStackTrace();
		}
		System.out.println(counterProteins);
	}
	
	/*
	 * Function that finds the motifs
	 */
	private static Map<String,Set<Integer>> FindMotifs(Map<String,String> ProteinNameToSequences){
		Map<String,Set<Integer>> ProteinsToMotifs = new HashMap<String, Set<Integer>>();
		for(Entry<String, String> e: ProteinNameToSequences.entrySet()){
			String sequence = e.getValue();
			Set<Integer> indices = new HashSet<Integer>(); 
			for(int i = 0; i <= (sequence.length()-windowSize); i++) { //"slide" through sequence in windows
				//D, E, S or K, with at least one K.  

				boolean k = false;
				double Dcount = 0.0;
				double Ecount = 0.0;
				double Scount = 0.0;
				double Kcount = 0.0;
				
				for(int j = i; j < i+windowSize; j++) { //count frequency of DESK amino acids
					if(sequence.charAt(j)=='D') {
						Dcount++;
					}
					if(sequence.charAt(j)=='E') {
						Ecount++;
					}
					if(sequence.charAt(j)=='S') {
						Scount++;
					}
					if(sequence.charAt(j)=='K') {
						Kcount++;
						k = true;
					}
				}
				double count = Dcount+Ecount+Scount+Kcount;
				if(k==true && count/windowSize >= aminoAcidRatio) { //save windows with aa ratios that pass threshold
					for(int j = i; j < i+windowSize; j++) {
						indices.add(j);
					}
					
					double newWindowSize = windowSize;
					int start = i;
					int end = i+(int)windowSize;
					while(count/newWindowSize>= aminoAcidRatio){ //expand the window (+1 in either direction at each step) to find more DESK residues, until aa ratio is no longer satisfied
						for(int j = start; j < end; j++) {
							indices.add(j);
						}
						if((start-1)<0) {
							start = 0;
						}
						else {
							start = start-1;
						}
						
						if((end+1) >=sequence.length()) {
							end = sequence.length()-1;
						}
						else {
							end = end+1;
						}
						newWindowSize = newWindowSize+2;
						
						Dcount = 0.0;
						Ecount = 0.0;
						Scount = 0.0;
						Kcount = 0.0;
							
						for(int j = start; j < end; j++) {
							if(sequence.charAt(j)=='D') {
								Dcount++;
							}
							if(sequence.charAt(j)=='E') {
								Ecount++;
							}
							if(sequence.charAt(j)=='S') {
								Scount++;
							}
							if(sequence.charAt(j)=='K') {
								Kcount++;
							}
						}
						count = Dcount+Ecount+Scount+Kcount;
					}
				}
			}
			ProteinsToMotifs.put(e.getKey(), indices); //save all DESK-heavy residue indices to the relevant protein
		}
		return ProteinsToMotifs;
	}
	
	/*
	 * Function reads a fasta file of protein sequences and store them in a HashMap
	 */
	private static Map<String,String> readProteinSequences(String fileName){
		Map<String,String> ProteinNameToSequences = new HashMap<String,String>();
		try {
			BufferedReader in = new BufferedReader(new FileReader(new File(fileName)));
			String s = in.readLine();
			while(s!=null){ //read from each > to next > (i.e. protein sequence across multiple lines) and save each seq to protein name
				String protein = s.split(">")[1];
				s = in.readLine();
				String seq = "";
				
				while(s.charAt(0)!='>'){
					seq = seq+s;
					s = in.readLine();
					if(s == null){
						break;
					}
				}
				ProteinNameToSequences.put(protein, seq);
			}
			in.close();
		}catch(Exception e) {
			e.printStackTrace();
		}
		return ProteinNameToSequences;
	}
}
