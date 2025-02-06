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

public class MotifFinder {
	
	public static final double windowSize = 20;
	
	public static double aminoAcidRatio;
	public static void main(String[] args){
		aminoAcidRatio = Double.parseDouble(args[1]);
		String inputfileName = args[0];
		Map<String,String> ProteinNameToSequences = readProteinSequences(inputfileName);
		Map<String,Set<Integer>> ProteinsToMotifs = FindMotifs(ProteinNameToSequences);
		PrintSeqsWithMotifs(ProteinsToMotifs, ProteinNameToSequences, inputfileName);
		
	}
	
	/*
	 * function printing the sequences with the motifs identified
	 */
	
	private static void PrintSeqsWithMotifs(Map<String,Set<Integer>> ProteinsToMotifs,Map<String,String> ProteinNameToSequences, String inputFileName) {
		String outputFileName = inputFileName.substring(0,inputFileName.length()-6) + "_Result_Ratio_"+aminoAcidRatio+".txt";
		try {
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputFileName)));
			for(Entry<String, Set<Integer>> e: ProteinsToMotifs.entrySet()){
				if(e.getValue().size()>0) {
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
	}
	
	/*
	 * Function that finds the motifs
	 */
	private static Map<String,Set<Integer>> FindMotifs(Map<String,String> ProteinNameToSequences){
		Map<String,Set<Integer>> ProteinsToMotifs = new HashMap<String, Set<Integer>>();
		for(Entry<String, String> e: ProteinNameToSequences.entrySet()){
			String sequence = e.getValue();
			Set<Integer> indices = new HashSet<Integer>();
			for(int i = 0; i <= (sequence.length()-windowSize); i++) {
				//D, E, S or K, with at least one K.  

				boolean k = false;
				double Dcount = 0.0;
				double Ecount = 0.0;
				double Scount = 0.0;
				double Kcount = 0.0;
				
				for(int j = i; j < i+windowSize; j++) {
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
				if(k==true && count/windowSize >= aminoAcidRatio) {
					for(int j = i; j < i+windowSize; j++) {
						indices.add(j);
					}
				}
			}
			ProteinsToMotifs.put(e.getKey(), indices);
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
			while(s!=null){
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
