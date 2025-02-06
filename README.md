# PASKMotifFinder

A PASK Motif as defined by [Bentley-DeSousa et al.](https://doi.org/10.1016/j.celrep.2018.02.104) is a protein sequence of 20 amino acid length that is made up of at least 75% aspartate (D), glutamate (E), serine (S), and lysine (K) residues and contains at least one lysine. Such motifs are a likely site of lysine polyphosphorylation: the non-covalent binding of inorganic polyphosphate (PolyP) chains on a lysine residue. These post-translational modifications have been linked to stress responses and cellular regulation in bacteria and yeast, and many PASK-containing proteins share homologs in higher eukaryotes. Identifying PASK motifs in humans and other organisms is a critical step for investigating the biological roles of polyphosphorylation and uncovering conserved regulatory mechanisms that can be exploited for therapeutic strategies. 
    
The PASKMotifFinder is a Java tool for identifying PASK motif regions in protein sequences. The program reads protein sequences from a FASTA file and uses a 20 amino acid size sliding window to find sequence regions where DESK amino acids make up 75% or more of the sequence and contain at least one K residue. With PolyPhosphoFinder_Expand, the motif region can additionally be dynamically expanded to include adjacent residues until the ratio of DESK amino acids no longer satisfies the threshold. Both programs output a text file of the protein sequences where a PASK motif was found, with the motif region highlighted. 
    
## Dependencies
    
- Java (version 6 or greater) 

    
## PolyPhosphoFinder vs PolyPhosphoFinder_Expand
    
Both tools find motifs of 20 amino acid length that have enrichments of DESK amino acids (i.e. where the ratio of DESK amino acids to motif length is greater than or equal to 0.75).
    
PolyPhosphoFinder_Expand will additionally use each motif as a seed and dynamically expand the motif length by one amino acid in both directions. At each expansion step, the tool will check whether the composition of DESK amino acids in the new window still satisfies the ratio threshold. When the DESK ratio threshold is no longer satisfied, the longest acceptable motif will be saved and reported.

## Running PolyPhosphoFinder

0. Ensure you have Java installed first. 
1. Save `PolyPhosphoFinder.jar` to your preferred location (e.g. a folder in “Documents”).
    - You can follow the same steps for `PolyPhospoFinder_Expand.jar` . Be sure to use the correct .jar file name when launching the tool in step 5. 
2. Prepare and save your input files (described below) in the same location. 
3. Open a terminal window (Mac/Linux) or command prompt (Windows).
4. Navigate to the location with your files.
   - Mac example: `cd "/Users/username/Documents/MyProjectFolder/"`
   - Windows example: `cd "C:\Users\username\Documents\MyProjectFolder\"`
6. Run the program with the following arguments:

```bash
Command: 
java -jar PolyPhosphoFinder.jar [INPUT_FILE] [AMINO_ACID_RATIO]
    
Variables:
  * [INPUT_FILE] (required) = the path and filename of the FASTA file
  * [AMINO_ACID_RATIO] (required) = the minimum ratio of DESK amino acids needed in the motif (equal to or greater than)
```
        
Example commands:
    
```bash
#In a Mac terminal window, first navigate to the directory with the .jar file and the FASTA file
cd "/Users/myusername/Documents/MyPaskProject1/"
    
#Launch the motif finder program, looking for regions in the E. coli proteome where the composition ratio of DESK amino acids is equal to or greater than 0.75
java -jar PolyPhosphoFinder.jar uniprot-proteome_UP000000625.fasta 0.75
    
#The output file will be in "/Users/myusername/Documents/MyPaskProject1/" named "uniprot-proteome_UP000000625_Result_Ratio_0.75.txt"
```

```Shell
#In a Windows console window, first navigate to the directory with the .jar file and the FASTA file
cd "C:\Users\myusername\Downloads"
    
#Launch the expanded motif finder program, looking for regions in the E. coli proteome where the composition ratio of DESK amino acids is equal to or greater than 0.80
java -jar PolyPhosphoFinder_Expand.jar uniprot-proteome_UP000000625.fasta 0.8
    
#The output file will be in "C:\Users\myusername\Downloads" named "uniprot-proteome_UP000000625_Result_Ratio_0.8.txt"
```
    
## Input

1. FASTA File
    - Please provide a FASTA file to search. This can be a [reference proteome](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/) or any collection of protein sequences of interest, as long as you follow the FASTA file format.
    - The FASTA file must have protein names or descriptions on one line, beginning with a greater-than character (”>”). These names are to be followed by the protein sequence on the next line. The protein sequence can be all on one line, or on a series of lines.
        
    
This is an acceptable FASTA format:
    
```
>sp|Q46821|UACT_ECOLI Uric acid transporter UacT OS=Escherichia coli (strain K12) OX=83333 GN=uacT PE=1 SV=2
MSAIDSQLPSSSGQDRPTDEVDRILSPGKLIILGLQHVLVMYAGAVAVPLMIGDRLGLSK
EAIAMLISSDLFCCGIVTLLQCIGIGRFMGIRLPVIMSVTFAAVTPMIAIGMNPDIGLLG
QPLLHSGIMLATLSAVVLNVFFNGYQHHADLVKESVSDKDLKVRTVRMWLLMRKLKKNEH
GE
    
>sp|P77444|SUFS_ECOLI Cysteine desulfurase OS=Escherichia coli (strain K12) OX=83333 GN=sufS PE=1 SV=1
MIFSVDKVRADFPVLSREVNGLPLAYLDSAASAQKPSQVIDAEAEFYRHGYAAVHRGIHT
ALLQEMPPWEGGGSMIATVSLSEGTTWTKAPWRFEAGTPNTGGIIGLGAALEYVSALGLN
NIAEYEQNLMHYALSQLESVPDLTLYGPQNRLGVIAFNLGKHHAYDVGSFLDNYGIAVRT
GHHCAMPLMAYYNVPAMCRASLAMYNTHEEVDRLVTGLQRIHRLLG
```
    
This is also an acceptable FASTA format:
    
```
>Protein1
MNLHEYQAKQLFARYGLPAPVGYACTTPREAEEAASKIGAGPWVVKCQVHAGGRGKAGGVKVVNSKEDIRAFAENWLGKRLVTYQTDANGQPVNQILVEAATDIAKELYLGAVVDRSSRR
    
>This is Protein B, with a different format
FALPKERLWVTVYESDDEAYEIWEKEVGIPRERIIRIGDNKGAPYASDNFWQMGDTGPCG
PCTEIFYDHGDHIWGGPPGSPEEDGDRYIEIWNIVFMQFNRQADGTMEPLPKPSVDTGMG
LERIAAVLQHVNSNYDIDLFRTLIQAVAKVTGATDLSNKSLRVIADHIRSCAFLIADGVM
    
>ProteinThree | Small description | ID####
MSQLTHINAAGEAHMVDVSAKAETVREARAEAFVTMRSETLAMIIDGRHHKGDVFATARIAGIQAAKRTWDLIPLCHPLM
```
    
## Output
The tool will output a text file named as `[INPUT_FILE_NAME]_Result_Ratio_[AMINO_ACID_RATIO].txt`, saved in the same directory as where you run the command. 
    
The output file will include the protein names and sequences where a PASK motif was found. The motif region will be highlighted with asterisks (”*”) above the corresponding amino acids.

```
sp|P0A8H6|YIHI_ECOLI Der GTPase-activating protein YihI OS=Escherichia coli (strain K12) OX=83333 GN=yihI PE=1 SV=1
	                                                  
0	MKPSSSNSRSKGHAKARRKTREELDQEARDRKRQKKRRGHAPGSRAAGGN
	                                                  
50	TTSGSKGQNAPKDPRIGSKTPIPLGVTEKVTKQHKPKSEKPMLSPQAELE
	                                         *********
100	LLETDERLDALLERLEAGETLSAEEQSWVDAKLDRIDELMQKLGLSYDDD
	***********        
150	EEEEEDEKQEDMMRLLRGN

sp|P77173|ZIPA_ECOLI Cell division protein ZipA OS=Escherichia coli (strain K12) OX=83333 GN=zipA PE=1 SV=3
	                                         *********
0	MMQDLRLILIIVGAIAIIALLVHGFWTSRKERSSMFRDRPLKRMKSKRDD
	***********                                       
50	DSYDEDVEDDEGVGEVRVHRVNHAPANAQEHEAARPSPQHQYQPPYASAQ
	                                                  
100	PRQPVQQPPEAQVPPQHAPHPAQPVQQPAYQPQPEQPLQQPVSPQVAPAP
	                                                  
150	QPVHSAPQPAQQAFQPAEPVAAPQPEPVAEPAPVMDKPKRKEAVIIMNVA
	                                                  
200	AHHGSELNGELLLNSIQQAGFIFGDMNIYHRHLSPDGSGPALFSLANMVK
	                                                  
250	PGTFDPEMKDFTTPGVTIFMQVPSYGDELQNFKLMLQSAQHIADEVGGVV
	                            
300	LDDQRRMMTPQKLREYQDIIREVKDANA
```
    
Troubleshooting possible issues:
- Empty output file → If an output file is created but empty, no motifs were found. You could try a less stringent amino acid ratio threshold.
- Motif region highlight alignment → If the asterisks that highlight the motif region are not aligning well with the protein sequence, please try changing to a monospaced font where all characters are the same width (e.g. Courier, Consolas, Lucida Console)
