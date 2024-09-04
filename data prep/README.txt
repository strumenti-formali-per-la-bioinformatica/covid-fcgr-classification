Lo script da lanciare per ottenere le rappresentazioni FCGR.
In particolare, abbiamo seguito questo schema:
1. Siamo partiti dai file vcf ;
2. Ciascun vcf è stato ordinato  e successivamente è ottenuto l’indice;
3. In seguito, abbiamo applicato la funzione di FastaAlternateGenerate per ottenere il file fasta per ciascun paziente;
4. Dal file fasta ottenuto sono stati estratti tanti file fasta;
5. Per ciascun file fasta per paziente si è applicata la metodologia necessaria all’ottenimento della rappresentazione fcgr.

Lo script è stato creato per funzionare in modo sequenziale rispetto agli step precedentemente spiegati.
 In particolare, i punti  1 - 4 vengono eseguiti da un file bash. L’ultimo punto ,invece ,  è stato creato tramite un file Python che viene eseguito e richiamato dal file bash in maniera automatica una volta eseguiti gli step precedenti.

Nel testare lo script vengono create le seguenti cartelle per immagazzinare gli output:
 
New_vcf_sorted: contiene i file vcf ordinati e gli indici;
Output_fasta: contiene tante cartelle quanti sono i pazienti , ciascuna di esse conserva il file fasta completo, quello relativo ai cromosomi e le immagini fcgr.

Per quel che concerne il valore k stabilito per il k-mers  è stato fissato a 7. 