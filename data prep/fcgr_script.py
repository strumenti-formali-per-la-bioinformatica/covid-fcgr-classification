
import sys
import os
from Bio import SeqIO
from complexcgr import FCGR
import matplotlib.pyplot as plt

# Si assicura che almeno un argomento sia stato passato (la directory di output)
if len(sys.argv) < 2:
    print("Usage: python script.py <directory_path>")
    sys.exit(1)

# La directory di output viene passata come primo argomento dallo script bash
directory_path = sys.argv[1]

# Inizializza FCGR
fcgr = FCGR(k=5) #numero di k-mers pari a 5

# Itera su tutti i file FASTA nella directory specificata
for filename in os.listdir(directory_path):
    if filename.startswith('chr') and filename.endswith('.fasta'):  # Controlla l'estensione e il prefisso del file
        file_path = os.path.join(directory_path, filename)
        print(f"{file_path} is being processed...")
        
        # Legge la sequenza dal file FASTA
        record = next(SeqIO.parse(file_path, 'fasta'))
        print(f"{record.id}: {len(record.seq)} nucleotides.")
        
        seq = str(record.seq)
        print(f"{seq[:200]}...")
        
        # Calcola FCGR per la sequenza
        chaos = fcgr(seq)
        
        # Genera e salva l'immagine
        fcgr.plot(chaos)
        image_path = os.path.join(directory_path, filename.replace('.fasta', '.jpeg'))
        fcgr.save_img(chaos, path=image_path)
        
        print("Immagine generata e salvata con successo.")

print("Totali immagini generate e salvate con successo.")
