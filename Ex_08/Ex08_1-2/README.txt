Al momento il programma esegue una simulazione MC con i valori ottimali di mu e sigma e stampa nella cartella Temp/Average i valori medi di energia totale, cinetica, potenziale e l’istogramma in funzione del numero di blocchi.

Per effettuare il VMC è necessario modificare alcuni parametri nel file input.dat (impostare un valore positivo per nvar (es. 150) e porre print_blocks = 0).
Il VMC viene effettuato automaticamente a partire da 15 configurazioni diverse di mu e sigma riportate nel file start_param.txt (colonna 1:mu colonna 2:sigma). 
Il programma stampa i 15 cammini nello spazio dei parametri nella cartella Temp/Paths.