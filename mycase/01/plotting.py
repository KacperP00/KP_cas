import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

# Tworzy docelowy katalog na wykresy
os.makedirs('plots', exist_ok=True)

# Pobiera listę wszystkich plików wyjściowych
files = glob.glob('*spray*')

for file in files:
    # Pomija pliki z metadanymi i parametrami stałymi
    if 'simulation_param' in file:
        continue

    try:
        # Ekstrakcja nagłówków z zakomentowanej linii (specyfika plików CAS)
        with open(file, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('#'):
                headers = first_line.replace('#', '').split()
            else:
                headers = None

        # Wczytanie danych numerycznych
        df = pd.read_csv(file, delim_whitespace=True, comment='#', header=None)
        
        # Przerywa jeśli plik jest pusty
        if df.empty or df.shape[1] < 2:
            continue

        # W spray_PL.out kolumna 1 to czas (X). W plikach data/omega kolumna 0 to oś Z (X).
        is_spray_pl = 'spray_water_PL' in file
        x_col_idx = 1 if is_spray_pl else 0
        x_data = df.iloc[:, x_col_idx]
        x_label = headers[x_col_idx] if headers and len(headers) > x_col_idx else f'Kolumna {x_col_idx}'

        # Generuje wykres dla każdej pozostałej zmiennej w pliku
        for y_col_idx in range(df.shape[1]):
            # Pomija oś X oraz kolumnę 0 dla spray_PL (iteracja solvera)
            if y_col_idx == x_col_idx or (is_spray_pl and y_col_idx == 0):
                continue

            y_label = headers[y_col_idx] if headers and len(headers) > y_col_idx else f'Kolumna {y_col_idx}'

            plt.figure(figsize=(10, 6))
            plt.plot(x_data, df.iloc[:, y_col_idx], label=y_label)
            plt.title(f'{file}: {y_label} vs {x_label}')
            plt.xlabel(x_label)
            plt.ylabel(y_label)
            plt.grid(True)
            
            # Formatuje nazwę pliku wynikowego i zapisuje wykres
            safe_y_label = "".join([c for c in y_label if c.isalnum() or c in ('_', '-')])
            save_path = os.path.join('plots', f"{file}_{safe_y_label}.png")
            plt.savefig(save_path)
            plt.close()

    except Exception as e:
        print(f"Błąd przetwarzania pliku {file}: {e}")

print("Generowanie wykresów zakończone.")