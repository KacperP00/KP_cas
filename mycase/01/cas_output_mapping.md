# Mapa plików wynikowych solvera CAS (1D Euler-Lagrange)

Dokument ten służy jako referencyjna mapa zmiennych dla plików wynikowych generowanych przez solver. Ponieważ większość plików `.out` ma nagłówki ukryte za znakiem komentarza `#`, a wartości przestrzenne są **znormalizowane (bezwymiarowe)**, poniższe zestawienie tłumaczy fizyczne i matematyczne znaczenie każdej kolumny.

---

## 1. Globalne parametry w czasie (`spray_PL.out` / `spray_water_PL.out`)
Plik ten opisuje ewolucję całej strugi w funkcji czasu. Otwiera się go w celu wyznaczenia głównych makroskopowych parametrów wtrysku.

| Kolumna | Zmienna | Fizyczne znaczenie i interpretacja |
| :--- | :--- | :--- |
| **1** | `step` | Numer kroku całkowania (iteracji) w pętli głównej solvera. |
| **2** | `time` | Aktualny czas symulacji [ms]. |
| **3** | `LPL` | **Liquid Penetration Length (Zasięg penetracji cieczy) [mm].** Maksymalna odległość od wylotu dyszy, na jakiej fizycznie istnieje faza ciekła (zanim całkowicie odparuje). |
| **4** | `VPL` | **Vapor Penetration Length (Zasięg chmury pary) [mm].** Maksymalna odległość, na jaką dotarł odparowany gaz z paliwa. W warunkach *flash-boilingu* LPL i VPL często się pokrywają. |
| **5** | `U_tip` | **Prędkość czoła strugi [m/s].** Prędkość, z jaką poruszają się krople znajdujące się na samym końcu zasięgu LPL. |
| **6** | `D_tip` | Średnica kropel (często Sautera - SMD) na samym czole strugi. |
| **7** | `Chi_st` | Stechiometryczna szybkość dyssypacji skalarnej. Parametr ważny wyłącznie przy włączonym module spalania (wskazuje na intensywność mieszania w strefie reakcji). |

---

## 2. Profile przestrzenne strugi (`data_*.out_XXXXXX`)
Główny plik diagnostyczny wyrzucany w określonych interwałach czasowych (np. `003267`). Opisuje stan strugi wzdłuż osi dyszy (1D). 
**UWAGA:** Zmienne dynamiczne i termodynamiczne są tu **bezwymiarowe**. Aby uzyskać ich fizyczne wartości (w m/s, K, mm), należy je przemnożyć przez wartości referencyjne z pliku `simulation_param.out`.

| Kolumna | Zmienna | Fizyczne znaczenie i interpretacja |
| :--- | :--- | :--- |
| **1** | `z` | Pozycja na osi wtrysku (odległość od dyszy) [mm]. |
| **2** | `rho` | Gęstość mieszaniny gazowej (powietrze + para paliwa) znormalizowana gęstością referencyjną. Rośnie drastycznie przy gwałtownym parowaniu. |
| **3** | `Y_l` | **Ułamek masy cieczy (Liquid mass fraction).** Od 1.0 (czysta ciecz) do 0.0 (brak cieczy). Wskazuje, ile procent masy w danej komórce siatki to niesparowane krople. |
| **4** | `Y_v` | **Ułamek masy pary (Vapor mass fraction).** Wskazuje, ile procent masy w komórce stanowi odparowane paliwo (gazowe). |
| **5** | `Y_a` | Ułamek masy gazu otoczenia (Ambient mass fraction). Wskazuje na stopień porwania powietrza przez strugę (entrainment). |
| **6** | `Y_g` | Całkowity ułamek fazy gazowej (`Y_g = Y_v + Y_a`). |
| **7** | `u_l` | **Bezwymiarowa prędkość kropel.** Należy pomnożyć przez `Injection velocity` z pliku parametrycznego, aby uzyskać wynik w [m/s]. Pokazuje utratę pędu pod wpływem oporu (drag). |
| **8** | `u_g` | **Bezwymiarowa prędkość gazu.** Prędkość, z jaką powietrze i para są porywane i wleczone przez uderzające krople. |
| **9** | `dm` | **Bezwymiarowa średnica kropel.** Należy pomnożyć przez `Effective jet diameter`, aby uzyskać wymiar w metrach. Kluczowy wskaźnik rozpadu strugi. |
| **13** | `Td` | **Bezwymiarowa temperatura kropel.** Należy pomnożyć przez początkową `Fuel temperature`. W warunkach *flash-boilingu* spada wzdłuż osi (utrata ciepła utajonego). |
| **14** | `Tg` | **Bezwymiarowa temperatura gazu.** Należy pomnożyć przez `Ambient temperature`. Pokazuje, jak gaz chłodzi się/nagrzewa od wtryskiwanych kropel. |
| **15** | `b` | Szerokość strugi (promień stożka) w danej odległości `z`. |

---

## 3. Analiza członów źródłowych - Fizyka Solver'a (`omega.out_XXXXXX`)
Pliki te zawierają surowe wartości tzw. *Source Terms* – równań przekazujących masę, pęd i energię pomiędzy kropelkami Lagrange'a a siatką Eulera. Służą do debugowania mechaniki modeli (np. dlaczego struga zwalnia lub paruje).

| Kolumna | Zmienna | Fizyczne znaczenie i interpretacja |
| :--- | :--- | :--- |
| **1** | `z` | Pozycja osiowa [mm]. |
| **2** | `omega_ent` | **Entrainment rate.** Szybkość, z jaką uwięzione wokół strugi powietrze otoczenia jest zasysane do wnętrza chmury kropel. |
| **3** | `omega_drag` | **Transfer pędu (Drag source term).** Siła hamująca wywierana przez gaz na krople. Ujemna wartość oznacza wytracanie prędkości przez ciecz i przekazywanie jej gazowi. |
| **4** | `omega_bre1` | Siła rozpadu wtórnego (aerodynamicznego KH-RT). Definiuje tempo zmniejszania się średnicy kropel. |
| **5** | `omega_bre2` | Siła rozpadu termicznego (*Flash-boiling*). Kluczowa zmienna w symulacjach przegrzanej cieczy. Jeśli jest równa zeru, zjawisko wrzenia nie niszczy kropel. |
| **6** | `omega_vap` | **Transfer masy (Vaporization rate).** Ilość odparowującego paliwa na jednostkę czasu. Dodawana do siatki gazowej, odejmowana z kropel. |
| **7** | `omega_T` | **Transfer ciepła.** Szybkość wymiany energii wewnętrznej. Pokazuje chłodzenie kropel w wyniku wymuszonego odparowania lub różnicy temperatur z otoczeniem. |

---

## 4. Dodatkowe pliki diagnostyczne

### `dsd.out_XXXXXX` oraz `bpdf.out_XXXXXX`
* **Przeznaczenie:** Reprezentacja funkcji rozkładu wielkości kropel (DSD - Droplet Size Distribution).
* **Struktura:** Brak nagłówków. Pierwsza kolumna to odległość od dyszy (`z`). Każda kolejna kolumna reprezentuje prawdopodobieństwo ułamkowe (lub gęstość prawdopodobieństwa) wystąpienia kropli o określonej średnicy wewnątrz rozpatrywanego przedziału (binu) klas wielkości.

### `simulation_param.out`
* **Przeznaczenie:** Plik generowany podczas pierwszego kroku `000000`. 
* **Znaczenie:** Zawiera kluczowe liczby bezwymiarowe (Liczba Reynoldsa, Webera) dla wtrysku wyliczone przez model wypływu z dyszy (Nozzle Flow Model) oraz **fizyczne wartości odniesienia** (Injection velocity, Effective jet diameter) niezbędne do zdekodowania znormalizowanych wyników z plików `data_*.out`.

### `roi_water.txt` (oraz `roi.txt`)
* **Przeznaczenie:** Wejściowy warunek brzegowy (Rate of Injection) wprowadzony przez użytkownika (lub wygenerowany z krzywej).
* **Struktura:** Kolumna 1 to czas, kolumna 2 to profil masowego natężenia przepływu lub znormalizowany mnożnik prędkości wypływu. Zmiany w tym pliku modelują kształt iglicy wtryskiwacza (otwieranie/zamykanie).
