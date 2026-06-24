# Next steps for the case

- wybrać co robić
- są skrypty liczące parametry wody i SUCHEJ pary wodnej 
- dokończyc input do casa + walidacja
- gasjet / CFD
- 3D / CFD z x-ray, czekanie na pliki - liczenie 
- bonus-pytanie: kod czy przegląd? KOD

- model na one-drive do wymiany plików od p. Chmielewskiego 
https://wutwaw-my.sharepoint.com/:f:/g/personal/01133018_pw_edu_pl/IgCj4Luc7GMNQJ7W8JrN_caPAWw-mb0nYJ00Isjrr_Ubmj4?e=IP1oxL

- nie dzielic się gasjet.docx, kim ma najdokłądniej (rozdział 2.2.2)
- w przeglądówce spojrzec na gasjet (Li et al)
- z SAE spróbować znaleźć publikacje w granicach rozsądku
- przejrzeć input_vapour_cantera_ig - zobaczyć czy coś można dodać
- 3. !!! obrobić dane - jak h2o2 wypada w trzech case'ach (ideal, real, coolprop)
- 4. !!! po tym jak będą wyniki, zacząc z cas od ideal gas
- 1. !!! dodać DiffusionCoefficient w skrypcie dla pary wodnej
- 2. !!! w tabeli dla vapor nie ma ciśnienia - znaleźć czy parametry mają być dla ciśnienia otoczenia czy płynu
- baseline to ideal gas
- czy skrypty vapor też wykrywają zmiany fazy?


!!! Zmiana: propertytables nie działają - trzeba wpisać wodę w kod (jak inne paliwa). Już nie kukułcze jajo, ale wpisujemy dokładnie wodę w kod:
"Myślę, że trzeba cofnąć się do źródła i w pliku pc_database.f90 zdefiniować parametry wody w podobny sposób jak zdefiniowane są pozostałe paliwa. Wtedy plik makefile doda pliki pc_database.o i pc_database.mod odpowiednio do katalogów obj i mod. Trzeba jeszcze zmodyfikować plik pc_func.f90. Po dopisaniu kolejnego case'a 'water', CAS nie wyrzuca błędu o braku paliwa w bazie. Pliki po zmianach zamieściłem tutaj"

Potem dodano:
"
Proszę zapoznać się ze zmienionym plikami. f90. Wystarczy, że wyszuka Pan "water" w plikach. Tylko w tych miejscach wprowadzałem zmiany 
 
Z rzeczy TO DO jest jeszcze dobranie współczynników A-D wielomianów opisujących parametey wody jak gęstość, lepkość itd. Trzeba jeszcze rozkminić definicję range 1 i range 2 
 "