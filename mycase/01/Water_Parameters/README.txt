Strategia na zdefiniowane parametrów wody do modelu CAS:

1. Skorzystać z pliku fuel_template.dat w folderze propertytables
2. Uzupełnić parametry wody w nieco szerszych od testowych zakresach temperatur i ciśnień (300 - 450K) oraz (1 - 8 bar). Wstępnie przygotowano skrypt Input_Liquid.py używając ct.water dla cieczy. "Suchą" parę wodną można przeliczyć przy użyciu h2o2.yaml zakładając gaz doskonały lub model Redlich-Kwong lub używając CoolProp (rekomendowane). Wstępnie napisano skrypty Input_Vapour_Cantera_IG oraz Input_Vapour_CoolProp do tego ale wymagają one dopracowania.
3. Zrobić post-processing danych. Parametry pary wodnej uzyskane różnymi metodami warto ze sobą porównać i obliczyć różnice.
4. W pliku fuel_template.dat zapisać parametry wody a następnie parametry "suchej" pary wodnej uzyskane w pierwszej kolejnośi z CoolProp.
5. Zmienić nazwę pliku na np. water.dat

Ad 2. 
Skrypt liczący parametry cieczy jest gotowy. Można dodać warunek kończący pętlę T gdy gęstość spadnie poniżej np. 100 [kg/m3].
* ct.water liczy parametry pary wodnej