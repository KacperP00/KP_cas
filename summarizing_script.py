import os
import glob

def pobierz_parametry_odniesienia():
    # Pobiera wartosci referencyjne z pliku wynikowego
    params = {'d_ref': 1.0, 'u_ref': 1.0}
    if os.path.exists('simulation_param.out'):
        with open('simulation_param.out', 'r') as f:
            for line in f:
                if 'Effective jet diameter' in line:
                    params['d_ref'] = float(line.split(':')[1].strip())
                elif 'Injection velocity' in line:
                    params['u_ref'] = float(line.split(':')[1].strip())
    return params

def pobierz_temperatury():
    # Pobiera temperatury wejsciowe z pliku konfiguracyjnego
    t_fuel, t_amb = 293.15, 293.15
    for cas_file in ['cas.in', '01cas.in']:
        if os.path.exists(cas_file):
            with open(cas_file, 'r') as f:
                for line in f:
                    if 'Fuel temperature' in line:
                        t_fuel = float(line.split(':')[1].strip())
                    elif 'Ambient temperature' in line:
                        t_amb = float(line.split(':')[1].strip())
            break
    return t_fuel, t_amb

def generuj_zestawienie_globalne(d_ref, u_ref):
    # Skaluje i zapisuje parametry zalezne od czasu
    mapa_czasu = {}
    pliki_spray = glob.glob('*spray*PL.out')
    if not pliki_spray:
        return mapa_czasu
    
    plik_we = pliki_spray[0]
    plik_wy = 'Skalowane_Wyniki_Globalne.csv'
    
    with open(plik_we, 'r') as we, open(plik_wy, 'w') as wy:
        wy.write('Krok,Czas_[ms],Zasieg_cieczy_LPL_[mm],Zasieg_pary_VPL_[mm],Predkosc_czola_[m/s],Srednica_czola_SMD_[um]\n')
        for linia in we:
            if linia.startswith('#'): continue
            dane = linia.split()
            if len(dane) < 6: continue
            
            krok = dane[0]
            czas_ms = float(dane[1]) * 1000.0
            lpl_mm = float(dane[2]) * d_ref * 1000.0
            vpl_mm = float(dane[3]) * d_ref * 1000.0
            u_tip = float(dane[4]) * u_ref
            d_tip = float(dane[5]) * d_ref * 1e6
            
            mapa_czasu[krok] = czas_ms
            wy.write(f'{krok},{czas_ms:.5f},{lpl_mm:.4f},{vpl_mm:.4f},{u_tip:.4f},{d_tip:.4f}\n')
            
    return mapa_czasu

def generuj_zestawienie_przestrzenne(d_ref, u_ref, t_fuel, t_amb, mapa_czasu):
    # Llaczy wszystkie kroki przestrzenne w jeden plik i skaluje wartosci
    pliki_data = glob.glob('data*.out_*')
    if not pliki_data:
        return
        
    plik_wy = 'Skalowane_Wyniki_Przestrzenne.csv'
    naglowek = 'Krok,Czas_[ms],Pozycja_Z_[mm],Gestosc_wzgledna_[-],Ulamek_cieczy_Yl_[-],Ulamek_pary_Yv_[-],' \
               'Predkosc_cieczy_[m/s],Predkosc_gazu_[m/s],Srednica_kropel_SMD_[um],' \
               'Temp_kropel_[K],Temp_gazu_[K],Promien_strugi_b_[mm]\n'
               
    with open(plik_wy, 'w') as wy:
        wy.write(naglowek)
        for plik in sorted(pliki_data):
            krok = plik.split('_')[-1]
            czas_ms = mapa_czasu.get(krok, 0.0)
            
            with open(plik, 'r') as we:
                for linia in we:
                    if linia.startswith('#'): continue
                    dane = linia.split()
                    if len(dane) < 15: continue
                    
                    z_mm = float(dane[0]) * d_ref * 1000.0
                    rho = float(dane[1])
                    y_l = float(dane[2])
                    y_v = float(dane[3])
                    u_l = float(dane[6]) * u_ref
                    u_g = float(dane[7]) * u_ref
                    dm_um = float(dane[8]) * d_ref * 1e6
                    t_d = float(dane[12]) * t_fuel
                    t_g = float(dane[13]) * t_amb
                    b_mm = float(dane[14]) * d_ref * 1000.0
                    
                    wy.write(f'{krok},{czas_ms:.5f},{z_mm:.4f},{rho:.4f},{y_l:.6f},{y_v:.6f},'
                             f'{u_l:.4f},{u_g:.4f},{dm_um:.4f},{t_d:.2f},{t_g:.2f},{b_mm:.4f}\n')

def main():
    # Glowna funkcja sterujaca
    parametry = pobierz_parametry_odniesienia()
    t_fuel, t_amb = pobierz_temperatury()
    
    mapa_czasu = generuj_zestawienie_globalne(parametry['d_ref'], parametry['u_ref'])
    generuj_zestawienie_przestrzenne(parametry['d_ref'], parametry['u_ref'], t_fuel, t_amb, mapa_czasu)

if __name__ == '__main__':
    main()