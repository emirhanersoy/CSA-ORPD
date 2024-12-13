import numpy as np
import math

class IEEE14BaraSistemi:
    def __init__(self):
        # Bara Parametreleri
        self.bara_sayisi = 14
        self.pv_baralari = [1, 2, 3, 4, 5]  # PV bara indeksleri
        self.slack_bara = 0  # Slack bara indeksi

        # Gerilim Sınırları
        self.v_min_pq = 0.95  # PQ baraları için min gerilim
        self.v_max_pq = 1.05  # PQ baraları için max gerilim
        self.v_min_pv = 1.0   # PV baraları için min gerilim
        self.v_max_pv = 1.1   # PV baraları için max gerilim

        # Düzeltilmiş Hat Parametreleri
        self.hat_verileri = {
            'from_bus': [0, 0, 1, 1, 1, 2, 3, 5, 5, 5, 8, 8, 9, 11, 12],
            'to_bus':   [1, 4, 2, 3, 4, 3, 4, 10, 11, 12, 9, 13, 10, 12, 13],
            'r': [0.01938, 0.05403, 0.04699, 0.06701, 0.01335, 0.01487, 0.01524, 
                  0.03181, 0.06414, 0.06113, 0.08130, 0.15010, 0.01480, 0.11850, 0.11200],
            'x': [0.05917, 0.22304, 0.19797, 0.17103, 0.04211, 0.04234, 0.05075,
                  0.08437, 0.15544, 0.17632, 0.17070, 0.41070, 0.06470, 0.39250, 0.37200],
            'b': [0] * 15  # Susceptance değerleri
        }

        # Reaktif Güç Sınırları
        self.reaktif_guc_sinir = {
            'Qmin': [-0.4] * 5,  # PV baralar için min reaktif güç
            'Qmax': [0.5] * 5    # PV baralar için max reaktif güç
        }

        # Transformatör Kademe Ayarları
        self.trafo_min = 0.9
        self.trafo_max = 1.1
        self.trafo_adim = 0.0125

def aktif_guc_kaybi_hesapla(sistem, v, theta):
    """Detaylı aktif güç kaybı hesaplaması."""
    Ploss = 0
    for i in range(len(sistem.hat_verileri['from_bus'])):
        from_bus = sistem.hat_verileri['from_bus'][i]
        to_bus = sistem.hat_verileri['to_bus'][i]
        r = sistem.hat_verileri['r'][i]
        x = sistem.hat_verileri['x'][i]
        
        # Kompleks empedans hesabı
        z = r + 1j * x
        
        # Gerilim ve açı
        v_from = v[from_bus]
        v_to = v[to_bus]
        theta_from = theta[from_bus]
        theta_to = theta[to_bus]
        
        # Kompleks akım hesabı
        I_mag = v_from / np.abs(z) * np.sin(theta_from - theta_to + np.angle(z))
        
        # Aktif güç kaybı
        Ploss += I_mag**2 * r
    
    return Ploss

def kisitlama_kontrol(sistem, v, theta, trafo_kademe, reaktif_guc):
    """Gelişmiş kısıtlama kontrolü."""
    # Gerilim sınırları kontrolü
    for i in range(sistem.bara_sayisi):
        if i in sistem.pv_baralari:
            if v[i] < sistem.v_min_pv or v[i] > sistem.v_max_pv:
                return False
        else:
            if v[i] < sistem.v_min_pq or v[i] > sistem.v_max_pq:
                return False
    
    # Transformatör kademe sınırları kontrolü
    for kademe in trafo_kademe:
        if kademe < sistem.trafo_min or kademe > sistem.trafo_max:
            return False
    
    # Reaktif güç sınırları kontrolü
    for i, q in enumerate(reaktif_guc):
        if q < sistem.reaktif_guc_sinir['Qmin'][i] or q > sistem.reaktif_guc_sinir['Qmax'][i]:
            return False
    
    return True

def csa_orpd(sistem, kristal_sayisi=10, maksimum_iterasyon=100):
    """Crystal Structure Algoritması ile ORPD çözümü"""
    # Boyutları tanımla
    pv_bara_sayisi = len(sistem.pv_baralari)
    boyut = sistem.bara_sayisi + len(sistem.hat_verileri['from_bus']) + pv_bara_sayisi + pv_bara_sayisi
    
    # Alt ve üst sınırları tanımla
    alt_sinir = np.concatenate([
        np.full(sistem.bara_sayisi, sistem.v_min_pq),          # Gerilimler
        np.full(len(sistem.hat_verileri['from_bus']), -np.pi/2),  # Açılar
        np.full(pv_bara_sayisi, sistem.trafo_min),             # Transformatör kademeleri
        np.array(sistem.reaktif_guc_sinir['Qmin'])             # Reaktif güç
    ])
    
    ust_sinir = np.concatenate([
        np.full(sistem.bara_sayisi, sistem.v_max_pq),          # Gerilimler
        np.full(len(sistem.hat_verileri['from_bus']), np.pi/2),   # Açılar
        np.full(pv_bara_sayisi, sistem.trafo_max),             # Transformatör kademeleri
        np.array(sistem.reaktif_guc_sinir['Qmax'])             # Reaktif güç
    ])
    
    # Başlangıç kristal popülasyonu
    kristaller = np.random.uniform(alt_sinir, ust_sinir, (kristal_sayisi, boyut))
    
    # En iyi çözüm ve uygunluk
    en_iyi_uygunluk = float('inf')
    en_iyi_cozum = None
    
    for iterasyon in range(maksimum_iterasyon):
        # Her kristal için
        for i in range(kristal_sayisi):
            # Mevcut kristali çözüm olarak al
            x_mevcut = kristaller[i]
            
            # Çözüm bileşenlerini ayır
            v = x_mevcut[:sistem.bara_sayisi]
            theta = x_mevcut[sistem.bara_sayisi:sistem.bara_sayisi+len(sistem.hat_verileri['from_bus'])]
            trafo_kademe = x_mevcut[sistem.bara_sayisi+len(sistem.hat_verileri['from_bus']):
                                     sistem.bara_sayisi+len(sistem.hat_verileri['from_bus'])+pv_bara_sayisi]
            reaktif_guc = x_mevcut[-pv_bara_sayisi:]
            
            # Kısıtlamaları kontrol et
            if kisitlama_kontrol(sistem, v, theta, trafo_kademe, reaktif_guc):
                # Uygunluk hesaplama
                uygunluk = aktif_guc_kaybi_hesapla(sistem, v, theta)
                
                # En iyi çözümü güncelle
                if uygunluk < en_iyi_uygunluk:
                    en_iyi_uygunluk = uygunluk
                    en_iyi_cozum = x_mevcut
            
            # Kristal yapısını değiştirme (kristalizasyon)
            for j in range(boyut):
                # Rastgele kristalizasyon mekanizması
                if np.random.rand() < 0.5:
                    # Küçük rastgele değişim
                    x_mevcut[j] += np.random.normal(0, 0.1) * x_mevcut[j]
                else:
                    # Global en iyi çözüme doğru yönelme
                    if en_iyi_cozum is not None:
                        x_mevcut[j] += np.random.rand() * (en_iyi_cozum[j] - x_mevcut[j])
            
            # Sınırları kontrol et
            x_mevcut = np.clip(x_mevcut, alt_sinir, ust_sinir)
            
            # Güncellenmiş kristali kaydet
            kristaller[i] = x_mevcut
    
    return en_iyi_cozum, en_iyi_uygunluk

def csa_performans_testi(sistem, deneme_sayisi=30, kristal_sayisi=10, maksimum_iterasyon=100):
    """
    Crystal Structure Algoritmasının performansını test eden fonksiyon
    
    Parametreler:
    - sistem: IEEE14BaraSistemi nesnesi
    - deneme_sayisi: Bağımsız çalıştırma sayısı
    - kristal_sayisi: Her çalıştırmadaki kristal sayısı
    - maksimum_iterasyon: Her çalıştırmadaki maksimum iterasyon sayısı
    
    Returns:
    Performans istatistikleri dictionary'si
    """
    sonuclar = []
    tum_cozumler = []
    
    for _ in range(deneme_sayisi):
        # Rastgele tohum kullanarak deterministik olmayan sonuçlar elde et
        np.random.seed(None)
        
        # Crystal Structure algoritması ile çözüm
        en_iyi_cozum, en_iyi_uygunluk = csa_orpd(
            sistem, 
            kristal_sayisi=kristal_sayisi, 
            maksimum_iterasyon=maksimum_iterasyon
        )
        
        sonuclar.append(en_iyi_uygunluk)
        tum_cozumler.append(en_iyi_cozum)
    
    # Performans istatistikleri
    performans = {
        'ortalama': np.mean(sonuclar),
        'standart_sapma': np.std(sonuclar),
        'en_iyi': np.min(sonuclar),
        'en_kotu': np.max(sonuclar),
        'en_iyi_cozum': tum_cozumler[np.argmin(sonuclar)]
    }
    
    return performans

# Ana çalıştırma
if __name__ == "__main__":
    # IEEE 14-Bara Sistemini oluştur
    sistem = IEEE14BaraSistemi()
    
    # Performans testini çalıştır
    performans_sonuclari = csa_performans_testi(sistem)
    
    # Sonuçları detaylı yazdır
    print("\n--- Crystal Structure Algoritması Performans Sonuçları ---")
    for anahtar, deger in performans_sonuclari.items():
        print(f"{anahtar}: {deger}")