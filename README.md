# AlphaFold3 modelio tikslumo pokytis didėjant baltymų konformacinei įvairovei

Šioje repozitorijoje pasiekiami visi skriptai apo ir holo prognozuotų baltymų struktūrų analizei.

Duomenų paruošimo skriptai yra `python\ scripts/` direktorijoje. Taip pat pasiekiami [Google Colab](https://colab.research.google.com/drive/16k_oZTqws-6uCptCGIRy0ykzQ-ltVHtL?usp=sharing). Statistikos skaičiavimo skriptas `Plots.R` pasiekiamas šioje repozitorijoje. Visos diagramos randamos `Stats/` direktorijoje. 

Duomenų failai yra `data/` direktorijoje. Didesni duomenų failai pasiekiami [Google Drive](https://drive.google.com/drive/folders/1LV_uK2zhjc5m4qGG6RYtLP5EhHihjW4A?usp=sharing).

## Vizualizacija

  - RMSD reikšmių tankio grafikai.
  - Apo ir holo struktūrų palyginimo sklaidos grafikai.
  - Baltymų, padalintų į šeimas ir lanksčias grupes, histogramos.
  - Konformacinės įvairovės įtaka plDDT įverčiui.
  - RMSF įverčių palyginimas su plDDT įverčiais.

## R kalbos reikalavimai

- R (≥ 4.0)
- R paketai:
  - `ggplot2`
  - `dplyr`
  - `tidyr`
  - `readr`
  - `ggpubr`
  - `gridExtra`

Kaip instaliuoti paketus:
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "ggpubr", "gridExtra"))
```

## Naudojimas

1. Repozitorijos klonavimas:
   ```bash
   git clone https://github.com/Komceks/Impact-of-protein-conformational-diversity-on-AlphaFold3-predictions.git
   ```

2. Paleisti `Plots.R` RStudio programoje ar kitoje R aplinkoje.

3. Nusistatykite darbinę direktoriją.

5. Paleidus skriptą grafikai išsaugomi `Stats/` direktorijoje.

## Gauti grafikai

### **Figure 0**: Tankio grafikas
- RMSD pasiskirstymas tarp APO ir HOLO struktūrų.

### **Figure 1A & 1B**: Sklaidos grafikai
- APO struktūrų panašumas į eksperimentines APO ir HOLO struktūras sklaidos grafikai (pilnas ir priartintas).

### **Figure 2A & 2B**: Sklaidos grafikai
- HOLO struktūrų panašumas į eksperimentines APO ir HOLO struktūras sklaidos grafikai (pilnas ir priartintas).

### **Figure 4A & 4B**: plDDT pokyčio, didėjant RMSD, sklaidos grafikas

### **Figure 5A-D**: Sklaidos grafikai 
- Maksimalių ir minimalių RMSD reikšmių tarp eksperimentinių ir prognozuotų struktūrų pokytis didėjant konformacinei įvairovei.

### **Figure 6A & 6B**: Histogramos
- RMSD pokytis baltymų šeimose (homogeniškose ir heterogeniškose).

### **Figure 7A & 7B**: Histogramos
- RMSD pokytis baltymų lankstumo grupėse (lankstūs ir nelankstūs).

### **Figure 8A & 8B**: Histogramos
- RMSF pokytis didėjant plDDT įverčiui baltymuose.
