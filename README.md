# Analiza Różnicowej Ekspresji Genów (DEG) i Adnotacja funkcjonalna dla Danych RNA-Seq z eksperymentu Pasilla

Ten dokument opisuje krok po kroku proces analizy danych RNA-Seq pochodzących z eksperymentu "Pasilla" (badanie wpływu knockdownu genu pasilla na komórki Drosophila melanogaster S2-DRSC), wykorzystując pakiety R: `DESeq2` do analizy różnicowej ekspresji genów (DEG) oraz `gprofiler2` do analizy wzbogacenia funkcjonalnego (Gene Ontology i ścieżki KEGG).

## 1. Wprowadzenie

Celem analizy jest identyfikacja genów, których poziom ekspresji zmienia się istotnie pomiędzy próbkami poddanymi działaniu siRNA przeciwko genowi *pasilla* (treated) a próbkami kontrolnymi (untreated). Następnie, dla zidentyfikowanej listy genów różnicowo wyeksponowanych (DEG), przeprowadzana jest analiza wzbogacenia funkcjonalnego, aby zrozumieć, w jakie procesy biologiczne, funkcje molekularne, komponenty komórkowe czy szlaki metaboliczne zaangażowane są te geny.

## 2. Wymagane Pakiety R

Przed rozpoczęciem analizy upewnij się, że masz zainstalowane wszystkie niezbędne pakiety R. Poniższy kod sprawdza obecność pakietów i instaluje je, jeśli są potrzebne. Wykorzystywany jest `BiocManager` do instalacji pakietów biokonduktorowych (`pasilla`, `DESeq2`) oraz standardowa funkcja `install.packages` dla pakietów z CRAN (`RColorBrewer`, `pheatmap`, `ggplot2`, `dplyr`, `gprofiler2`).

```r
# Sprawdzenie i instalacja BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Lista pakietów Bioconductor
bioc_packages <- c("pasilla", "DESeq2")
# Lista pakietów CRAN
cran_packages <- c("RColorBrewer", "pheatmap", "ggplot2", "dplyr", "gprofiler2")

# Instalacja pakietów Bioconductor, jeśli brakuje
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

# Instalacja pakietów CRAN, jeśli brakuje
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

# Ładowanie bibliotek
library("pasilla")
library("DESeq2")
library("RColorBrewer")
library("pheatmap")
library("ggplot2")
library("dplyr")
library("gprofiler2")
```
## 3. Wczytywanie i Przygotowanie Danych
Dane do analizy pochodzą z pakietu pasilla. Zawiera on:Macierz zliczeń: Plik pasilla_gene_counts.tsv z surowymi zliczeniami odczytów przypisanych do poszczególnych genów (identyfikatory FlyBase).Metadane: Plik pasilla_sample_annotation.csv opisujący próbki (warunek eksperymentalny, typ biblioteki itp.).
### 3.1. Wczytywanie Danych
   Wczytujemy oba pliki. Macierz zliczeń jest konwertowana do formatu matrix, a z metadanych wybieramy interesujące nas kolumny (condition i type) i konwertujemy je na typ factor.
   ```r
   # Ścieżki do plików danych
   pasCtsFile <- system.file("extdata", "pasilla_gene_counts.tsv", package="pasilla", mustWork=TRUE)
   pasAnnoFile <- system.file("extdata", "pasilla_sample_annotation.csv", package="pasilla", mustWork=TRUE)
   
   # Wczytanie macierzy zliczeń
           cts <- as.matrix(read.csv(pasCtsFile, sep="\t", row.names="gene_id"))

           # Wczytanie i przygotowanie metadanych
           coldata_raw <- read.csv(pasAnnoFile, row.names=1)
           coldata <- coldata_raw[, c("condition", "type")]
           coldata$condition <- factor(coldata$condition)
           coldata$type <- factor(coldata$type)
   ```
### 3.2. Dopasowanie Nazw Próbek
Kluczowe jest, aby nazwy kolumn w macierzy zliczeń (cts) odpowiadały nazwom wierszy w tabeli metadanych (coldata) i były w tej samej kolejności. W tym zbiorze danych konieczne jest usunięcie prefiksu "fb" z nazw            wierszy w coldata i uporządkowanie kolumn w cts.
```r
# Ujednolicenie nazw próbek
 rownames(coldata) <- sub("fb", "", rownames(coldata))

# Sprawdzenie zgodności nazw (ważne!)
all(rownames(coldata) %in% colnames(cts)) # Czy wszystkie nazwy z coldata są w cts?
all(rownames(coldata) == colnames(cts))   # Czy są w tej samej kolejności?

# Uporządkowanie kolumn macierzy zliczeń, jeśli kolejność nie jest zgodna
if (!all(rownames(coldata) == colnames(cts))) {
  cts <- cts[, rownames(coldata)]
}

# Ostateczne sprawdzenie zgodności kolejności
all(rownames(coldata) == colnames(cts))
```
## 4. Analiza Różnicowej Ekspresji z DESeq2

Pakiet `DESeq2` jest standardowym narzędziem bioinformatycznym służącym do identyfikacji genów różnicowo wyeksponowanych (Differentially Expressed Genes - DEG) w danych pochodzących z sekwencjonowania RNA (RNA-Seq). Opiera się on na modelowaniu surowych zliczeń odczytów dla każdego genu przy użyciu rozkładu ujemnego dwumianowego (Negative Binomial distribution). Rozkład ten dobrze opisuje dane zliczeniowe, które charakteryzują się nadmierną dyspersją (wariancja większa niż średnia), co jest typowe dla danych RNA-Seq.

### 4.1. Tworzenie Obiektu `DESeqDataSet`

Pierwszym krokiem jest stworzenie obiektu klasy `DESeqDataSet`. Jest to centralna struktura danych w pakiecie `DESeq2`, która integruje:
* **Macierz zliczeń (`countData`):** Surowe, nieprzetworzone zliczenia odczytów dla każdego genu (wiersze) w każdej próbce (kolumny). Ważne jest, aby były to liczby całkowite.
* **Tabelę metadanych (`colData`):** Informacje opisujące próbki (kolumny macierzy zliczeń). Musi zawierać zmienne, które będą używane w analizie (np. warunek eksperymentalny, typ próbki, batch). Nazwy wierszy w `colData` muszą odpowiadać nazwom kolumn w `countData`.
* **Formułę projektową (`design`):** Formuła R określająca, które zmienne z `colData` mają być użyte do modelowania ekspresji genów. W najprostszym przypadku, jak tutaj (`~ condition`), modelujemy wpływ jednej zmiennej (warunku eksperymentalnego) na poziom ekspresji. Bardziej złożone formuły mogą uwzględniać dodatkowe zmienne lub interakcje.

```r
# Tworzenie obiektu DESeqDataSet
# countData: macierz surowych zliczeń (cts)
# colData: tabela metadanych (coldata)
# design: formuła określająca model statystyczny (~ condition oznacza modelowanie wpływu zmiennej 'condition')
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)

# Wyświetlenie podsumowania obiektu (informacje o wymiarach, zmiennych w colData itp.)
print(dds)
```
### 4.2. Ustawienie Poziomu Referencyjnego
DESeq2 porównuje warunki alfabetycznie. Jeśli mamy warunki "treated" i "untreated", porównanie będzie "treated vs untreated". Jednak dla jasności i spójności interpretacji wyników (szczególnie log2 Fold Change), dobrą praktyką jest jawne zdefiniowanie, który poziom faktora (w tym przypadku condition) ma być poziomem referencyjnym (kontrolnym). Używamy funkcji relevel, aby ustawić "untreated" jako referencję. Dzięki temu log2FoldChange będzie oznaczał zmianę ekspresji w grupie "treated" względem grupy "untreated".
```r
# Ustawienie 'untreated' jako poziomu referencyjnego dla faktora 'condition'
dds$condition <- relevel(dds$condition, ref = "untreated")

# Sprawdzenie poziomów faktora (pierwszy powinien być teraz "untreated")
levels(dds$condition)
```
### 4.3. Uruchomienie Głównej Analizy DESeq
Funkcja DESeq() jest sercem pakietu i wykonuje kompleksową analizę w kilku krokach:
Normalizacja (Estimation of size factors): Oblicza tzw. "size factors" dla każdej próbki. Są to współczynniki normalizacyjne, które korygują różnice w całkowitej liczbie odczytów (głębokości sekwencjonowania) i składzie RNA między próbkami. DESeq2 używa metody mediany stosunków (median-of-ratios) do estymacji tych czynników. Normalizacja jest kluczowa, aby móc porównywać ekspresję genów między próbkami.
Estymacja dyspersji (Estimation of dispersion): Ocenia zmienność ekspresji dla każdego genu. Dane RNA-Seq charakteryzują się tym, że wariancja zliczeń dla danego genu jest często większa niż jego średnia ekspresja (nadmierna dyspersja). DESeq2 modeluje tę zależność, estymując dyspersję dla każdego genu, a następnie "kurcząc" (shrinking) te oszacowania w kierunku średniej zależności dyspersji od ekspresji dla wszystkich genów. Poprawia to stabilność oszacowań, szczególnie dla genów o niskiej ekspresji.
Testowanie statystyczne (Fitting and Testing): Dla każdego genu dopasowuje uogólniony model liniowy (GLM) z rozkładem ujemnym dwumianowym, zgodnie z formułą design. Następnie przeprowadza testy statystyczne (domyślnie test Walda) w celu oceny istotności współczynników modelu, które odpowiadają log2 krotności zmiany (log2 Fold Change) między porównywanymi grupami.
```r
# Uruchomienie pełnej analizy DESeq2 (normalizacja, estymacja dyspersji, testowanie GLM)
dds <- DESeq(dds)
```
4.4. Pobranie Wyników
Po zakończeniu obliczeń przez DESeq(), wyniki analizy różnicowej ekspresji można wyodrębnić za pomocą funkcji results(). Domyślnie zwraca ona wyniki dla porównania ostatniego poziomu faktora w formule design względem poziomu referencyjnego (w naszym przypadku: "treated" vs "untreated").
Tabela wyników zawiera kluczowe informacje dla każdego genu:
* baseMean: Średnia znormalizowana ekspresja genu we wszystkich próbkach.
* log2FoldChange: Logarytm o podstawie 2 krotności zmiany ekspresji między porównywanymi grupami (np. treated / untreated). Wartość dodatnia oznacza wzrost ekspresji w grupie "treated", ujemna - spadek.
* lfcSE: Błąd standardowy oszacowania log2FoldChange.
* stat: Wartość statystyki testowej (Wald statistic).
* pvalue: Nominalna wartość p z testu Walda, oceniająca istotność zaobserwowanej zmiany log2FoldChange.
* padj: Wartość p skorygowana o problem testowania wielokrotnego (ponieważ testujemy tysiące genów jednocześnie). Domyślnie używana jest metoda Benjamini-Hochberg (BH), która kontroluje odsetek fałszywych odkryć (False Discovery Rate - FDR). Jest to główna wartość używana do oceny istotności statystycznej.
```r
# Pobranie tabeli wyników dla domyślnego porównania (treated vs untreated)
res <- results(dds)

# Wyświetlenie podsumowania wyników (liczba genów up/down-regulowanych przy różnych progach padj)
summary(res)

# Wyświetlenie pierwszych kilku wierszy tabeli wyników
head(res)
```
## 5. Filtrowanie Wyników i Identyfikacja DEG

Po uzyskaniu tabeli wyników z `DESeq2` (obiekt `res`), kolejnym krokiem jest jej przetworzenie w celu zidentyfikowania genów, które wykazują *istotną statystycznie* i *biologicznie znaczącą* zmianę ekspresji.

### 5.1. Wstępne Filtrowanie i Konwersja

Najpierw konwertujemy obiekt wyników `DESeqResults` na standardową ramkę danych (`data.frame`) w R, co ułatwia dalsze manipulacje. Dodajemy również kolumnę z identyfikatorami genów, które domyślnie są przechowywane jako nazwy wierszy.

Następnie usuwamy geny, dla których `DESeq2` nie było w stanie obliczyć skorygowanej wartości p (`padj`). Dzieje się tak zazwyczaj w przypadku genów o bardzo niskiej liczbie zliczeń lub o ekstremalnie wysokiej dyspersji (outliers), dla których testy statystyczne nie są wiarygodne. Filtrowanie to usuwa "szum" i pozwala skupić się na genach, dla których mamy pełne wyniki statystyczne.

```r
# Konwersja obiektu DESeqResults na data.frame
res_tab <- as.data.frame(res)

# Dodanie kolumny z identyfikatorami genów (z nazw wierszy)
res_tab$gene_id <- rownames(res_tab)

# Filtrowanie: usunięcie wierszy (genów) z brakującą wartością padj (NA)
# is.na(res_tab$padj) zwraca TRUE dla wierszy z NA
# !is.na(res_tab$padj) zwraca TRUE dla wierszy BEZ NA
res_sig <- res_tab[!is.na(res_tab$padj), ]

# Sprawdzenie liczby genów po filtrowaniu NA
# nrow(res_sig)
```
### 5.2. Definiowanie Progów Istotności
Aby uznać gen za różnicowo wyeksponowany (DEG), stosujemy dwa główne kryteria:
Próg istotności statystycznej (padj): Określa, jak bardzo jesteśmy pewni, że zaobserwowana zmiana ekspresji nie jest dziełem przypadku. Używamy skorygowanej wartości p (padj), która kontroluje odsetek fałszywych identyfikacji (FDR) przy testowaniu tysięcy genów jednocześnie. Typowy próg to padj < 0.05, co oznacza, że akceptujemy maksymalnie 5% fałszywych identyfikacji wśród zidentyfikowanych DEG.
Próg biologicznej istotności zmiany (log2FoldChange): Określa minimalną wielkość zmiany ekspresji, którą uznajemy za biologicznnie znaczącą. Samo uzyskanie niskiej wartości padj nie oznacza, że zmiana jest duża – może być statystycznie istotna, ale bardzo niewielka. Próg dla log2FoldChange pozwala odfiltrować geny o małych, choć statystycznie istotnych, zmianach. Często stosowany próg to abs(log2FoldChange) > 1, co odpowiada co najmniej 2-krotnej zmianie ekspresji (wzrostowi lub spadkowi). Wartość progu może zależeć od kontekstu biologicznego eksperymentu.

### 5.3. Identyfikacja Genów DEG
Łącząc oba kryteria, wybieramy z przefiltrowanej tabeli res_sig te geny, które spełniają oba warunki jednocześnie.
```r
# Definicja progów
log2FC_threshold <- 1.0 # Odpowiada 2-krotnej zmianie (2^1 = 2)
padj_threshold <- 0.05   # Standardowy próg FDR

# Wybór genów spełniających oba kryteria
# abs(res_sig$log2FoldChange) > log2FC_threshold : Bezwzględna wartość zmiany > próg
# res_sig$padj < padj_threshold : Skorygowana wartość p < próg
DEG <- res_sig[abs(res_sig$log2FoldChange) > log2FC_threshold & res_sig$padj < padj_threshold, ]

# Wyświetlenie liczby zidentyfikowanych genów DEG
cat("Liczba zidentyfikowanych genów DEG:", nrow(DEG), "\n")

# Wyświetlenie pierwszych kilku wierszy tabeli DEG
head(DEG)
```
## 6. Wizualizacja Wyników

Wizualizacja danych i wyników jest kluczowym elementem analizy RNA-Seq. Pozwala na ocenę jakości danych, zrozumienie globalnych wzorców ekspresji oraz graficzną prezentację wyników analizy różnicowej ekspresji. Poniżej przedstawiono kilka standardowych typów wykresów używanych w tym kontekście.

### 6.1. MA Plot

MA Plot jest standardową wizualizacją wyników analizy różnicowej ekspresji. Przedstawia on:
* Na osi Y: Logarytm o podstawie 2 krotności zmiany (M = log2FoldChange).
* Na osi X: Średnią znormalizowaną ekspresję genu we wszystkich próbkach (A = baseMean), zazwyczaj w skali logarytmicznej.

Każdy punkt na wykresie reprezentuje jeden gen. Wykres pozwala ocenić zależność między wielkością zmiany ekspresji a średnim poziomem ekspresji genu. Geny uznane za istotnie różnicowo wyeksponowane (na podstawie progu `padj`) są zazwyczaj wyróżnione innym kolorem (domyślnie czerwonym w funkcji `plotMA`). Oczekujemy, że większość genów (szare punkty) będzie skupiona wokół linii M=0 (brak zmiany), podczas gdy geny DEG będą rozrzucone powyżej i poniżej tej linii. Ograniczenie osi Y (`ylim`) może pomóc w lepszej wizualizacji większości punktów.

```r
# Generowanie MA Plot
# res: obiekt wyników z results()
# ylim: ograniczenie zakresu osi Y dla lepszej czytelności
plotMA(res, ylim=c(-2,2))
title("MA Plot") # Dodanie tytułu
```
### 6.2. Dispersion Plot
Ten wykres jest ważnym narzędziem diagnostycznym do oceny, jak DESeq2 poradziło sobie z modelowaniem dyspersji (zmienności) danych zliczeń. Pokazuje on:
* Czarne punkty: Oszacowania dyspersji dla każdego genu (gene-wise estimates).
* Czerwoną linię: Dopasowana zależność dyspersji od średniej znormalizowanej ekspresji.
* Niebieskie punkty (lub okręgi wokół czarnych): Ostateczne oszacowania dyspersji użyte w modelu GLM, po "skurczeniu" (shrinkage) w kierunku dopasowanej czerwonej linii. Shrinkage poprawia stabilność oszacowań, szczególnie dla genów o niskiej ekspresji.


Oczekujemy, że czarne punkty będą generalnie rozłożone wokół czerwonej linii, a niebieskie punkty będą bliżej niej niż czarne. Odstające punkty (wysoko ponad czerwoną linią) mogą wskazywać na geny o nietypowo dużej zmienności.
```r
# Generowanie Dispersion Plot
plotDispEsts(dds)
title("Dispersion Estimates Plot") # Dodanie tytułu
```
### 6.3. Transformacja Danych (VST)
Wiele standardowych metod wizualizacji i analizy wielowymiarowej (jak PCA czy klastrowanie) działa lepiej na danych homoskedastycznych, czyli takich, gdzie wariancja nie zależy silnie od średniej. Surowe dane zliczeń RNA-Seq są heteroskedastyczne (wariancja rośnie wraz ze średnią). Aby ustabilizować wariancję, DESeq2 oferuje dwie główne transformacje:

* ***VST*** (Variance Stabilizing Transformation): Szybsza, szczególnie dla większych zbiorów danych.
* ***rlog*** (Regularized Logarithm Transformation): Może działać lepiej dla małych zbiorów danych (mała liczba próbek), ale jest wolniejsza obliczeniowo.

Tutaj używamy VST. Argument blind=FALSE oznacza, że transformacja uwzględnia informacje z formuły design, co jest zalecane, gdy celem jest wizualizacja różnic między grupami zdefiniowanymi w projekcie. blind=TRUE jest używane, gdy chcemy ocenić podobieństwo próbek bez uwzględniania projektu (np. do kontroli jakości).
```r
# Przeprowadzenie transformacji stabilizującej wariancję (VST)
# blind=FALSE: transformacja uwzględnia projekt eksperymentu (design)
vsd <- vst(dds, blind=FALSE)

# Alternatywnie można użyć rlog:
# rld <- rlog(dds, blind=FALSE)
```
Wynikiem jest obiekt DESeqTransform, którego macierz przetransformowanych wartości można uzyskać za pomocą funkcji assay().
### 6.4. PCA Plot (Principal Component Analysis)
Analiza głównych składowych (PCA) jest techniką redukcji wymiarowości, która pozwala zwizualizować dane wielowymiarowe (ekspresja tysięcy genów) w przestrzeni o niższej wymiarowości (zazwyczaj 2D). PCA identyfikuje główne osie zmienności w danych. PCA Plot pokazuje, jak próbki (punkty) układają się względem siebie w tej zredukowanej przestrzeni.

Oczekujemy, że próbki pochodzące z tej samej grupy eksperymentalnej (np. ten sam condition lub type) będą grupować się blisko siebie, a grupy różniące się warunkami będą odseparowane. Jest to ważny sposób na wizualną ocenę, czy największa zmienność w danych jest związana z badanymi warunkami eksperymentalnymi oraz czy nie ma wyraźnych próbek odstających (outliers). Wykres generowany jest na danych przetransformowanych (np. VST).
```r
# Generowanie PCA Plot na danych po transformacji VST
# intgroup: zmienne z colData używane do kolorowania/kształtowania punktów
pca_plot <- plotPCA(vsd, intgroup=c("condition", "type"))

# Dodanie tytułu i poprawa estetyki wykresu (opcjonalnie, wymaga ggplot2)
pca_plot <- pca_plot + ggtitle("PCA Plot (VST Data)") + theme_bw()

# Wyświetlenie wykresu
print(pca_plot)
```
### 6.5. Sample Distance Heatmap
Innym sposobem oceny podobieństwa między próbkami jest obliczenie macierzy odległości (np. Euklidesowych) między nimi na podstawie przetransformowanych danych ekspresji (np. VST). Macierz ta jest następnie wizualizowana jako mapa ciepła (heatmap).

Kolory na mapie reprezentują odległość między parami próbek – im ciemniejszy kolor, tym mniejsza odległość (większe podobieństwo). Próbki są zazwyczaj grupowane (klastrowane) hierarchicznie na podstawie ich odległości. Oczekujemy, że próbki z tej samej grupy eksperymentalnej utworzą klastry (bloki o ciemnym kolorze na diagonali tych bloków).
```r
# Obliczenie macierzy odległości Euklidesowych między próbkami
# t(assay(vsd)): transpozycja macierzy VST, aby próbki były w wierszach
sampleDists <- dist(t(assay(vsd)))

# Konwersja obiektu 'dist' na macierz
sampleDistMatrix <- as.matrix(sampleDists)

# Poprawa nazw wierszy/kolumn dla czytelności
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL # Usunięcie nazw kolumn (są takie same jak wierszy)

# Definicja palety kolorów
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# Generowanie heatmapy (wymaga pheatmap)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists, # Użyj obliczonych odległości do klastrowania
         clustering_distance_cols=sampleDists,
         col=colors,
         main="Sample Distance Heatmap (VST Data)") # Tytuł
```
### 6.6. Volcano Plot
Volcano Plot jest jedną z najpopularniejszych metod wizualizacji wyników analizy różnicowej ekspresji. Łączy on informacje o istotności statystycznej i wielkości zmiany ekspresji:
* ***Oś X:*** Log2 krotność zmiany (log2FoldChange).
* ***Oś Y:*** Ujemny logarytm o podstawie 10 skorygowanej wartości p (-log10(padj)). Wyższe wartości na osi Y oznaczają większą istotność statystyczną (niższe padj).

Każdy punkt to gen. Wykres ma kształt przypominający wulkan. Geny o dużej zmianie ekspresji (daleko od zera na osi X) i wysokiej istotności statystycznej (wysoko na osi Y) znajdują się w górnych rogach wykresu ("erupcja wulkanu"). Linie progowe dla padj i log2FoldChange są często dodawane, aby wyraźnie oddzielić geny DEG od genów nieistotnych. Punkty można kolorować w zależności od tego, czy gen jest up-regulowany, down-regulowany, czy nieistotny.
```r
# Przygotowanie danych (dodanie kolumny 'significant' - wymaga dplyr)
# Używamy tabeli res_sig (po usunięciu NA w padj)
res_sig <- res_sig %>%
  mutate(significant = ifelse(padj < padj_threshold & abs(log2FoldChange) > log2FC_threshold,
                                     ifelse(log2FoldChange > log2FC_threshold, "Up-regulated", "Down-regulated"),
                                     "Not Significant"))

# Tworzenie Volcano Plot (wymaga ggplot2)
volcano_plot <- ggplot(res_sig, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) + # Punkty
  scale_color_manual(values = c("Up-regulated" = "red", "Down-regulated" = "blue", "Not Significant" = "grey")) + # Kolory
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold), linetype = "dashed", color = "grey40") + # Linie progowe pionowe
  geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "grey40") + # Linia progowa pozioma
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value",
       color = "Significance") + # Etykiety
  theme_bw() + # Styl
  theme(legend.position = "bottom") # Legenda na dole

# Wyświetlenie wykresu
print(volcano_plot)
```
### 6.7. Counts Plot for Single Gene
Czasami chcemy zobaczyć, jak wygląda rozkład znormalizowanych zliczeń dla konkretnego, interesującego nas genu (np. genu o największej zmianie ekspresji lub znanego markera) we wszystkich próbkach. Funkcja plotCounts z DESeq2 ułatwia przygotowanie danych do takiego wykresu. Zwraca ona znormalizowane zliczenia dla wybranego genu, wraz z informacjami z colData (np. condition, type). Następnie można użyć ggplot2 do stworzenia wykresu punktowego (lub słupkowego, skrzynkowego itp.), pokazującego poziomy ekspresji w poszczególnych próbkach, pogrupowane według warunków. Skala logarytmiczna na osi Y jest często pomocna dla danych zliczeń.
```r
# Sprawdzenie, czy znaleziono jakiekolwiek geny DEG
if (nrow(DEG) > 0) {
  # Wybór genu o najwyższej up-regulacji
  top_gene <- DEG[order(-DEG$log2FoldChange), "gene_id"][1]

  # Pobranie znormalizowanych zliczeń dla wybranego genu
  # returnData=TRUE: zwraca ramkę danych zamiast generować wykres
  gene_counts <- plotCounts(dds, gene = top_gene, intgroup = c("condition", "type"), returnData = TRUE)

  # Tworzenie wykresu punktowego (wymaga ggplot2)
  counts_plot <- ggplot(gene_counts, aes(x = condition, y = count, color = type)) +
    geom_point(position = position_jitter(width = 0.1), size = 3) + # Punkty z lekkim rozrzutem
    scale_y_log10() + # Skala log10 na osi Y
    labs(title = paste("Normalized Counts for Gene:", top_gene),
         x = "Condition",
         y = "Normalized Counts (log10 scale)",
         color = "Library Type") + # Etykiety
    theme_bw() # Styl

  # Wyświetlenie wykresu
  print(counts_plot)
} else {
  print("Nie znaleziono genów DEG, pominięto wykres zliczeń.")
}
```
## 7. Analiza Funkcjonalna i strukturalna (GO i KEGG) z gProfiler2

Po zidentyfikowaniu listy genów różnicowo wyeksponowanych (DEG), naturalnym kolejnym krokiem jest próba zrozumienia, jakie procesy biologiczne, funkcje molekularne lub szlaki metaboliczne są nadreprezentowane w tej grupie genów. Do tego celu służy analiza funkcjonalna (Functional Analysis). Jednym z popularnych narzędzi do jej przeprowadzenia jest `gProfiler` (dostępny w R jako pakiet `gprofiler2`).

Metoda ta, w najprostszej formie znana jako Over-Representation Analysis (ORA), porównuje naszą listę genów DEG z predefiniowanymi zestawami genów (np. geny związane z konkretnym terminem Gene Ontology lub ścieżką KEGG). Statystycznie ocenia, czy dany zestaw genów (termin GO, ścieżka KEGG) występuje w naszej liście DEG częściej, niż można by się tego spodziewać losowo, biorąc pod uwagę całe tło (np. wszystkie geny analizowane w eksperymencie).

### 7.1. Przygotowanie Danych i Uruchomienie Analizy

Najpierw sprawdzamy, czy w ogóle zidentyfikowaliśmy jakieś geny DEG. Jeśli tak, wyodrębniamy ich identyfikatory. Następnie używamy funkcji `gost()` z pakietu `gprofiler2`.

Kluczowe argumenty funkcji `gost()`:
* `query`: Wektor identyfikatorów genów z naszej listy zainteresowania (tutaj: `deg_gene_ids`).
* `organism`: Nazwa organizmu, dla którego przeprowadzamy analizę. Musi być zgodna z typem identyfikatorów genów w `query` i bazami danych `gProfiler`. Dla danych Pasilla (Drosophila melanogaster, identyfikatory FlyBase) używamy `"dmelanogaster"`.
* `sources`: Wektor określający, które bazy danych (źródła adnotacji) mają być przeszukiwane. Typowe wybory to:
    * `"GO:BP"`: Gene Ontology - Biological Process (procesy biologiczne)
    * `"GO:MF"`: Gene Ontology - Molecular Function (funkcje molekularne)
    * `"GO:CC"`: Gene Ontology - Cellular Component (komponenty komórkowe)
    * `"KEGG"`: Kyoto Encyclopedia of Genes and Genomes (ścieżki metaboliczne, sygnałowe itp.)
    * Inne bazy, np. `"REAC"` (Reactome), `"WP"` (WikiPathways), `"HP"` (Human Phenotype Ontology).
* `evcodes`: `TRUE` dołącza kody dowodowe dla adnotacji GO (informują o typie dowodu eksperymentalnego/obliczeniowego dla danej adnotacji).
* `user_threshold`: Próg istotności (skorygowanej wartości p) dla wyników analizy funkcjonalnej (np. 0.05).
* `correction_method`: Metoda korekcji dla testowania wielokrotnego (domyślna i zalecana `"g_SCS"`).

```r
# Sprawdzenie, czy istnieją geny DEG do analizy
if (nrow(DEG) > 0) {
  # Wyodrębnienie identyfikatorów genów DEG
  deg_gene_ids <- DEG$gene_id

  # Uruchomienie analizy funkcjonalnej za pomocą gProfiler2
  gost_results <- gost(query = deg_gene_ids,
                       organism = "dmelanogaster", # Ważne: dopasuj do organizmu i typu ID
                       sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG"), # Wybrane bazy danych
                       evcodes = TRUE, # Dołącz kody dowodowe GO
                       user_threshold = 0.05, # Próg istotności padj
                       correction_method = "g_SCS") # Metoda korekcji

  # Sprawdzenie, czy uzyskano wyniki
  if (!is.null(gost_results) && nrow(gost_results$result) > 0) {
    # Wyodrębnienie tabeli wyników
    functional_analysis_results <- gost_results$result # Zmieniono nazwę zmiennej dla jasności

    # Wyświetlenie pierwszych kilku wierszy istotnych wyników analizy funkcjonalnej
    # Sortowanie wg skorygowanej wartości p (p_value) i wybór kluczowych kolumn
    print(head(functional_analysis_results[order(functional_analysis_results$p_value),
                           c("source", "term_id", "term_name", "p_value", "intersection_size")]))

    # Opcjonalnie: Wizualizacja wyników za pomocą gostplot
    # Tworzy wykres słupkowy (lub inny) pokazujący termy/ścieżki istotne w analizie funkcjonalnej
    # interactive = FALSE generuje statyczny wykres ggplot
    gost_plot_static <- gostplot(gost_results, capped = TRUE, interactive = FALSE)
    print(gost_plot_static)

  } else {
    print("Nie znaleziono istotnych wyników analizy funkcjonalnej dla podanej listy DEG i progów.")
  }

} else {
  print("Pominięto analizę funkcjonalną, ponieważ nie znaleziono genów DEG.")
}
```
### 7.2. Interpretacja Wyników
Tabela wyników (functional_analysis_results) zawiera listę terminów GO i ścieżek KEGG, które są istotnie nadreprezentowane w naszej liście DEG. Kluczowe kolumny to:
* source: Baza danych, z której pochodzi termin (np. GO:BP, KEGG).
* term_id: Unikalny identyfikator terminu/ścieżki.
* term_name: Nazwa terminu/ścieżki.
* p_value: Skorygowana wartość p dla testu analizy funkcjonalnej. Niższe wartości wskazują na większą istotność statystyczną.
* intersection_size: Liczba genów z naszej listy DEG, które należą do danego terminu/ścieżki.

Analiza tych wyników pozwala sformułować hipotezy dotyczące pełnionych funkcji zaobserwowanych zmian ekspresji genów. Na przykład, jeśli wiele istotnych terminów GO:BP dotyczy odpowiedzi immunologicznej, sugeruje to, że knockdown genu pasilla wpływa na procesy odpornościowe w komórkach. Wykres gostplot dostarcza wizualnego podsumowania najbardziej istotnych wyników analizy funkcjonalnej. 
